"""AutoDock Vina docking wrapper: single-molecule and batch scoring."""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import Optional
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm



def cross_dock_smiles_list(
    smiles_list: list[str],
    docker_wt,
    docker_mut: "VinaDocker",
    output_csv: str,
    source_label: str = "",
    checkpoint_every: int = 25,
) -> pd.DataFrame:
    """Dock each SMILES against both receptors (WT + Mutant) and save score pairs.

    Resumable: if `output_csv` exists, mol_ids already present are skipped.

    Output columns: mol_id, source, smiles, vina_WT, vina_Mut, status_WT, status_Mut
    """
    # Load existing progress if any
    done_ids: set[int] = set()
    existing = None
    if os.path.exists(output_csv):
        existing = pd.read_csv(output_csv)
        done_ids = set(existing["mol_id"].tolist())
        print(f"[resume] {len(done_ids)} molecules already in {output_csv}, skipping them")

    new_records = []
    pending = [(i, s) for i, s in enumerate(smiles_list) if i not in done_ids]
    if not pending:
        print(f"[resume] all {len(smiles_list)} molecules already docked.")
        return pd.read_csv(output_csv) if existing is not None else pd.DataFrame()

    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    pbar = tqdm(pending, desc=f"{source_label} cross-docking")
    for i, smi in pbar:
        mol = Chem.MolFromSmiles(smi)
        rec = {"mol_id": i, "source": source_label, "smiles": smi,
               "vina_WT": None, "vina_Mut": None,
               "status_WT": "parse_error", "status_Mut": "parse_error"}
        if mol is not None:
            score_wt = docker_wt.dock_mol(mol)
            rec["vina_WT"] = score_wt
            rec["status_WT"] = "ok" if score_wt is not None else "dock_failed"

            score_mut = docker_mut.dock_mol(mol)
            rec["vina_Mut"] = score_mut
            rec["status_Mut"] = "ok" if score_mut is not None else "dock_failed"

        new_records.append(rec)
        pbar.set_postfix({"vina_WT": rec["vina_WT"], "vina_Mut": rec["vina_Mut"]})

        # Checkpoint append
        if len(new_records) % checkpoint_every == 0:
            _flush_records(output_csv, existing, new_records)

    # Final flush
    _flush_records(output_csv, existing, new_records)
    out = pd.read_csv(output_csv)
    n_ok = ((out["status_WT"] == "ok") & (out["status_Mut"] == "ok")).sum()
    print(f"\n[done] {n_ok}/{len(out)} molecules have both WT + Mut scores -> {output_csv}")
    return out


def _flush_records(output_csv: str, existing: Optional[pd.DataFrame],
                   new_records: list[dict]) -> None:
    """Append the in-memory new_records to the output CSV (overwrite-merge)."""
    new_df = pd.DataFrame(new_records)
    if existing is not None and len(existing):
        merged = pd.concat([existing, new_df], ignore_index=True)
    else:
        merged = new_df
    merged = merged.drop_duplicates(subset=["mol_id", "source"], keep="last")
    merged.to_csv(output_csv, index=False)


class VinaDocker:
    def __init__(
        self,
        receptor_pdbqt: str,
        center: tuple[float, float, float],
        box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
        exhaustiveness: int = 8,
        num_modes: int = 9,
        vina_bin: str = "vina",
        cpu: int | None = None,   # Vina --cpu; None = auto-detect all cores
    ):
        self.receptor = receptor_pdbqt
        self.center = center
        self.box_size = box_size
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.vina_bin = vina_bin
        self.cpu = cpu if cpu is not None else (os.cpu_count() or 1)

    def _mol_to_pdbqt(self, mol: Chem.Mol, tmp_dir: str) -> Optional[str]:
        """Convert RDKit Mol to PDBQT via Open Babel."""
        sdf_path = os.path.join(tmp_dir, "ligand.sdf")
        pdbqt_path = os.path.join(tmp_dir, "ligand.pdbqt")

        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) < 0:
            return None
        AllChem.MMFFOptimizeMolecule(mol)

        writer = Chem.SDWriter(sdf_path)
        writer.write(mol)
        writer.close()

        result = subprocess.run(
            ["obabel", sdf_path, "-O", pdbqt_path, "--partialcharge", "gasteiger"],
            capture_output=True,
        )
        if result.returncode != 0 or not os.path.exists(pdbqt_path):
            return None
        return pdbqt_path

    def dock_mol(self, mol: Chem.Mol) -> Optional[float]:
        """Dock a single RDKit Mol; return best Vina score or None on failure."""
        with tempfile.TemporaryDirectory() as tmp:
            ligand_pdbqt = self._mol_to_pdbqt(mol, tmp)
            if ligand_pdbqt is None:
                return None

            out_path = os.path.join(tmp, "out.pdbqt")
            cmd = [
                self.vina_bin,
                "--receptor", self.receptor,
                "--ligand", ligand_pdbqt,
                "--center_x", str(self.center[0]),
                "--center_y", str(self.center[1]),
                "--center_z", str(self.center[2]),
                "--size_x", str(self.box_size[0]),
                "--size_y", str(self.box_size[1]),
                "--size_z", str(self.box_size[2]),
                "--exhaustiveness", str(self.exhaustiveness),
                "--num_modes", str(self.num_modes),
                "--cpu", str(self.cpu),
                "--out", out_path,
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                return None

            for line in result.stdout.splitlines():
                if line.strip().startswith("1 "):
                    try:
                        return float(line.split()[1])
                    except (IndexError, ValueError):
                        pass
        return None


def dock_molecules(
    sdf_file: str,
    docker: VinaDocker,
    output_csv: str,
) -> pd.DataFrame:
    """Dock all molecules in an SDF file; checkpoint every 50 molecules."""
    supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
    records = []

    for i, mol in enumerate(tqdm(supplier, desc="Docking")):
        if mol is None:
            records.append({"idx": i, "smiles": None, "vina_score": None, "status": "parse_error"})
            continue

        smiles = Chem.MolToSmiles(mol)
        score = docker.dock_mol(mol)
        records.append({
            "idx": i,
            "smiles": smiles,
            "vina_score": score,
            "status": "ok" if score is not None else "docking_failed",
        })

        # Checkpoint every 50 molecules
        if (i + 1) % 50 == 0:
            pd.DataFrame(records).to_csv(output_csv + ".tmp", index=False)  # crash-recovery checkpoint

    df = pd.DataFrame(records)
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    print(f"Docking results saved → {output_csv}  ({df['vina_score'].notna().sum()}/{len(df)} succeeded)")
    return df
