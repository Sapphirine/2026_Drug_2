"""Shared helpers: config loading, SDF I/O, results management."""

from __future__ import annotations
import yaml
import pandas as pd
from pathlib import Path
from rdkit import Chem


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def sdf_to_smiles_list(sdf_file: str) -> list[str]:
    """Read an SDF file and return a list of canonical SMILES (None for failures)."""
    supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
    smiles = []
    for mol in supplier:
        if mol is None:
            smiles.append(None)
        else:
            smiles.append(Chem.MolToSmiles(mol))
    return smiles


def smiles_to_sdf(smiles_list: list[str], output_sdf: str):
    """Write a list of SMILES to an SDF file with 3D coordinates."""
    Path(output_sdf).parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(output_sdf)
    from rdkit.Chem import AllChem
    for smi in smiles_list:
        if smi is None:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) >= 0:
            AllChem.MMFFOptimizeMolecule(mol)
            writer.write(mol)
    writer.close()
    print(f"SDF written → {output_sdf}")


def save_results(df: pd.DataFrame, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"Results saved → {path}  ({len(df)} rows)")


def merge_vina_and_metrics(vina_csv: str, metrics_csv: str) -> pd.DataFrame:
    """Join Vina scores with drug-likeness metrics on SMILES column."""
    vina = pd.read_csv(vina_csv)
    metrics = pd.read_csv(metrics_csv)
    merged = vina.merge(metrics, on="smiles", how="inner")
    return merged


def filter_top_candidates(
    df: pd.DataFrame,
    vina_cutoff: float = -8.0,
    qed_cutoff: float = 0.5,
    sa_cutoff: float = 4.0,
    require_lipinski: bool = True,
    top_k: int = 10,
) -> pd.DataFrame:
    """Apply multi-criteria filter and return top-k ranked by Vina score."""
    mask = df["vina_score"] <= vina_cutoff
    if "qed" in df:
        mask &= df["qed"] >= qed_cutoff
    if "sa_score" in df:
        mask &= df["sa_score"] <= sa_cutoff
    if require_lipinski and "Lipinski" in df:
        mask &= df["Lipinski"]

    filtered = df[mask].copy()
    filtered = filtered.sort_values("vina_score").head(top_k)
    print(f"Top-{top_k} candidates: {len(filtered)} pass all filters (from {mask.sum()} qualifying)")
    return filtered
