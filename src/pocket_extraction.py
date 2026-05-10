"""Extract binding pocket residues from PDB around a reference ligand."""

import numpy as np
from pathlib import Path
from typing import Optional
from Bio import PDB


def get_ligand_center(pdb_file: str, ligand_resname: str) -> np.ndarray:
    """Return the geometric center of the named ligand in the PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_resname:
                    for atom in residue:
                        coords.append(atom.get_vector().get_array())

    if not coords:
        raise ValueError(f"Ligand '{ligand_resname}' not found in {pdb_file}.")
    return np.mean(coords, axis=0)


def extract_pocket(
    pdb_file: str,
    ligand_resname: str,
    output_file: str,
    radius: float = 10.0,
    include_ligand: bool = True,
) -> np.ndarray:
    """
    Extract protein residues within `radius` Å of the ligand center.

    Returns the ligand center coordinates (used for Vina box definition).
    """
    center = get_ligand_center(pdb_file, ligand_resname)

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    class PocketSelect(PDB.Select):
        def accept_residue(self, residue):
            for atom in residue:
                if np.linalg.norm(atom.get_vector().get_array() - center) <= radius:
                    return 1
            if include_ligand and residue.get_resname().strip() == ligand_resname:
                return 1
            return 0

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file, PocketSelect())

    print(f"Pocket saved → {output_file}  (center: {center.round(2)})")
    return center


def prepare_receptor_pdbqt(
    pdb_file: str,
    output_pdbqt: str,
    mgltools_prepare: Optional[str] = None,
) -> str:
    """Convert receptor PDB to PDBQT for AutoDock Vina.

    Strips heteroatoms (waters, co-crystal ligands, ions) before conversion
    to avoid spurious Vina contacts in the binding pocket.
    """
    import subprocess

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    class KeepProtein(PDB.Select):
        def accept_residue(self, residue):
            return residue.id[0] == " "  # blank hetflag = standard amino acid

    clean_pdb = str(pdb_file).replace(".pdb", "_clean.pdb")
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(clean_pdb, KeepProtein())

    if mgltools_prepare:
        cmd = [
            "python", mgltools_prepare,
            "-r", clean_pdb,
            "-o", output_pdbqt,
            "-A", "hydrogens",
            "-U", "nphs_lps",
        ]
    else:
        cmd = ["obabel", clean_pdb, "-O", output_pdbqt, "-xr",
               "--partialcharge", "gasteiger"]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Receptor preparation failed:\n{result.stderr}")
    print(f"Receptor PDBQT saved → {output_pdbqt}")
    return output_pdbqt
