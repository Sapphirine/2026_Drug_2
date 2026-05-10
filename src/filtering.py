"""Drug-likeness filters: PAINS, Lipinski violations, Tanimoto similarity to reference set."""

from __future__ import annotations
from typing import Iterable
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from src.evaluation import compute_sa_score, lipinski_details



def get_pains_catalog() -> FilterCatalog:
    """Return RDKit PAINS filter catalog (A + B + C subsets)."""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


def count_pains_alerts(mol: Chem.Mol, catalog: FilterCatalog | None = None) -> int:
    """Number of PAINS substructure matches; 0 means clean."""
    if mol is None:
        return -1
    if catalog is None:
        catalog = get_pains_catalog()
    return len(catalog.GetMatches(mol))



def lipinski_violations(mol: Chem.Mol) -> int:
    """Count violations of the four Ro5 rules (0 = perfect, 4 = all violated)."""
    mw = Descriptors.ExactMolWt(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    logp = Descriptors.MolLogP(mol)
    return int(mw > 500) + int(hbd > 5) + int(hba > 10) + int(logp > 5)



def morgan_fp(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, n_bits)


def tanimoto_max(mol: Chem.Mol, ref_fps: list, ref_names: list[str] | None = None
                 ) -> tuple[float, str | None]:
    """Return (max Tanimoto, name of closest reference) for a single molecule."""
    if mol is None or not ref_fps:
        return 0.0, None
    fp = morgan_fp(mol)
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    idx = max(range(len(sims)), key=lambda i: sims[i])
    name = ref_names[idx] if ref_names else None
    return float(sims[idx]), name



def compute_all_metrics(
    mols: Iterable[Chem.Mol],
    baseline_fps: list | None = None,
    baseline_names: list[str] | None = None,
) -> pd.DataFrame:
    """Compute physicochemical properties and filter flags for a list of molecules.

    If baseline_fps is provided, also adds MaxTanimotoToBaseline and ClosestBaseline.
    Skips None entries; returns only successfully parsed molecules.
    """
    catalog = get_pains_catalog()
    records = []
    for i, mol in enumerate(mols):
        if mol is None:
            continue
        rec: dict = {"mol_id": i, "smiles": Chem.MolToSmiles(mol)}
        rec.update(lipinski_details(mol))
        rec["QED"] = round(QED.qed(mol), 4)
        rec["SA"] = compute_sa_score(mol)
        rec["LipinskiViolations"] = lipinski_violations(mol)
        rec["PAINS_Alerts"] = count_pains_alerts(mol, catalog)
        if baseline_fps is not None:
            sim, name = tanimoto_max(mol, baseline_fps, baseline_names)
            rec["MaxTanimotoToBaseline"] = round(sim, 3)
            rec["ClosestBaseline"] = name
        records.append(rec)
    return pd.DataFrame(records)



def apply_filters(
    df: pd.DataFrame,
    qed_min: float = 0.3,
    sa_max: float = 5.0,
    lipinski_max_violations: int = 1,
    pains_max: int = 0,
    mw_range: tuple[float, float] = (200, 600),
    logp_range: tuple[float, float] = (-2, 6),
) -> tuple[pd.DataFrame, dict]:
    """Apply drug-likeness hard filters; return (kept_df, per-stage counts dict).

    QED floor is 0.3 (not the typical 0.4) so that covalent drugs like Osimertinib
    (QED=0.31 due to its acrylamide warhead) are not discarded.
    """
    summary = {"input": len(df)}
    mask = pd.Series(True, index=df.index)

    rules = [
        (df["QED"]                  >= qed_min,                f"QED >= {qed_min}"),
        (df["SA"]                   <= sa_max,                 f"SA <= {sa_max}"),
        (df["LipinskiViolations"]   <= lipinski_max_violations, f"Lipinski violations <= {lipinski_max_violations}"),
        (df["PAINS_Alerts"]         <= pains_max,              f"PAINS alerts <= {pains_max}"),
        (df["MW"].between(*mw_range),                          f"MW in {mw_range}"),
        (df["LogP"].between(*logp_range),                      f"LogP in {logp_range}"),
    ]
    for rule_mask, label in rules:
        summary[label] = int(rule_mask.sum())
        mask &= rule_mask

    summary["all combined"] = int(mask.sum())
    return df[mask].copy().reset_index(drop=True), summary



def load_sdf_mols(sdf_path: str) -> list[Chem.Mol]:
    """Read an SDF, return only successfully sanitized molecules."""
    return [m for m in Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=False)
            if m is not None]
