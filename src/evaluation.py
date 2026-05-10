"""Drug-likeness metrics: Validity, Uniqueness, Diversity, QED, SA Score, Lipinski."""

from __future__ import annotations
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, AllChem
from rdkit.Chem.rdMolDescriptors import CalcTPSA

# SA Score (Ertl & Schuffenhauer 2009) — shipped with RDKit contrib
try:
    from rdkit.Contrib.SA_Score import sascorer
    _SA_AVAILABLE = True
except ImportError:
    _SA_AVAILABLE = False


def compute_sa_score(mol: Chem.Mol) -> float | None:
    if not _SA_AVAILABLE:
        return None
    try:
        return sascorer.calculateScore(mol)
    except Exception:
        return None


def lipinski_pass(mol: Chem.Mol) -> bool:
    """Return True if molecule satisfies all four Lipinski Ro5 rules."""
    mw = Descriptors.ExactMolWt(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    logp = Descriptors.MolLogP(mol)
    return mw <= 500 and hbd <= 5 and hba <= 10 and logp <= 5


def lipinski_details(mol: Chem.Mol) -> dict:
    return {
        "MW": round(Descriptors.ExactMolWt(mol), 2),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(CalcTPSA(mol), 2),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Lipinski": lipinski_pass(mol),
    }


class MolEvaluator:
    """Evaluate a list of SMILES strings and return per-molecule metrics."""

    def evaluate(self, smiles_list: list[str]) -> pd.DataFrame:
        records = []
        for smi in smiles_list:
            row = {"smiles": smi}
            mol = Chem.MolFromSmiles(smi) if smi else None
            row["valid"] = mol is not None
            if mol:
                row["qed"] = round(QED.qed(mol), 4)
                row["sa_score"] = compute_sa_score(mol)
                row.update(lipinski_details(mol))
            else:
                row.update({"qed": None, "sa_score": None,
                            "MW": None, "HBD": None, "HBA": None,
                            "LogP": None, "TPSA": None, "RotBonds": None,
                            "Lipinski": False})
            records.append(row)
        return pd.DataFrame(records)

    @staticmethod
    def validity(df: pd.DataFrame) -> float:
        return df["valid"].mean()

    @staticmethod
    def uniqueness(df: pd.DataFrame) -> float:
        valid_smiles = df.loc[df["valid"], "smiles"].tolist()
        if not valid_smiles:
            return 0.0
        canonical = set()
        for smi in valid_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                canonical.add(Chem.MolToSmiles(mol))
        return len(canonical) / len(valid_smiles)

    @staticmethod
    def diversity(df: pd.DataFrame, sample_n: int = 500) -> float:
        """Mean pairwise Tanimoto distance (1 - similarity) on Morgan FP."""
        valid_smiles = df.loc[df["valid"], "smiles"].tolist()
        if len(valid_smiles) < 2:
            return 0.0

        rng = np.random.default_rng(42)
        if len(valid_smiles) > sample_n:
            valid_smiles = rng.choice(valid_smiles, sample_n, replace=False).tolist()

        fps = []
        for smi in valid_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048))

        distances = []
        for i in range(len(fps)):
            for j in range(i + 1, len(fps)):
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                distances.append(1.0 - sim)

        return float(np.mean(distances)) if distances else 0.0


def evaluate_batch(smiles_list: list[str], output_csv: str | None = None) -> pd.DataFrame:
    """Convenience wrapper: evaluate + print summary."""
    evaluator = MolEvaluator()
    df = evaluator.evaluate(smiles_list)

    print(f"Validity:   {evaluator.validity(df):.1%}")
    print(f"Uniqueness: {evaluator.uniqueness(df):.1%}")
    print(f"Diversity:  {evaluator.diversity(df):.3f} (mean Tanimoto distance)")
    if df["qed"].notna().any():
        print(f"QED:        {df['qed'].mean():.3f} ± {df['qed'].std():.3f}")
    if df["sa_score"].notna().any():
        print(f"SA Score:   {df['sa_score'].mean():.3f} ± {df['sa_score'].std():.3f}")
    if "Lipinski" in df:
        print(f"Lipinski:   {df['Lipinski'].mean():.1%} pass rate")

    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"Metrics saved → {output_csv}")
    return df
