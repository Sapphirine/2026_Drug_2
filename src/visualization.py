"""Plotting: Vina score distributions, UMAP chemical space, and radar charts."""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import rdDepictor

PALETTE = {"WT": "#4C72B0", "T790M": "#DD8452", "Baseline": "#55A868"}


def _savefig(fig: plt.Figure, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"Figure saved → {path}")
    plt.close(fig)


def plot_score_distribution(
    wt_df: pd.DataFrame,
    mut_df: pd.DataFrame,
    baseline_scores: dict[str, float],
    output_path: str = "results/figures/fig2_vina_distribution.png",
):
    """Vina score distributions for WT vs T790M with baseline drug reference lines."""
    fig, ax = plt.subplots(figsize=(9, 5))

    wt_scores = wt_df["vina_score"].dropna()
    mut_scores = mut_df["vina_score"].dropna()

    ax.hist(wt_scores, bins=40, alpha=0.6, color=PALETTE["WT"], label="EGFR WT (1M17)")
    ax.hist(mut_scores, bins=40, alpha=0.6, color=PALETTE["T790M"], label="EGFR T790M (4I22)")

    colors_b = ["#2ca02c", "#d62728", "#9467bd"]
    for (drug, score), c in zip(baseline_scores.items(), colors_b):
        ax.axvline(score, color=c, linestyle="--", linewidth=1.5, label=f"{drug} ({score:.1f})")

    ax.set_xlabel("Vina Score (kcal/mol)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Generated Molecule Vina Score Distribution vs. Clinical Baselines", fontsize=13)
    ax.legend(fontsize=9)
    ax.invert_xaxis()

    _savefig(fig, output_path)
    return fig


def plot_chemical_space(
    wt_smiles: list[str],
    mut_smiles: list[str],
    baseline_smiles: dict[str, str],
    output_path: str = "results/figures/fig5_umap.png",
):
    """UMAP of Morgan fingerprints for WT, T790M, and baseline drugs."""
    try:
        import umap
    except ImportError:
        print("umap-learn not installed. Run: pip install umap-learn")
        return

    def smiles_to_fp(smiles_list):
        fps = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(list(fp))
        return np.array(fps)

    wt_fp = smiles_to_fp(wt_smiles)
    mut_fp = smiles_to_fp(mut_smiles)
    base_fp = smiles_to_fp(list(baseline_smiles.values()))

    all_fp = np.vstack([wt_fp, mut_fp, base_fp])
    labels = (["WT"] * len(wt_fp) + ["T790M"] * len(mut_fp) + ["Baseline"] * len(base_fp))

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric="jaccard", random_state=42)
    embedding = reducer.fit_transform(all_fp)

    fig, ax = plt.subplots(figsize=(9, 7))
    for label, color in PALETTE.items():
        mask = np.array(labels) == label
        size = 60 if label == "Baseline" else 5
        marker = "*" if label == "Baseline" else "."
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   c=color, s=size, marker=marker, alpha=0.6, label=label)

    base_start = len(wt_fp) + len(mut_fp)
    for i, name in enumerate(baseline_smiles.keys()):
        ax.annotate(name, embedding[base_start + i], fontsize=8, fontweight="bold")

    ax.set_title("Chemical Space (UMAP, Morgan FP)", fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    _savefig(fig, output_path)
    return fig


def plot_radar(
    df_top: pd.DataFrame,
    baseline_df: pd.DataFrame,
    output_path: str = "results/figures/fig3_radar.png",
):
    """Radar chart comparing top-10 generated molecules against baseline drugs."""
    metrics = ["vina_score", "qed", "sa_score", "MW", "LogP"]
    labels = ["Vina\n(neg better)", "QED", "SA Score\n(low better)", "MW", "LogP"]

    def normalize(series, invert=False):
        mn, mx = series.min(), series.max()
        n = (series - mn) / (mx - mn + 1e-8)
        return 1 - n if invert else n

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True))
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]

    for label, df, color in [("Top-10 WT", df_top, PALETTE["WT"]),
                              ("Baselines", baseline_df, PALETTE["Baseline"])]:
        values = []
        for i, m in enumerate(metrics):
            s = df[m].dropna()
            if s.empty:
                values.append(0.5)
                continue
            mean_val = s.mean()
            invert = m in ("vina_score", "sa_score", "MW")
            values.append(float(normalize(pd.Series([mean_val, s.min(), s.max()]), invert)[0]))
        values += values[:1]
        ax.plot(angles, values, "o-", linewidth=2, color=color, label=label)
        ax.fill(angles, values, alpha=0.15, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, size=10)
    ax.set_title("Generated vs. Clinical Drugs (Normalized Metrics)", size=12, pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))
    _savefig(fig, output_path)
    return fig


def render_3d_pose(smiles: str, title: str = "", size: tuple = (400, 300)) -> object:
    """Return a py3Dmol view object for inline Jupyter display."""
    try:
        import py3Dmol
        from rdkit.Chem import AllChem
    except ImportError:
        print("py3Dmol not installed: pip install py3Dmol")
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    mb = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=size[0], height=size[1])
    view.addModel(mb, "mol")
    view.setStyle({"stick": {}})
    view.setBackgroundColor("white")
    view.zoomTo()
    return view
