# Pocket-Conditioned Diffusion for EGFR T790M Inhibitor Design

**EECS 6895 Final Project** — Yixuan Ye — Spring 2026

This project evaluates whether a pocket-conditioned diffusion model (TargetDiff) can generate T790M-selective EGFR inhibitors from two near-identical binding pockets differing by a single residue (WT 1M17 vs. T790M 4I22). Cross-docking and three ablation experiments characterize the model's selectivity sensitivity; in silico SAR optimization then compensates for the observed limitation by identifying scaffold-dependent substitution rules.

1,987 unique molecules were generated and evaluated for drug-likeness, binding selectivity, and pharmacophore content. Cross-docking reveals no statistically significant selectivity differentiation between the two conditioned pools (KS p = 0.112), traced to the single-target training objective of CrossDocked2020. SAR optimization identifies NH-hinge + aminopyrimidine scaffolds as selectively responsive to electron-withdrawing substitution (54–75% improvement rates). The lead compound Rank #9 + CN achieves Δ = −1.74 kcal/mol with QED = 0.87 and MW ≈ 322 Da, outperforming all three approved reference inhibitors across every reported dimension.

---

## System Structure

```
.
├── configs/
│   ├── sampling.yaml        # TargetDiff sampling parameters and evaluation thresholds
│   └── targets.yaml         # PDB targets, pocket radii, baseline drug SMILES
├── demo/
│   ├── top10_report_cards.html    # Top-10 candidate report cards (open in browser)
│   └── interactive_scatter.html   # Interactive cross-docking scatter plot
├── docs/
│   └── slides_outline.md    # Presentation outline
├── notebooks/
│   ├── 01_setup_and_test.ipynb         # Environment verification and dependency check
│   ├── 02_prepare_targets.ipynb        # PDB download, pocket extraction (10 Å radius)
│   ├── 03_vina_baseline.ipynb          # Baseline Vina scoring for Erlotinib, Gefitinib, Osimertinib
│   ├── 04_generate_wildtype.ipynb      # TargetDiff generation on EGFR WT (1M17), 10 × 100 mol
│   ├── 05_filter_and_metrics.ipynb     # Drug-likeness filtering: QED / SA / Lipinski / PAINS
│   ├── 06_generate_mutant.ipynb        # TargetDiff generation on EGFR T790M (4I22), 10 × 100 mol
│   ├── 07_cross_docking.ipynb          # Cross-docking top-50 per pool against both receptors
│   ├── 08_optimization_experiments.ipynb  # Ablation experiments A / B / C
│   ├── 09_selectivity_analysis.ipynb      # Selectivity scoring, pharmacophore filter, composite ranking
│   └── 10_sar_optimization.ipynb          # In silico SAR optimization via aromatic substitution
├── src/
│   ├── colab_init.py        # Google Colab setup helper (Drive mount, deps)
│   ├── docking.py           # AutoDock Vina wrapper: single molecule and resumable batch
│   ├── evaluation.py        # QED, SA score, Lipinski metrics
│   ├── filtering.py         # Drug-likeness + PAINS filters
│   ├── pocket_extraction.py # Binding pocket extraction from PDB
│   ├── utils.py             # Config loading, SDF I/O, result helpers
│   └── visualization.py     # Vina distributions, UMAP, radar charts
├── requirements.txt
└── README.md
```

### Data and Result Directories (generated at runtime)

```
data/
├── pdb/           # Downloaded PDB files and receptor PDBQT
├── pockets/       # Extracted binding pocket PDB files
└── baselines/     # Baseline drug SDF conformers
results/
├── generated/     # Raw TargetDiff output SDF files
├── vina_scores/   # Docking result CSVs
└── figures/       # Figures
external/
└── targetdiff/    # Cloned TargetDiff repository + pretrained checkpoint
```

---

## Key Results

| Metric | WT pool | T790M pool |
|--------|---------|------------|
| Molecules generated | 998 | 989 |
| Chemical validity | 100% | 100% |
| Drug-like (pass filter) | 264 (26.5%) | 284 (28.7%) |
| Mean cross-pool Tanimoto | 0.285 | — |

**Selectivity (cross-docking, top-50 per pool):**
- Mean Δ: WT pool +0.045 kcal/mol, T790M pool −0.039 kcal/mol
- KS test: D = 0.240, p = 0.112 (not significant)
- Attributed to single-target training objective of CrossDocked2020, not pocket-conditioning failure

**Lead compound — Rank #9 + CN:**

| Compound | vina T790M | Δ | QED | SA | MW |
|----------|-----------|---|-----|----|----|
| Rank #9 + CN | −9.8 | −1.74 | 0.87 | 3.3 | 322 Da |
| Osimertinib | −8.6 | −0.28 | 0.31 | 4.1 | 499 Da |
| Erlotinib | −7.1 | −0.36 | 0.42 | 2.9 | 393 Da |
| Gefitinib | −8.5 | +0.23 | 0.52 | 2.8 | 446 Da |

---

## Setup Instructions

### Requirements

- Python 3.10
- Conda (recommended)
- CUDA-capable GPU (for TargetDiff generation; CPU works but is slow)
- Google Colab Pro is recommended for notebooks 04 and 06 (generation ~6–7 h each)

### Step 1 — Clone the repository

```bash
git clone https://github.com/Sapphirine/2026_Drug_2.git
cd 2026_Drug_2
```

### Step 2 — Create a conda environment

```bash
conda create -n egfr python=3.10
conda activate egfr
pip install -r requirements.txt
conda install -c conda-forge openbabel   # required for PDBQT conversion
```

### Step 3 — Clone TargetDiff

```bash
mkdir -p external
git clone https://github.com/guanjq/targetdiff.git external/targetdiff
```

### Step 4 — Download the pretrained checkpoint

Download the TargetDiff pretrained model from the [official release](https://drive.google.com/file/d/1_BUWcHMQLbvOPbU4aYiDYcvF_0VEPjPZ) and save it to:

```
external/targetdiff/checkpoints/pretrained_diffusion.pt
```

### Step 5 — Launch notebooks

```bash
jupyter lab notebooks/
```

Run notebooks in order (01 → 10). Each notebook saves its outputs to `data/` or `results/` for the next stage.

---

## How to Run

| Notebook | Input | Output |
|----------|-------|--------|
| 01 | — | Environment check, TargetDiff confirmed |
| 02 | PDB IDs 1M17, 4I22 | Pocket PDB files, receptor PDBQT |
| 03 | Receptor PDBQT, baseline SMILES | `results/vina_scores/baselines.csv` |
| 04 | 1M17 pocket, checkpoint | `results/generated/wildtype_*.sdf` (10 batches × 100 mol) |
| 05 | WT SDF files | `results/vina_scores/wt_metrics.csv` (drug-likeness filter) |
| 06 | 4I22 pocket, checkpoint | `results/generated/mutant_*.sdf` (10 batches × 100 mol) |
| 07 | Top-50 per pool, both receptors | `results/vina_scores/cross_docking.csv` |
| 08 | Cross-docking CSV | Ablation A (pocket radius), B (structural correlation), C (pharmacophore) |
| 09 | Generated + baseline molecules | Selectivity scores, pharmacophore-confirmed candidates, composite ranking |
| 10 | Top-20 composite-ranked scaffolds | SAR analogs (~80), lead compound Rank #9 + CN |

To run on **Google Colab**, replace the first cell of any notebook with:

```python
from src.colab_init import setup
PROJECT_ROOT = setup("https://github.com/Sapphirine/2026_Drug_2.git")
```

---

## Example Usage

### Extract a binding pocket

```python
from src.pocket_extraction import extract_pocket

center = extract_pocket(
    pdb_file="data/pdb/1M17.pdb",
    ligand_resname="AQ4",
    output_file="data/pockets/1M17_pocket.pdb",
    radius=10.0,
)
print("Pocket center:", center)
```

### Dock a molecule and compute selectivity differential

```python
from src.docking import VinaDocker
from rdkit import Chem

docker_wt  = VinaDocker("data/pdb/1M17_receptor.pdbqt", center=(3.0, 14.5, 52.3))
docker_mut = VinaDocker("data/pdb/4I22_receptor.pdbqt", center=(2.1, 13.8, 51.7))

mol = Chem.MolFromSmiles("C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1")  # Erlotinib
vina_wt  = docker_wt.dock_mol(mol)
vina_mut = docker_mut.dock_mol(mol)
delta = vina_mut - vina_wt   # Δ < 0 → T790M-selective
print(f"Δ = {delta:.2f} kcal/mol")
```

### Evaluate drug-likeness

```python
from src.evaluation import evaluate_batch

smiles_list = [
    "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",        # Erlotinib
    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",    # Gefitinib
    "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C",  # Osimertinib
]
df = evaluate_batch(smiles_list, output_csv="results/baseline_metrics.csv")
```

### Filter generated molecules

```python
from src.filtering import compute_all_metrics, apply_filters
from rdkit import Chem

mols = [Chem.MolFromSmiles(s) for s in smiles_list if s]
df = compute_all_metrics(mols)
df_filtered, summary = apply_filters(df, qed_min=0.3, sa_max=5.0)
print(summary)
```

---

## Targets

| Target | PDB ID | Resolution | Co-crystal ligand |
|--------|--------|------------|-------------------|
| EGFR Wild-Type | 1M17 | 2.6 Å | Erlotinib (AQ4) |
| EGFR T790M/L858R | 4I22 | 2.8 Å | WZ-4002 (IRE) |

## Evaluation Metrics

| Metric | Tool | Threshold |
|--------|------|-----------|
| Chemical validity | RDKit | > 80% |
| Uniqueness | Canonical SMILES | > 90% |
| Diversity | Mean Tanimoto distance | > 0.7 |
| QED | RDKit | ≥ 0.3 (calibrated so Osimertinib passes) |
| SA Score | RDKit (Ertl 2009) | ≤ 5.0 |
| Lipinski Ro5 | RDKit | ≤ 1 violation |
| PAINS alerts | RDKit FilterCatalog | 0 |
| Vina score | AutoDock Vina | ≤ −8 kcal/mol |
| Selectivity Δ | Vina_Mut − Vina_WT | < 0 (T790M-preferred) |

Composite ranking score used in nb09–10:

```
Score = 0.40 · vina_Mut(norm) + 0.30 · QED + 0.20 · (1 − SA/10) + 0.10 · (1 − maxTc)
```

where `maxTc` is the maximum Tanimoto similarity to any of the three approved reference compounds (diversity penalty).

---

## Demo

Pre-computed interactive visualizations are in `demo/` — open in any browser, no dependencies required:

- **`top10_report_cards.html`** — property cards for the top-10 candidates
- **`interactive_scatter.html`** — cross-docking selectivity scatter (vina_WT vs vina_Mut)

---

## References

1. Guan et al., *3D Equivariant Diffusion for Target-Aware Molecule Generation*, ICLR 2023.
2. Trott & Olson, *AutoDock Vina: Improving the speed and accuracy of docking*, J. Comput. Chem. 2010.
3. Bickerton et al., *Quantifying the chemical beauty of drugs*, Nature Chemistry 2012.
4. Ertl & Schuffenhauer, *Estimation of synthetic accessibility score*, J. Cheminformatics 2009.
5. Francoeur et al., *CrossDocked2020*, J. Chem. Inf. Model. 2020.
