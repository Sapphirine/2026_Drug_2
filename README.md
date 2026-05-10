# EGFR Inhibitor Design via Target-Conditioned 3D Molecular Generation

**EECS 6895 Final Project** — Spring 2026

This project builds a computational drug discovery pipeline for EGFR inhibitors. Starting from experimental crystal structures, the pipeline extracts binding pockets, generates novel 3D molecules using the TargetDiff diffusion model, scores them with AutoDock Vina, and applies multi-criteria drug-likeness filters to identify top candidates.

---

## System Structure

```
.
├── configs/
│   ├── sampling.yaml        # TargetDiff sampling + evaluation thresholds
│   └── targets.yaml         # PDB targets, pocket radii, baseline drug SMILES
├── docs/
│   └── slides_outline.md    # Presentation outline
├── notebooks/
│   ├── 01_setup_and_test.ipynb       # Environment verification and dependency check
│   ├── 02_prepare_targets.ipynb      # PDB download, pocket extraction
│   ├── 03_vina_baseline.ipynb        # AutoDock Vina baseline scoring
│   ├── 04_generate_wildtype.ipynb    # TargetDiff generation on EGFR WT (1M17)
│   ├── 05_filter_and_metrics.ipynb   # Validity, QED, SA, Lipinski filtering
│   ├── 06_generate_mutant.ipynb      # TargetDiff generation on EGFR T790M (4I22)
│   ├── 07_cross_docking.ipynb        # Cross-docking WT vs T790M candidates
│   ├── 08_optimization_experiments.ipynb  # Optimization experiments A/B/C
│   ├── 09_selectivity_analysis.ipynb      # Selectivity scoring and pharmacophore filter
│   └── 10_sar_optimization.ipynb          # In silico SAR optimization
├── src/
│   ├── colab_init.py        # Google Colab setup helper (Drive mount, deps)
│   ├── docking.py           # AutoDock Vina wrapper (single mol + batch)
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
└── figures/       # Publication-quality figures
external/
└── targetdiff/    # Cloned TargetDiff repository + checkpoint
```

---

## Setup Instructions

### Requirements

- Python 3.10
- Conda (recommended)
- CUDA-capable GPU (for TargetDiff generation; CPU works but is slow)
- Google Colab Pro is recommended for notebooks 04 and 06 (generation takes ~6–7 h each)

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

Run notebooks in order (01 → 10). Each notebook mounts the previous notebook's outputs.

---

## How to Run

The pipeline is designed to run as a sequence of notebooks. Each notebook is self-contained and saves its outputs to `data/` or `results/` for the next stage.

| Notebook | Input | Output |
|----------|-------|--------|
| 01 | — | Environment check, TargetDiff clone confirmed |
| 02 | PDB IDs (1M17, 4I22) | Pocket PDB files, receptor PDBQT |
| 03 | Receptor PDBQT, baseline SMILES | `results/vina_scores/baselines.csv` |
| 04 | 1M17 pocket, TargetDiff checkpoint | `results/generated/wildtype_*.sdf` |
| 05 | WT SDF files | `results/vina_scores/wt_metrics.csv` |
| 06 | 4I22 pocket, TargetDiff checkpoint | `results/generated/mutant_*.sdf` |
| 07 | WT + T790M SDF, both receptors | `results/vina_scores/cross_docking.csv` |
| 08 | Cross-docking CSV | Optimization experiment results |
| 09 | Generated + baseline molecules | Selectivity scores, pharmacophore-filtered candidates |
| 10 | Selectivity CSV | SAR-optimized Top candidates |

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

### Dock a molecule

```python
from src.docking import VinaDocker

docker = VinaDocker(
    receptor_pdbqt="data/pdb/1M17_receptor.pdbqt",
    center=(3.0, 14.5, 52.3),   # from pocket extraction
    exhaustiveness=8,
)
from rdkit import Chem
mol = Chem.MolFromSmiles("C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1")
score = docker.dock_mol(mol)
print(f"Vina score: {score:.2f} kcal/mol")
```

### Evaluate drug-likeness

```python
from src.evaluation import evaluate_batch

smiles_list = [
    "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1",   # Erlotinib
    "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",  # Gefitinib
]
df = evaluate_batch(smiles_list, output_csv="results/baseline_metrics.csv")
```

### Filter generated molecules

```python
from src.filtering import compute_all_metrics, apply_filters
from rdkit import Chem

mols = [Chem.MolFromSmiles(s) for s in smiles_list if s]
df = compute_all_metrics(mols)
df_filtered, summary = apply_filters(df)
print(summary)
```

---

## Targets

| Target | PDB ID | Resolution | Co-crystal ligand |
|--------|--------|------------|-------------------|
| EGFR Wild-Type | 1M17 | 2.6 Å | Erlotinib (AQ4) |
| EGFR T790M mutant | 4I22 | 2.8 Å | IRE |

## Evaluation Metrics

| Metric | Tool | Threshold |
|--------|------|-----------|
| Validity | RDKit | > 80% |
| Uniqueness | Canonical SMILES | > 90% |
| Diversity | Mean Tanimoto distance | > 0.7 |
| Vina score | AutoDock Vina | ≤ −8 kcal/mol |
| QED | RDKit | ≥ 0.5 |
| SA Score | RDKit (Ertl 2009) | ≤ 4.0 |
| Lipinski Ro5 | RDKit | Pass |

---

## References

1. Guan et al., *3D Equivariant Diffusion for Target-Aware Molecule Generation*, ICLR 2023.
2. Trott & Olson, *AutoDock Vina: Improving the speed and accuracy of docking*, J. Comput. Chem. 2010.
3. Bickerton et al., *Quantifying the chemical beauty of drugs*, Nature Chemistry 2012.
4. Ertl & Schuffenhauer, *Estimation of synthetic accessibility score*, J. Cheminformatics 2009.
