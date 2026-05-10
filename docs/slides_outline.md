# Presentation Slides Outline (12–15 pages)

## Slides 1–2: Problem Definition
- Lung cancer statistics; EGFR mutation prevalence in Asian populations
- Drug resistance: T790M/C797S — unmet clinical need
- Hypothesis: generative AI can explore chemical space beyond existing libraries

## Slide 3–4: Technical Background
- Diffusion models in 3D molecule generation
- SE(3) equivariance: why it matters for molecular geometry
- TargetDiff: E(3)-equivariant diffusion, pocket-conditioned generation

## Slide 5–6: Pipeline
- Figure 1: full pipeline diagram (PDB → pocket → diffusion → Vina → ranking)
- Key design choices: why TargetDiff, why EGFR, why Vina

## Slide 7–10: Results
- Figure 2: Vina score distribution (WT vs T790M vs baselines)
- Figure 5: Chemical space UMAP
- Figure 3: Radar chart — generated Top-10 vs Erlotinib/Gefitinib/Osimertinib
- Figure 4: Top-10 3D binding poses (PyMOL renders)

## Slide 11–12: Discussion & Limitations
- What the results mean scientifically
- Limitations: Vina ≠ real affinity, no MD, no ADMET, training data bias
- What would come next: MM/GBSA, ADMET filters, synthesis

## Slide 13: Summary
- Three takeaways in one sentence each

## Anticipated Q&A
- "Why not train from scratch?" — two-week / single-GPU constraint
- "Is Vina Score predictive?" — no, it's a first-pass filter; cite Vina paper
- "What is SE(3) equivariance?" — rotation/translation invariance baked into architecture
