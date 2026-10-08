# A21 pre-registration: PP3 x curated PM1 on functional truth, expanded to MaveDB, with gene as unit and a negative control

Written and committed before any R or credit was computed on MaveDB data (8 October 2026). Already done: the list of
MaveDB score sets per gene (work/mavedb/candidates.tsv: 16 genes with curated PM1 residues have score sets), curated PM1
residues (841 in 67 genes, as in drafts 15-18), and the v16 ProteinGym result (7 assays, 3 genes, R 0.49).

## Data and definitions
- Curated PM1: protein positions where any ClinGen expert panel applied PM1 to a missense variant (Evidence Repository,
  this export). AlphaMissense: protein-level scores (AlphaMissense_aa_substitutions), UniProt canonical isoform.
- Score sets: human MaveDB score sets for these genes with protein-level variants (hgvs_pro). Positions are mapped to the
  UniProt sequence: the offset (0..2000) maximising reference-residue agreement is used; a score set is kept if >= 90% of
  its variants' reference residues match. Curated residues are checked the same way (>= 90% reference agreement).
- Truth (damaging): orientation and threshold from the score set itself. Requires >= 10 nonsense (Ter) and >= 10
  synonymous (=) variants; damaging = missense score closer to the nonsense median than to the synonymous median
  (threshold at their midpoint). Score sets without enough nonsense/synonymous variants are excluded.
- Informative score set: >= 10 assayed missense at curated residues, >= 10 damaging and >= 10 not damaging.
- PP3: AlphaMissense >= 0.792 (Supporting) or >= 0.972 (Strong) (Bergquist et al. 2025 thresholds as used in drafts 15-18).

## Analysis
- Per score set: R = LR(PP3 AND F)/(LR(PP3) LR(F)) and PP3's share of log-evidence kept at curated residues
  (ln LR(PP3 | F) / ln LR(PP3)), position-clustered bootstrap (500).
- Gene is the unit: per gene, score sets are combined by inverse-variance mean of ln R; across genes, DerSimonian-Laird
  random effects; leave-one-gene-out.
- Primary: random-effects R across genes for AlphaMissense Supporting x curated PM1. Dependence confirmed if the upper
  95% bound < 1.
- Negative control (pre-specified): in each score set, a random residue set of the same size as the curated set, drawn
  from assayed residues not curated (200 draws, seed 20261008); R(PP3, random flag) summarised by its median across draws.
  The control is passed if the random-effects control R interval includes 1 or lies above the curated estimate.
- Secondary: Strong band; share kept; ProteinGym-only vs MaveDB-only.

## Failure
If the primary interval includes 1, or the negative control shows the same R as curated PM1, we report that functional
truth does not confirm the dependence beyond the three genes of draft 16.
