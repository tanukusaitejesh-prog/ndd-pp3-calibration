# snRNA-VEP v2

A selection-supervised, ClinGen-calibrated variant-effect framework for spliceosomal snRNA genes. This folder is separate from the PP3/PM1 study at the repository root.

| File | What it is |
|---|---|
| `PLAN.md` | The study plan (13 October 2026) |
| `PREREG.md` | The frozen, executable pre-registration. It also holds the deviation log (section 14) and the proposed amendments (section 15) |
| `prereg_manifest.txt` + `.ots` | SHA-256 of `PREREG.md` and `PLAN.md`, OpenTimestamps-stamped before any endpoint was computed |
| `code_manifest.txt` + `.ots` | Same for the analysis code, stamped before the week-0 and G1 runs |
| `snrna.py` | The pipeline (all commands below) |
| `evo2_score.py` | Evo 2 scoring for gate G1 (needs a GPU) |
| `g2_power.py` | Simulated power and false-positive rate of the G2 tests (no real data) |
| `tests/` | 13 known-answer self-tests: `python -m pytest -q tests` (about 2.5 min) |
| `table2_regions.tsv` | Guidance Table 2 regions, verified from the full text |
| `results/` | All outputs, committed |

## Status (13 October 2026)

| Step | Status | Result |
|---|---|---|
| Pre-registration | **Frozen and timestamped** (OpenTimestamps). OSF deposit not yet done | `PREREG.md` |
| Week 0: gene universe | **Done** | 111 genes, which reproduces the guidance; 87 pass coverage, also reproducing the guidance; 55 after paralog QC. All target genes are retained |
| Week 0: sites | **Done** | 17,385 SNVs, 0 REF conflicts, Roulette rates for 98.2% |
| Week 0: SGE checks | **All pass** | 539/151/521/12/346; replicate r 0.835–0.865; CADD column equals CADD v1.7; ΔΔG reproduced exactly |
| Week 1: G1 floors and baselines | **Done** | See below |
| Week 1: **G1 decision** | **Pending** | Evo 2 needs a GPU (none in this environment) |
| G2 test design | Simulated | The frozen blind-hotspot test is underpowered; amendment A1 is proposed, not adopted |

### G1 so far (RNU4-2 SGE, 521 in-transcript variants)

| Predictor | ρ with function score [95% CI] | AUROC, depleted vs not | ReNU vs gnomAD controls (12 vs 200) |
|---|---|---|---|
| CADD v1.7 | 0.145 [0.061, 0.230] | 0.598 [0.547, 0.652] | 0.621 [0.529, 0.711] |
| ViennaRNA U4/U6 ΔΔG | 0.269 [0.187, 0.350] | 0.656 [0.602, 0.707] | 0.719 [0.564, 0.858] |
| Allele-count floor (gnomAD) | 0.204 [0.122, 0.284] | 0.629 [0.588, 0.668] | (undefined) |
| Position-only floor (leave-one-position-out on SGE) | 0.570 [0.508, 0.627] | 0.867 [0.832, 0.898] | 0.922 (circular, see below) |
| Evo 2 7B | pending | pending | pending |

- **Reproduction check passes.** On the paper's 346 population SNVs, CADD gives 0.649, ΔΔG 0.718 and SGE 0.934, against the published 0.65, 0.72 and 0.93.
- **The position-only floor is the bar to beat.** A label-free model has to beat ρ 0.57, not the 0.27 of the best published baseline. The floor is fitted on SGE scores and uses the critical region, so it is a supervised reference. Its ReNU AUROC is circular, because the critical region is defined from the ReNU positions.
- **The frozen G1 floor rule does not cover this case.** It ends the deep-model framing only if the floor is *within* 0.05 of the best predictor; here the floor is 0.30 *above* every predictor. This is the authors' call (deviation log, 13 October 2026).

## How to run

```bash
pip install numpy scipy pandas openpyxl scikit-learn pysam ViennaRNA==2.7.0 pytest
python -m pytest -q tests        # self-tests, no network
python snrna.py fetch            # SGE table, GENCODE v49, GIAB (about 100 MB)
python snrna.py genes            # coverage + paralog QC (about 6 min; gnomAD API, cached)
python snrna.py sites            # all possible SNVs + Roulette/CADD/gnomAD (remote tabix, cached)
python snrna.py sge              # SGE known-fact checks; stops on any mismatch
python snrna.py vienna           # U4/U6 ddG QC
python snrna.py evo2-prep        # 8,192-nt windows for Evo 2
python snrna.py g1               # G1 report: results/g1_report.txt
```

Remote files are read by tabix over HTTPS and cached under `data/cache/`, so re-runs are offline.

### Finishing G1 (needs a GPU with at least 24 GB)

```bash
python snrna.py evo2-prep
python evo2_score.py work/evo2_inputs.tsv.gz work/evo2_scores.tsv   # 521 variants x 4 forward passes of 8,192 nt
python snrna.py g1                                                   # applies the frozen G1 rule
```

## Author actions

1. **Deposit on OSF.** `PREREG.md`, `PLAN.md` and both manifests with their `.ots` files. The OpenTimestamps proofs already predate the analyses; the OSF copy will be later and should be described as such.
2. **Upgrade the timestamps.** After a few hours, run `ots upgrade prereg_manifest.txt.ots code_manifest.txt.ots` and commit the result, so the proofs are Bitcoin-anchored and verifiable offline.
3. **Run Evo 2.** On a GPU, as above, to decide G1.
4. **Decide on amendment A1** (`PREREG.md` section 15) before G2 is run. The frozen hotspot test has about 57% power even for a 50-fold depleted 8-nt hotspot.
5. **Note the RNU2-2 region type.** Guidance Table 2 was verified from the Europe PMC full text (`table2_regions.tsv`); both blind hotspots match `PLAN.md`. But the RNU2-2 region is a recessive region, so its population depletion may be modest. This weighs on decision 4.
6. **Tally and seal the blind labels.** Tally patient variants, and have someone other than the analyst compile `sealed/blind_labels.tsv` for RNU2-2 and RNU5B-1. Seal it with `python snrna.py seal sealed/blind_labels.tsv`, then record both SHA-256 values in a dated amendment before any selection model is trained. RNUdb now serves RNU2-2 variants at `/api/genes/RNU2-2/variants`; the analyst must not open it.
7. **Send the outreach note.** It is in `PLAN.md`, and should go with the G1 result.
