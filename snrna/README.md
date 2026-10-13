# snRNA-VEP v2

A selection-supervised, ClinGen-calibrated variant-effect framework for spliceosomal snRNA genes. This folder is separate from the PP3/PM1 study at the repository root.

| File | What it is |
|---|---|
| `PLAN.md` | The study plan (13 October 2026) |
| `PREREG.md` | The frozen, executable pre-registration. It also holds the deviation log (section 14) and the proposed amendments (section 15) |
| `prereg_manifest.txt` + `.ots` | SHA-256 of `PREREG.md` and `PLAN.md`, OpenTimestamps-stamped before any endpoint was computed |
| `code_manifest.txt` + `.ots` | Same for the analysis code, stamped before the week-0 and G1 runs |
| `snrna.py` | The pipeline (all commands below) |
| `evo2_score.py` | Evo 2 scoring for gate G1 (needs a 24 GB GPU) |
| `evo2_run/` | Evo 2 input windows (committed) and, once run, the scores |
| `g2_power.py` | Simulated power and false-positive rate of the G2 tests (no real data) |
| `model_features.py` | Builds `model_inputs/`: per-gene graphs (snRNA + partner RNAs), Rfam profiles, ΔΔG, paralog clusters |
| `evo2_features.py` | Evo 2 embeddings and all-SNV Δ log-likelihoods for the model (GPU) |
| `model.py`, `train.py` | The selection model and every pre-registered fit (amendment A2) |
| `evaluate.py` | Aim 2, ablations, negative controls, G3, exploratory endpoints |
| `labels.py` | Non-blind patient labels (`labels/`) and sealing of the blind genes' labels |
| `run_a100.sh` | One command for the GPU machine: seal → Evo 2 features → train → evaluate → pack |
| `tests/` | Known-answer self-tests: `python -m pytest -q tests` (pipeline, model on synthetic counts, mock Evo 2, sealing) |
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
| Week 1: **G1 decision** | **Done: PASSES, continue** | Evo 2 7B (scored on the authors' A100): ρ 0.303 [0.22, 0.38], AUROC 0.670, against the bar of ρ ≥ 0.60 and AUROC ≥ 0.85. Best label-free tool, but far from enough; T-loop ρ 0.000 |
| G2 test design | Amendment A1 adopted before any G2 data | Stronger blind-hotspot test (simulation-based) |
| **G2 (selection signal)** | **Done: PARTIAL, continue** | Mutation model M0 frozen on 42 training genes before testing. RNU4-2 critical region O/E 0.13 (p 9×10⁻⁹): recovered. Blind RNU2-2 n.1–67 O/E 0.62 (p 2×10⁻²⁶): recovered. Blind RNU5B-1 n.37–44 O/E 0.16 (p 2×10⁻⁸) missed on circular p 0.051 vs 0.05. Post hoc: every placement that tied it overlaps the hotspot (`results/g2_report.txt`) |
| Selection model (amendment A2) | **Code complete; tested on synthetic counts and a mock Evo 2. Awaiting the A100 run** | Patient labels for non-blind genes tallied (`labels/L2_nonblind.tsv`: 77 labelled SNVs, P/LP and B/LB, all mapped to model variants). Blind labels are sealed on the GPU machine as stage 1 |

### G1 so far (RNU4-2 SGE, 521 in-transcript variants)

| Predictor | ρ with function score [95% CI] | AUROC, depleted vs not | ReNU vs gnomAD controls (12 vs 200) |
|---|---|---|---|
| CADD v1.7 | 0.145 [0.061, 0.230] | 0.598 [0.547, 0.652] | 0.621 [0.529, 0.711] |
| ViennaRNA U4/U6 ΔΔG | 0.269 [0.187, 0.350] | 0.656 [0.602, 0.707] | 0.719 [0.564, 0.858] |
| Allele-count floor (gnomAD) | 0.204 [0.122, 0.284] | 0.629 [0.588, 0.668] | (undefined) |
| Position-only floor (leave-one-position-out on SGE) | 0.570 [0.508, 0.627] | 0.867 [0.832, 0.898] | 0.922 (circular, see below) |
| Evo 2 7B | 0.303 [0.220, 0.381] | 0.670 [0.619, 0.720] | 0.658 [0.552, 0.749] |

- **Reproduction check passes.** On the paper's 346 population SNVs, CADD gives 0.649, ΔΔG 0.718 and SGE 0.934, against the published 0.65, 0.72 and 0.93.
- **The position-only floor is the bar to beat.** A label-free model has to beat ρ 0.57, not the 0.27 of the best published baseline. The floor is fitted on SGE scores and uses the critical region, so it is a supervised reference. Its ReNU AUROC is circular, because the critical region is defined from the ReNU positions.
- **The frozen G1 floor rule does not cover this case.** It ends the deep-model framing only if the floor is *within* 0.05 of the best predictor; here the floor is 0.30 *above* every predictor. This is the authors' call (deviation log, 13 October 2026).

### G2 in one paragraph

Healthy-population data (gnomAD v4.1) show strong depletion exactly where dominant disease variants cluster, and none where recessive ones do. In RNU4-2, the T-loop and Stem III have O/E 0.10 and 0.22, while the recessive Stem II, k-turn and Sm site have 0.80–1.40. The RNU4ATAC recessive regions sit at 0.83–1.09. Both blind dominant or mixed regions are depleted: RNU2-2 n.1–67 at 0.62 and RNU5B-1 n.37–44 at 0.16. By the pre-registered rule, G2 is PARTIAL: RNU5B-1 misses the circular-shift criterion by 0.001. A post-hoc check shows the gene's most depleted stretch is n.33–44, overlapping the hotspot, so the miss reflects a limitation of the shift test for short hotspots, not an absent signal. That explanation is post hoc and does not change the verdict.

## How to run

```bash
pip install numpy scipy pandas openpyxl scikit-learn pysam ViennaRNA==2.7.0 pytest
python -m pytest -q tests        # self-tests, no network
python snrna.py fetch            # SGE table, GENCODE v49, GIAB (about 100 MB)
python snrna.py genes            # coverage + paralog QC (about 6 min; gnomAD API, cached)
python snrna.py sites            # all possible SNVs + Roulette/CADD/gnomAD (remote tabix, cached)
python snrna.py sge              # SGE known-fact checks; stops on any mismatch
python snrna.py vienna           # U4/U6 ddG QC
python snrna.py evo2-prep        # 8,192-nt windows for Evo 2 -> evo2_run/ (already committed)
python snrna.py g1               # G1 report: results/g1_report.txt
INFERNAL_BIN=/path/to/infernal/bin python snrna.py structure   # Rfam paired/unpaired labels (Infernal 1.1.5)
python snrna.py g2-fit           # choose + freeze the mutation model on training genes; commit its output
python snrna.py g2-test          # one-shot G2 tests: results/g2_report.txt
```

Remote files are read by tabix over HTTPS and cached under `data/cache/`, so re-runs are offline.

### Finishing G1 on an RTX 4090 (Windows + WSL2 or Linux)

The 521 input windows are already committed (`evo2_run/evo2_inputs.tsv.gz`), so only the scoring step runs on your machine.

1. **Check the GPU.** Update the NVIDIA Windows driver; `nvidia-smi` inside WSL should list the RTX 4090. Keep about 25 GB of disk free (the weights are about 15 GB).
2. **Get this branch in a fresh folder** (your PP3 working copy is untouched):
   ```bash
   git clone -b claude/snrna-vep-v2-framework-b0abng https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration.git ~/snrna-run
   cd ~/snrna-run/snrna
   ```
3. **Set up the environment.** This is Evo 2's "light install" for the 7B model; it needs no Transformer Engine:
   ```bash
   conda create -n evo2 python=3.12 -y && conda activate evo2
   pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
   pip install packaging ninja psutil pandas
   MAX_JOBS=4 pip install flash-attn==2.8.0.post2 --no-build-isolation   # if it compiles for long or fails, send the error
   pip install evo2
   python -m evo2.test.test_evo2_generation --model_name evo2_7b      # downloads the weights, checks the install
   ```
4. **Smoke test, then the full run.** The script resumes, so the second command only scores the remaining 516:
   ```bash
   python evo2_score.py evo2_run/evo2_inputs.tsv.gz evo2_run/evo2_scores.tsv --limit 5
   python evo2_score.py evo2_run/evo2_inputs.tsv.gz evo2_run/evo2_scores.tsv
   ```
   The smoke test prints seconds per variant and minutes left; expect very roughly 30–90 minutes in total. Close other GPU-heavy apps. On an out-of-memory error, stop and report it; the 8,192-nt window is frozen (PREREG section 6.1).
5. **Push the scores.** Then `python snrna.py g1` (or Claude) applies the frozen G1 rule. `g1` ignores an incomplete scores file:
   ```bash
   git add evo2_run/evo2_scores.tsv evo2_run/evo2_scores.tsv.json
   git commit -m "Evo 2 7B scores for G1" && git push
   ```

## The selection model on the A100 (amendment A2)

Everything runs from one script. It skips finished work, so after any interruption you re-run the same command.

```bash
cd ~/snrna-run && git pull && cd snrna
conda activate evo2
bash run_a100.sh --dry-run          # 20-30 min: checks everything on the GPU; fits no real count
```

Send back the last 40 lines of the dry run (or upload `snrna_dryrun_*.tar.gz`). The dry run:

- runs the fast tests;
- runs real Evo 2 on one gene and prints its speed;
- times the largest model configuration and prints an estimate of the full run;
- trains and evaluates on **synthetic** counts, so its numbers mean nothing.

Then start the full run inside tmux, so it survives a closed browser tab:

```bash
tmux new -s snrna 'bash run_a100.sh; bash'
# detach: Ctrl-b then d.   Re-attach later: tmux attach -t snrna
```

1. **Seal (asks once, about 1 min).** Choose a passphrase and keep it; the analyst must not know it. The script fetches the RNU2-2 and RNU5B-1 patient labels straight into `sealed/blind_labels.tsv.enc` and prints only counts and SHA-256 values. Nothing is ever shown in clear text. To run fully unattended, `export SNRNA_SEAL_PASS=...` first.
2. **Evo 2 features (about 5–10 h).** All 17,385 SNVs on both strands, plus per-position embeddings. The script prints minutes left.
3. **Training (about 3–8 h).** It runs, in order:
   - the hyperparameter grid (48 fits);
   - strict and loose fits, and the 14 leave-cluster-out fits;
   - 12 ablations;
   - 3 negative controls.

   Each scenario uses 5 seeds. Every fit is saved as it finishes.
4. **Evaluate and pack (about 10 min).** Upload the resulting `snrna_a100_results_*.tar.gz` (a few MB) to the chat, or paste `results/model_report.txt`.

If anything stops, the log in `logs/` and the last lines on screen say where. Fix it and re-run the same command.

## After the A100 run

Upload the results tarball to the chat. The analysis then runs in this order; each step only reads earlier outputs:

```bash
python evaluate.py                     # Aim 2, ablations, controls, G3 (done on the A100 as part of run_a100.sh)
python comparators_report.py           # GPN-MSA, phyloP, ReMM, Evo 2 on the same endpoints
python calibrate.py                    # preliminary Aim 4 on the non-blind genes
python release.py --stage scores       # dated Zenodo bundle of all scores: upload it, record the DOI in PREREG s14
SNRNA_SEAL_PASS=... python final_eval.py --reason "final evaluation, PREREG s9"   # on the machine with the sealed file
python release.py --stage final        # scores + calibrated PP3/BP4 thresholds and evidence labels
```

`final_eval.py` unseals the blind labels exactly once (logged) and refuses to run on a dry run.

## Author actions

1. **Deposit on OSF.** `PREREG.md`, `PLAN.md` and both manifests with their `.ots` files. The OpenTimestamps proofs already predate the analyses; the OSF copy will be later and should be described as such.
2. **Upgrade the timestamps.** After a few hours, run `ots upgrade prereg_manifest.txt.ots code_manifest.txt.ots` and commit the result, so the proofs are Bitcoin-anchored and verifiable offline.
3. **Run Evo 2.** On the RTX 4090, as above, to decide G1.
4. **Decide on amendment A1** (`PREREG.md` section 15) before G2 is run. The frozen hotspot test has about 57% power even for a 50-fold depleted 8-nt hotspot.
5. **Note the RNU2-2 region type.** Guidance Table 2 was verified from the Europe PMC full text (`table2_regions.tsv`); both blind hotspots match `PLAN.md`. But the RNU2-2 region is a recessive region, so its population depletion may be modest. This weighs on decision 4.
6. **Seal the blind labels.** `run_a100.sh` does this as stage 1 (amendment A2.10), before any model sees real counts, and `train.py` refuses to start otherwise. Keep the passphrase; it is needed once, at the final evaluation. RNUdb serves RNU2-2 variants at `/api/genes/RNU2-2/variants`; the analyst must not open it. Literature-only patient variants for non-blind genes can be added later in `labels/L2_curated_additions.tsv` (same columns, dated).
7. **Send the outreach note.** It is in `PLAN.md`, and should go with the G1 result.
