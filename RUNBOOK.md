# Complete run commands: `ndd-pp3-calibration`, branch `claude/snrna-vep-v2-framework-b0abng`

Every command comes from reading each script in the branch: its argparse options, docstring and runner scripts, plus the
deviation logs that record the order the author actually ran things in. The argparse command-line interfaces (every
snRNA script, `gate.py` and `fullstudy.py`) were checked with `--help`. Two exceptions: `train.py` needs torch, so its
options were read from the source, and the root `make_*.py` scripts take plain positional arguments.

The branch holds **two separate studies**:

- **Part A: snRNA-VEP v2** in `snrna/`. This is what the branch is named for.
- **Part B: the NDD de novo PP3/PM1 study** at the repo root.

## Tags

| Tag | Meaning |
|---|---|
| `[CPU]` | Any Linux machine |
| `[GPU]` | NVIDIA GPU: 24 GB (RTX 4090) for G1 scoring; 80 GB A100/H100 for the Evo 2 features and training |
| `[NET]` | Needs internet (remote tabix, APIs, downloads) |
| `[SEALED]` | Needs the decrypted blind labels `sealed/blind_labels.tsv` (exists only after the unsealing step) |
| `[FROZEN]` | Overwrites a committed or hash-frozen output. Run `git diff` afterwards: any difference means the reproduction failed |

## Run everything in WSL Ubuntu, not Windows

pysam, Infernal, bcftools and Evo 2 don't run on native Windows. Keep the clone in your Linux home (`~/`), not under
`/mnt/c`, which is much slower.

## Files not in git (needed for an *exact* reproduction)

These exist only on the author's DGX or laptop:

- **snRNA study:**
  - `snrna/results_v4/` (the hash-frozen v4 model)
  - `snrna/sealed/` (encrypted blind labels and the unsealing log)
  - `snrna/model_runs/` and `snrna/evo2_run/features/`
  - `snrna/a7/`, plus the A9–A11 fetched tables
  - `snrna/release/`
- **Root study:** the external data folder `/mnt/x/ndd_gate_data` and the column store `work/fullstudy/col`.

Without these files, the commands below *re-run* the study rather than replay it. That matters most for the sealing step: the patient labels are public now, so a fresh sealing is no longer a blind test.

---

# Part A: snRNA-VEP v2 (`snrna/`)

## A0. Setup

```bash
git clone --single-branch -b claude/snrna-vep-v2-framework-b0abng https://github.com/tanukusaitejesh-prog/ndd-pp3-calibration.git ~/snrna-run
cd ~/snrna-run/snrna
git config pull.rebase false
```

**CPU environment** (everything except Evo 2, training and SpliceAI):

```bash
conda create -y -n snrna -c conda-forge -c bioconda python=3.12 infernal=1.1.5 openssl
conda activate snrna
pip install numpy scipy pandas openpyxl scikit-learn pysam ViennaRNA==2.7.0 pytest matplotlib pyBigWig pyliftover py2bit opentimestamps-client
export INFERNAL_BIN=$CONDA_PREFIX/bin          # Infernal 1.1.5 (cmalign, cmemit), used by structure/model_features/homology/a7
pip install torch --index-url https://download.pytorch.org/whl/cpu   # optional: only for the CPU dry run and tests/test_model.py
```

**GPU environment** (Evo 2 7B "light install", no Transformer Engine). The author used one `evo2` environment for every later step too:

```bash
conda create -y -n evo2 -c conda-forge -c bioconda python=3.12 infernal=1.1.5 openssl
conda activate evo2
pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
pip install packaging ninja psutil numpy pandas scipy scikit-learn pytest pysam matplotlib openpyxl ViennaRNA==2.7.0 pyBigWig pyliftover py2bit opentimestamps-client
MAX_JOBS=4 pip install flash-attn==2.8.0.post2 --no-build-isolation
pip install evo2
python -m evo2.test.test_evo2_generation --model_name evo2_7b     # downloads the ~15 GB weights and checks the install
export INFERNAL_BIN=$CONDA_PREFIX/bin
```

**SpliceAI environments.** `run_next.sh part5` and `spliceai_parallel.sh` create these automatically; to create them by hand:

```bash
conda create -y -n spliceai -c conda-forge python=3.10 && conda run -n spliceai pip install "tensorflow-cpu==2.15.1" "setuptools<70" spliceai
conda create -y -n spliceai-gpu -c conda-forge python=3.10 && conda run -n spliceai-gpu pip install "tensorflow[and-cuda]==2.15.1" "setuptools<70" spliceai
```

**Optional overrides.** `snrna.py` honours `SNRNA_DATA`, `SNRNA_WORK` and `SNRNA_RESULTS`; the defaults are `data/`, `work/` and `results/`.

## A1. Self-tests and code checks `[CPU]`

```bash
python -m pytest -q tests/test_snrna.py tests/test_labels.py tests/test_homology.py tests/test_final_eval.py   # no torch (~2 min; 21 passed, 1 skipped on a test machine)
python -m pytest -q tests                              # every test (test_model.py and test_evo2_features.py need torch)
SNRNA_TEST_CPU=1 bash run_a100.sh --dry-run            # the whole GPU pipeline on CPU with a mock Evo 2 and synthetic counts
python g2_power.py > results/g2_power_simulation.txt   # [FROZEN] G2 power/false-positive simulation, no real data
python comparators_report.py --shuffle                 # code check with permuted labels (numbers meaningless; re-run without --shuffle later)
```

## A2. Week 0: gene universe, sites and SGE checks `[CPU][NET]`

Remote files are read by tabix over HTTPS and cached under `data/cache/`, so re-runs work offline.

```bash
python snrna.py fetch     # SGE table, GENCODE v49, GIAB stratification (~100 MB)
python snrna.py genes     # 111-gene universe -> 87 pass coverage -> 55 after paralog QC (~6 min, gnomAD API)
python snrna.py sites     # all 17,385 SNVs + Roulette, CADD v1.7, gnomAD v4.1
python snrna.py sge       # SGE known-fact checks (stops on any mismatch)
python snrna.py vienna    # U4/U6 ddG recomputation, QC against the paper (ViennaRNA 2.7.0)
```

## A3. Gate G1: Evo 2 vs baselines on RNU4-2 SGE

```bash
python snrna.py evo2-prep     # [FROZEN] 8,192-nt windows -> evo2_run/evo2_inputs.tsv.gz (already committed)
```

`[GPU]` 24 GB, `evo2` environment. The scorer resumes, and the committed scores file is complete, so to re-score from scratch, write to a new file:

```bash
python evo2_score.py evo2_run/evo2_inputs.tsv.gz evo2_run/evo2_scores.tsv --limit 5      # smoke test + seconds per variant
python evo2_score.py evo2_run/evo2_inputs.tsv.gz evo2_run/evo2_scores.tsv                # all 521 (resumes; no-op on the committed file)
python evo2_score.py evo2_run/evo2_inputs.tsv.gz /tmp/evo2_rescore.tsv                   # independent full re-score (30-90 min on a 4090)
```

```bash
python snrna.py g1            # G1 metrics, bootstrap CIs, reproduction check, decision -> results/g1_report.txt
python snrna.py g1 --boot 2000
```

## A4. Gate G2: population selection signal `[CPU]`

```bash
python snrna.py structure                  # Rfam paired/unpaired per position via cmalign (needs INFERNAL_BIN)
python snrna.py g2-fit                     # [FROZEN] choose + freeze the mutation model on training genes -> results/g2_mutation_model.json
python snrna.py g2-test                    # one-shot G2 tests; refuses unless the frozen model file is committed
python snrna.py g2-test --allow-uncommitted   # tests only: skips the commit check
```

## A5. Model inputs, labels and comparators `[CPU][NET]`

All of these outputs are committed. `run_a100.sh` (full mode) checks code **and inputs** against `amendment_A2b_manifest.txt`, so if a rebuild changes any byte, the full run refuses.

```bash
python model_features.py                   # [FROZEN] model_inputs/: gene graphs + partner RNAs, Rfam profiles, ddG, paralog clusters, Evo 2 windows, sge_eval
python homology_transfer.py                # [FROZEN] RNU4-2 SGE -> RNU4ATAC transfer score h -> model_inputs/homology_rnu4atac.*
python labels.py nonblind                  # [FROZEN] labels/L2_nonblind.tsv (ClinVar + RNUdb, non-blind genes); --force overwrites
python curate_literature.py lit/tally_2026-10-14.tsv   # [FROZEN] labels/L2_curated_{additions,indels,rejected}.tsv
python comparators.py                      # [FROZEN] phyloP100, phyloP447, ReMM, GPN-MSA for every SNV -> comparators/comparators.tsv.gz
python v3_counts.py cohorts/ukb/*.csv > /dev/null      # gnomAD + UK Biobank counts -> model_inputs_v3/ (needed by v4, v3, gene_discovery_v3)
```

Don't pass `lit/a8_*.tsv` to `curate_literature.py`. Those files are blind-gene rows, and the script stops on them by design.

## A6. Seal the blind labels (before any model sees real counts) `[CPU][NET]`

```bash
export SNRNA_SEAL_PASS='choose-a-passphrase'        # or leave unset and you will be prompted
python labels.py seal                               # RNU2-2 + RNU5B-1 labels -> sealed/blind_labels.tsv.enc (never printed)
ots stamp sealed/SEALED_MANIFEST.jsonl              # OpenTimestamps proof of the sealing time
unset SNRNA_SEAL_PASS
```

Low-level form, for sealing any file: `python snrna.py seal <file> --pass-env SNRNA_SEAL_PASS`.

## A7. Pre-registered v2 selection model `[GPU]` (A100/H100 80 GB)

```bash
python evo2_features.py --limit-genes 1 --out evo2_run/features_dryrun   # smoke test on one gene
python evo2_features.py                                                  # every SNV, both strands, + per-position embeddings (5-10 h; resumable)
python evo2_features.py --batch 2                                        # if memory is short (the batch also halves itself on OOM)
SNRNA_EVO2_MOCK=1 python evo2_features.py --out evo2_run/features_mock   # CPU, fake Evo 2 (development only)
```

```bash
python train.py --benchmark                                              # time the largest configuration, estimate the full run
python train.py --synthetic-counts 1 --quick --max-lofo 3 --evo2 evo2_run/features_dryrun   # dry run on planted counts -> model_runs_synthetic/
python train.py --stage hp                                               # 48-fit hyperparameter grid
python train.py --stage main                                             # strict, loose, 14 leave-paralog-cluster-out fits
python train.py --stage ablations                                        # 12 ablations
python train.py --stage controls                                         # 3 negative controls
python train.py --stage all                                              # all four in order (3-8 h); 5 seeds each
python train.py --stage all --lofo all                                   # also hold out every unnamed paralog cluster
python train.py --stage all --device cpu                                 # force CPU
```

`train.py` refuses real counts until the labels are sealed (A6).

```bash
python evaluate.py                                                       # Aim 2, ablations, controls, gate G3 -> results/model_*
python evaluate.py --runs model_runs_synthetic --results results_dryrun --evo2 evo2_run/features_dryrun --boot 200   # dry-run evaluation
```

The one-command version of A6 and A7 (seal, Evo 2 features, train, evaluate, pack):

```bash
bash run_a100.sh --dry-run                  # 20-30 min; real Evo 2 on one gene, timing, synthetic training; no passphrase
tmux new -s snrna 'bash run_a100.sh; bash'  # full run, resumable (detach Ctrl-b d; reattach: tmux attach -t snrna)
```

## A8. v3 (exploratory): v2 retrained with gnomAD + UK Biobank counts `[GPU]`

```bash
python train.py --inputs model_inputs_v3 --out model_runs_v3 --stage all
python evaluate.py --runs model_runs_v3 --inputs model_inputs_v3 --results results_v3
```

One-command version: `tmux new -s v3 'bash run_v3.sh; bash'` (~3-4 h).

## A9. Non-blind follow-ups: comparators, preliminary calibration, gene discovery, Geuvadis `[CPU][NET]`

```bash
python comparators_report.py                                             # results/comparators_report.txt, comparators_aim2.tsv, comparators_units.tsv
python calibrate.py                                                      # preliminary Aim 4 PP3/BP4 (non-blind genes only)
python release.py --stage scores                                         # dated score bundle -> release/ (upload to Zenodo, record the DOI)
python gene_discovery.py                                                 # gnomAD only -> results/gene_discovery.{tsv,json}
python gene_discovery.py --inputs model_inputs_v3 --out results/gene_discovery_v3   # + UK Biobank
python gene_discovery.py --window 20                                     # window size is an option (20 nt is the default)
python gene_discovery_robustness.py > results/gene_discovery_robustness.txt   # run from snrna/; needs model_inputs_v3
mkdir -p geuv_data
curl -fsS --retry 4 -o geuv_data/junc.txt.gz https://ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1/analysis_results/GD462.JunctionQuantCount.45N.50FN.samplename.resk10.txt.gz
curl -fsS --retry 4 -o geuv_data/hg19.2bit https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
python geuvadis_carriers.py --junctions geuv_data/junc.txt.gz --twobit geuv_data/hg19.2bit   # 1000G genotypes streamed
```

One-command version: `tmux new -s aim4 'bash run_aim4.sh'` (needs `results/model_scores.tsv.gz` from `evaluate.py`).

## A10. v4: the stacked model, frozen before unsealing (amendment A5) `[CPU]`

```bash
bash run_v4.sh                                   # refuses if results_v4/v4_manifest.txt already exists
python v4.py                                     # the same, directly -> results_v4/{v4_cv.tsv,v4_report.txt,v4_model.json,v4_scores.tsv.gz,v4_manifest.txt}
python v4.py --inputs-counts model_inputs_v3     # explicit counts source (the default uses model_inputs_v3 if present)
sha256sum results_v4/v4_scores.tsv.gz results_v4/v4_model.json   # compare with the 62-character hash suffixes in amendment_A5_manifest.txt
```

v4 matches the frozen hashes only with the same Evo 2 deltas (`evo2_run/features/deltas.tsv`) and the same `model_inputs_v3/`. The author says: never re-run v4 after the freeze.

## A11. The single unsealing `[SEALED]`

```bash
SNRNA_SEAL_PASS='passphrase' python final_eval.py --reason "final evaluation, PREREG s9"
```

This one command unseals, scores v2 on RNU2-2 and RNU5B-1, then runs `v4_blind_eval.py`, `inheritance_blind_eval.py` and `calibrate.py --blind`. The pieces, if one of them stops (the author re-ran some by hand on the already-decrypted labels):

```bash
python snrna.py unseal sealed/blind_labels.tsv.enc --reason "final evaluation, PREREG s9" --pass-env SNRNA_SEAL_PASS
python v4_blind_eval.py --blind sealed/blind_labels.tsv           # hash-checks results_v4, then blind AUROC v4 vs CADD
python inheritance_blind_eval.py --blind sealed/blind_labels.tsv  # A6: dominant vs recessive from frozen local constraint
python calibrate.py --blind sealed/blind_labels.tsv               # final Aim 4 calibration
```

`final_eval.py` options: `--boot 2000 --calib-boot 1000 --runs --inputs --labels --results --evo2 --sealed`.

## A12. After the unsealing `[SEALED]`

```bash
python release.py --stage final                                 # scores + calibrated PP3/BP4 thresholds and evidence labels -> release/
python v4_ablation.py --blind sealed/blind_labels.tsv           # drop-one / single-feature v4 on blind RNU2-2 -> results/v4_ablation.tsv
```

One-command version: `tmux new -s afk 'bash run_afk.sh'`. It runs both of the above, then `run_v3.sh` (~5-6 h).

## A13. Post hoc checks on blind RNU2-2 `[SEALED][NET]`

```bash
python circularity_check.py --blind none --fetch-only           # any machine: ClinVar submission evidence -> sealed/blind_evidence.tsv
python circularity_check.py --blind sealed/blind_labels.tsv --fetch   # or fetch for the blind P/LP set itself
python circularity_check.py --blind sealed/blind_labels.tsv     # reuses the evidence -> results/circularity_check.{txt,tsv}
python blind_robustness.py --blind sealed/blind_labels.tsv      # position-only baselines, within-region AUROC, AC>=20/50/100, rarity-only null
python v4_calibrate.py --blind sealed/blind_labels.tsv          # out-of-sample v4 PP3/BP4 (pooled + per gene) + 279 ClinVar VUS
```

One-command version: `bash run_next.sh part1`. It runs these three, then `release.py --stage final` and `a7_ncrna.py score`.

## A14. Amendment A7: zero-shot on non-snRNA RNA genes (RMRP, TERC, SNORD118, RNU7-1) `[NET]`

```bash
python a7_ncrna.py score                                        # label-free scores -> a7/a7_scores.tsv.gz + a7/a7_manifest.txt
git add a7/ && git commit -m "A7 scores frozen" && git push     # eval refuses until the manifest is committed
ots stamp a7/a7_manifest.txt
python a7_ncrna.py eval --acmg                                  # -> results/a7_report.txt (--acmg also reads ClinVar ACMG codes; slow)
python a7_ncrna.py eval --boot 2000
```

`run_next.sh part2` runs the `eval --acmg` step.

## A15. Amendment A8: second blind test on literature-only variants `[SEALED]`

```bash
python a8_literature.py map lit/a8_RNU2-2_2026-10-14.tsv lit/a8_RNU5B-1_2026-10-14.tsv   # [FROZEN] -> labels/A8_literature{,_rejected,_indels}.tsv
python a8_literature.py eval --blind sealed/blind_labels.tsv    # refuses unless the A8 labels are committed
python a8_literature.py eval --blind sealed/blind_labels.tsv --boot 2000
python extras.py a8_sensitivity --blind sealed/blind_labels.tsv # A8 under the curator's STRICT rule and the union
```

`run_next.sh` has no step for `a8_literature.py eval`, so run it by hand, after part 2 and before part 3.

## A16. Amendment A9: genome-wide rarity inflation, non-coding `[NET]`

```bash
python a9_inflation.py select                                   # [FROZEN] a9/positives.tsv (ClinVar fileDate 2026-10-11)
python a9_inflation.py select --refresh                         # re-downloads today's ClinVar; the selection will differ from the frozen one
python a9_inflation.py fetch --workers 8 --timeout 180          # gnomAD/CADD/GPN-MSA/ReMM (+phyloP447 via pyBigWig) in +-100 bp windows (1-2 h)
python a9_inflation.py eval                                     # -> results/a9_report.txt, a9_results.tsv
python a9_inflation.py eval --boot 1000
python extras.py a9_phylop                                      # phyloP447 from the UCSC REST API, frozen A9 evaluation re-run -> results/a9_phylop/
python extras.py mechanism --blind sealed/blind_labels.tsv      # [SEALED] per-gene hypermutability vs score failure
```

`run_next.sh part3` runs fetch and eval. `part4` runs `a8_sensitivity`, `mechanism`, `a9_phylop` and `figures.py`.

## A17. Amendment A10: missense (AlphaMissense) and SpliceAI inflation `[NET]`

```bash
python a10_extend.py select                                     # [FROZEN] a10/missense_positives.tsv, a10/blb_*.tsv
python a10_extend.py missense_fetch --workers 8 --timeout 180   # downloads data/AlphaMissense_hg38.tsv.gz if missing, then the windows
python a10_extend.py spliceai_prep                              # -> a10/spliceai_in.vcf
curl -sSL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip > data/hg38.fa   # ~950 MB gz
```

Run SpliceAI with **one** of these three:

```bash
conda run -n spliceai spliceai -I a10/spliceai_in.vcf -O a10/spliceai_out.vcf -R data/hg38.fa -A grch38 -D 50   # serial, CPU
bash spliceai_parallel.sh cpu 16                                # 16 parallel CPU chunks, merged into a10/spliceai_out.vcf
bash spliceai_parallel.sh gpu                                   # 4 chunks on the GPU (env spliceai-gpu)
```

```bash
python a10_extend.py spliceai_parse                             # max delta score per SNV
python a10_extend.py eval --boot 1000                           # -> results/a10_report.txt, a10_results.tsv
```

`run_next.sh part5` runs all of A17, using serial SpliceAI.

## A18. Amendment A11: re-score the published NCBoost benchmark `[NET]`

```bash
python a11_ncboost.py fetch --workers 8 --timeout 180           # liftover GRCh37 -> GRCh38 (pyliftover) + scores
python a11_ncboost.py eval --boot 1000                          # -> results/a11_report.txt, a11_results.tsv
```

`run_next.sh part6` runs both.

## A18b. Amendment A12: prospective blind tests on new patient variants (added 15 Oct 2026)

```bash
python a12_new_variants.py scan                                  # monthly: ClinVar -> a12/clinvar_scan_<fileDate>.tsv, a12/scan_log.tsv
python a12_new_variants.py map lit/a12_*.tsv                     # literature tallies made with lit/A12_TALLY_PROMPT.md
python a12_new_variants.py build --blind sealed/blind_labels.tsv # [SEALED] -> labels/A12_labels.tsv; commit it
python a12_new_variants.py status                                # positives per unit vs minimum (no scores)
python a12_new_variants.py eval --unit U1_RNU6_adRP --blind sealed/blind_labels.tsv   # once per unit; also U2_U5, U3_RNU2-2, U4_other, U5_trained
```

One-command version: `bash run_next.sh part7`.

## A18c. Amendment A13: other published benchmarks (TraitGym first) and RarityBench

```bash
python a13_benchmarks.py freeze traitgym_mendelian               # already frozen and committed
python a13_benchmarks.py freeze mybench --tsv my.tsv --build GRCh37 --group match_col   # any other benchmark
python a13_benchmarks.py fetch traitgym_mendelian                # [NET] needs pysam (Linux)
python a13_benchmarks.py eval traitgym_mendelian
python raritybench.py                                            # one table across A9, A10, A11, A13 + fig9
python -m pytest -q tests/test_a12_a13.py
```

One-command version: `bash run_next.sh part8`.

## A18d. Amendment A14: PP3/BP4 recalibration with rare benign variants

```bash
python a14_recalibration.py select --vcf data/a12_clinvar_<date>.vcf.gz   # or any ClinVar GRCh38 VCF; commit a14/ before fetch
python a14_recalibration.py fetch                                # [NET] gnomAD AC, CADD, GPN-MSA, AlphaMissense (~2 h)
python a14_recalibration.py fetch --revel revel_with_transcript_ids.csv   # with REVEL v1.3, if downloaded
python a14_recalibration.py eval                                 # -> results/a14_report.txt
```

One-command version: `bash run_next.sh part9` (set `REVEL=<path>` first to include REVEL).

## A18e. Amendment A15: null validity (simulation) and bootstrap intervals

```bash
python a15_null_validity.py simulate                             # CPU, no data (already run and committed)
python a15_null_validity.py ci a9                                # window bootstrap; also a10, a11, a13:traitgym_mendelian
python a15_null_validity.py ci a11 --boot 200 --draws 50
```

One-command version: `bash run_next.sh part10`.

## A18g. Amendment A16: PP3 evidence against the rare-variant population

```bash
python a16_population_calibration.py sample                      # 3,000 random 1-kb windows; commit a16/ before fetch
python a16_population_calibration.py fetch --vcf data/a12_clinvar_<date>.vcf.gz   # [NET] ~2-3 h; add --revel <csv> if A14 had REVEL
python a16_population_calibration.py eval                        # -> results/a16_report.txt, a16_bands.tsv
python paper_tables.py                                           # paper/tables.md + fig10 forest plot from results/
```

One-command versions: `bash run_next.sh part11` (A16), then `bash run_next.sh part12` (RarityBench, tables, figures).

## A18f. Put the DGX-only outputs into git

```bash
bash sync_dgx_outputs.sh --list                                  # what would be added; prints full v4 hashes
bash sync_dgx_outputs.sh                                         # commit + push, or a tarball if push fails
```

## A19. Figures, release and lookup page

```bash
python figures.py --blind sealed/blind_labels.tsv               # 8 figures (PDF + PNG) -> paper/figures/; draws only those whose inputs exist
python release.py --stage final                                 # Zenodo bundle release/snrna-vep-v2-final-<date>.zip
python make_lookup.py                                           # release/lookup.html (static variant lookup page)
python make_lookup.py --release release/snrna-vep-v2-final-<date>/scores.tsv.gz --thresholds results/v4_calibration_thresholds.tsv --out release/lookup.html
```

## A20. The rarity-only null on your own benchmark

```bash
python rarity_null.py my_scores.tsv --score cadd_phred --draws 1000 --seed 0   # TSV: the score column + class in {pos, neg, pool}
```

From Python: `from rarity_null import rarity_null; rarity_null(score_pos, score_neg, score_pool, draws=1000, seed=0)`.

## A21. Integrity: timestamps and frozen hashes

```bash
ots upgrade *.ots                                               # Bitcoin-anchor the OpenTimestamps proofs (a few hours after stamping); also sealed/*.ots and a7/*.ots once they exist
ots verify prereg_manifest.txt.ots                              # check the pre-registration timestamp
ots info code_manifest.txt.ots
for m in *manifest*.txt; do echo "== $m"; grep -E '^[0-9a-f]{64}  ' "$m" | grep -v '\.md$' | sha256sum -c --quiet; done
```

On the current commit, the expected mismatches are files edited after an earlier freeze. Each edit is logged in PREREG s14:

- `release.py` and `curate_literature.py` against the A2b manifest
- `evo2_score.py` against `code_manifest.txt`
- the A2 code, superseded by A2b
- the `a7/` and `results_v4/` files, which are not in git

## Runner scripts (shortcuts for the steps above)

| Command | What it runs |
|---|---|
| `bash run_a100.sh --dry-run` / `bash run_a100.sh` | A6 seal → A7 Evo 2 features → train → evaluate → pack |
| `SNRNA_TEST_CPU=1 bash run_a100.sh --dry-run` | The same on CPU with a mock Evo 2 |
| `bash run_v3.sh` | A8 |
| `bash run_aim4.sh` | A9 (comparators, calibrate, release scores, gene discovery, Geuvadis) |
| `bash run_v4.sh` | A10 |
| `bash run_afk.sh` | A12 + v3 |
| `bash run_next.sh part1` | A13 + `release.py --stage final` + `a7_ncrna.py score` |
| `bash run_next.sh part2` | `a7_ncrna.py eval --acmg` |
| `bash run_next.sh part3` | A9 fetch + eval |
| `bash run_next.sh part4` | `extras.py a8_sensitivity`, `mechanism`, `a9_phylop` + `figures.py` |
| `bash run_next.sh part5` | A10 missense + SpliceAI |
| `bash run_next.sh part6` | A11 |
| `bash spliceai_parallel.sh cpu 16` / `gpu` | Parallel SpliceAI for A10 |

**Order the author used:** A1 → A2 → A3 → A4 → A5 → `run_a100.sh` (A6–A7) → A9 → A10 → A11 → A12 → A13 → A14 → A15 → A16 → A17 → A18 → A19.

---

# Part B: NDD de novo PP3/PM1 study (repo root)

## B0. Setup and path fixes

```bash
cd ~/snrna-run
conda env create -f environment.yml          # env "ndd-gate": python 3.12, numpy, pandas, scipy, matplotlib, bcftools, htslib, pyliftover
conda activate ndd-gate
pip install scikit-learn openpyxl            # NOT in environment.yml, but fullstudy.py (aim4*, clinvar, t13) and the .xlsx inputs need them
```

The code hard-codes the author's own paths: an external drive `/mnt/x/ndd_gate_data`, a OneDrive folder and `/root/ndd_pp3_gate`. Point them at your machine:

```bash
mkdir -p ~/ndd_ext
sudo mkdir -p /mnt/x && sudo ln -sfn ~/ndd_ext /mnt/x/ndd_gate_data      # EXT folder read by fullstudy.py and assemble_public3.py
sed -i "s#/mnt/c/Users/saite/OneDrive/Desktop/Neuro/G-Test-LOCKIN/ndd_pp3_gate/ndd_pp3_gate#$PWD#g; s#/root/ndd_pp3_gate#$PWD#g" fullstudy.py make_*.py results_fullstudy/fig_aims.py
mkdir -p results_fullstudy_backup work/cvt work/pg work/erepo work_erepo work_af work_clinvar
```

If Windows has an X: drive, WSL already mounts it at `/mnt/x`; then put `ndd_gate_data` on that drive instead.

## B1. Self-test (synthetic data with known answers, ~4 min) `[CPU]`

```bash
bash selftest.sh                             # must end with ALL SELF-TESTS PASSED (clean -> GO, biased -> KILL, CpG -> rescued)
```

## B2. Feasibility gate `[CPU][NET]`

```bash
bash get_data.sh                             # DeNovoWEST @ pinned commit, AlphaMissense hg38 (md5-checked), GENCODE v32, hg19->hg38 chain (~0.8 GB)
python gate.py check-inputs                  # every line must say OK
python gate.py prep-dnv                      # lift 45,221 de novos GRCh37 -> GRCh38
python gate.py prep-am                       # split AlphaMissense (~5 min); on MemoryError: --force --chunksize 500000
python gate.py prep-gtf                      # coding-sequence models + region files for Roulette
CHROMS="21 22" bash get_roulette.sh          # smoke test: 2 chromosomes (~1.9 GB)
python gate.py sites --chroms 21-22
python gate.py run --allow-partial --qc-only # pipeline QC only; computes no calibration number
PARALLEL=2 bash get_roulette.sh              # all 22 autosomes (~73 GB downloaded; resumable; PARALLEL=1 if disk < 25 GB; KEEP_VCF=1 keeps VCFs)
python gate.py sites                         # all chromosomes (~20-40 min; --force rebuilds)
python gate.py run                           # verdict -> results/gate_report.txt (--n-boot N to change the bootstrap count)
```

Steps `get_data.sh` through `gate.py run` in one command: `PARALLEL=2 bash run_all.sh`.

## B3. External data for the full study (download by hand into `~/ndd_ext` = `/mnt/x/ndd_gate_data`)

| Put it at | What it is |
|---|---|
| `dbnsfp/dbNSFP4.1a_variant.chr{1..22}.gz` | dbNSFP 4.1a (Zenodo 4323592; academic licence) |
| `dbnsfp54/dbNSFP5.4a_grch38.gz` | dbNSFP 5.4a, md5 `d2ad53e1a974208fea6050117be8fee1` (or the byte-range `part.*` pieces for `scores54_stream`) |
| `aims/gnomAD_v2.1.1_transcripts_with_rmc.tsv`, `aims/gnomAD_v2.1.1_transcripts_without_rmc.tsv` | gnomAD v2.1.1 regional missense constraint |
| `aims/ccrs.autosomes.v2.20180420.bed.gz` | Constrained coding regions (Havrilla 2019) |
| `aims/bm_s003_labels.csv`, `aims/bm_s005.xlsx` | Badonyi & Marsh mechanism labels (supplementary S3, S5; Aim 2) |
| `gnomad_constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` | gnomAD v2.1.1 gene constraint |
| `replication/zhou2022_dnv_combined.tsv` | Zhou 2022 autism de novos (Supplementary Data 1+2; columns IID, Chrom, Position, Ref, Alt, Pheno, Cohort) |
| `jin2017/jin2017_supplementary_tables.xlsx`, `denovodb/denovo-db.non-ssc-samples.variants.v1.6.1.tsv.gz` | Inputs to `assemble_public3.py` (third cohort) |

In the repo folder:

| Put it at | What it is |
|---|---|
| `work_clinvar/clinvar.vcf.gz` | ClinVar GRCh38 VCF (ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) |
| `work/cvt/submission_summary_2026-10.txt.gz`, `work/cvt/variant_summary_2026-10.txt.gz`, `work/cvt/submission_summary_2023-01.txt.gz` | ClinVar tab-delimited files, current and the Jan 2023 archive |
| `work_erepo/erepo_2026-10-07.tsv` (copy also to `work/erepo/`) | ClinGen Evidence Repository export |
| `work/pg/dms.zip` | ProteinGym v1.3 substitutions: `https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip` |
| `work_af/af_structure_v6.tsv.gz`, `work_af/uniprot_human_sp.tsv.gz` | Built in B4 (AlphaFold + UniProt) |

## B4. Full study: engineering steps (build the site universe and column store) `[CPU]`

Every `fullstudy.py` command accepts `--data --work --fs --results`. The ones below also take `--chroms 1-22 --force`.

```bash
python fullstudy.py universe                 # per-site coding universe (must equal the gate's aggregates)
python fullstudy.py scores                   # dbNSFP 4.1a rankscores at missense sites
python fullstudy.py regions                  # gnomAD RMC + CCR lifted to GRCh38, per site
python fullstudy.py colstore                 # one .npy per column per chromosome
python fullstudy.py validate                 # engine must reproduce the gate's AlphaMissense table
python fullstudy.py scores54                 # dbNSFP 5.4a A3 tools into the column store (whole file), OR:
python fullstudy.py scores54_stream          #   stream the downloaded byte-range pieces (md5-checked in the same pass)
python fullstudy.py scores54_tail --chroms 22 --from-piece part.7   # re-extract late chromosomes from trailing pieces only
python fullstudy.py dom54_stream             # InterPro domain per site from dbNSFP 5.4a (amendment A4.2)
curl -s "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,xref_ensembl&format=tsv&compressed=true" -o work_af/uniprot_human_sp.tsv.gz
curl -sL https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar | python af_structure.py work_af/af_structure_v6.tsv.gz   # streamed, never stored; check the tar's current name in that folder (the code expects v6)
python fullstudy.py af_cols                  # AlphaFold pLDDT / burial into the column store
```

## B5. Full study: analyses, in the order they were run (dependencies respected) `[CPU]`

```bash
python fullstudy.py aim5                     # data-derived priors
python fullstudy.py aim2                     # mechanism-specific calibration (AlphaMissense)
python fullstudy.py aim2 --tools             # optional: every Aim 3 tool too
python fullstudy.py aim1                     # headline: PP3 within gene / RMC / CCR
python fullstudy.py aim1x                    # exploratory pseudo-region control
python fullstudy.py aim3                     # all 36 dbNSFP 4.1a tools
python fullstudy.py aim3_report              # ranking + thresholds in raw-score units
python fullstudy.py aim1tools                # A1 nested levels for every tool
python fullstudy.py a1_report                # -> a1_nested.tsv (needed by clinvar, a4ctrl, a10_post)
python fullstudy.py replicate                # A2 autism replication + sibling negative control (needs zhou2022 file)
python fullstudy.py v54                      # A3 dbNSFP 5.4a tools + popEVE (DDD only)
python fullstudy.py aim4                     # Aim 4 v1 null-contrast model
python fullstudy.py aim4v2                   # v1 refit + v2 GBDT/MLP/ensemble (A5; ~8 GB RAM)
python fullstudy.py m4_nested                # Aim 4 models through levels A-E on DDD
python fullstudy.py aim4v3                   # A7 gene-matched null model
python fullstudy.py aim4v4                   # A8 all-centre model, autism-only test
python assemble_public3.py                   # builds EXT/public3/public3_dnv_snv_grch37.tsv
python fullstudy.py public3                  # A9 third public cohort
python fullstudy.py clinvar --clinvar work_clinvar/clinvar.vcf.gz   # A4 ClinVar impact -> a4_clinvar_matched.tsv.gz
python make_clinvar_resource.py results_fullstudy   # Supplementary Data 1 (needed by a10_post)
python fullstudy.py overlap                  # A10 PP3 x PM1/PP2 overlap -> a10_overlap.tsv
python fullstudy.py a10_post                 # A10.2 ClinVar sensitivity + A10.3 model equivalence
python fullstudy.py a11 --part strata --clinvar work_clinvar/clinvar.vcf.gz
python fullstudy.py a11 --part sweep --clinvar work_clinvar/clinvar.vcf.gz
python fullstudy.py a11 --part cpg_down --clinvar work_clinvar/clinvar.vcf.gz
python fullstudy.py a11 --part cpg_up --clinvar work_clinvar/clinvar.vcf.gz
python fullstudy.py a11 --part noclinvar --clinvar work_clinvar/clinvar.vcf.gz
python fullstudy.py a11sim --reps 200        # A11.5 planted-signal simulation (the author used 100)
python fullstudy.py a12 --part decomp
python fullstudy.py a12 --part classify
python fullstudy.py a12 --part lto
python fullstudy.py a12 --part allr
python fullstudy.py a12 --part plddt
python fullstudy.py a12 --part burial
python fullstudy.py a12 --part gclust
python fullstudy.py a12 --part gclust2
python fullstudy.py a13                      # overlap within disease genes
python fullstudy.py aim4v5                   # A15 label-free-input model
python fullstudy.py a14                      # A14 ClinVar within-gene AUROC of the de novo models (needs aim4v5)
python fullstudy.py a6                       # A6 autism models, sibling control, bootstrap reruns (needs aim4v2 + aim4v3)
python fullstudy.py a6_report                # -> a6_retention_ci.tsv
python fullstudy.py a4ctrl --part all        # or --part nodom, then dom, then report
python fullstudy.py hier                     # A16 hierarchical gene -> window -> PP3 model
python make_erepo_truth.py work_erepo/erepo_2026-10-07.tsv            # expert truth -> erepo_truth_sites.tsv (needed by a16val)
python make_erepo_truth.py work_erepo/erepo_2026-10-07.tsv popfreq    # sensitivity without BA1/BS1/BS2/PM2
python fullstudy.py a16val --count-only      # feasibility count only
python fullstudy.py a16val                   # A16 temporal ClinVar validation (needs work/cvt files)
python make_mave_truth.py                    # ProteinGym truth -> mave_truth_sites.tsv (needed by t13, benchval)
python fullstudy.py t13                      # A17 T13 within-window Poisson GBDT
python fullstudy.py bench2                   # gene-clustered benchmark -> bench2.tsv
python fullstudy.py benchval                 # benchmark validity vs ClinVar + ProteinGym
python fullstudy.py oracle                   # oracle ceiling for retention (needs hier)
python fullstudy.py pm1pp2                   # PM1 x PP2 independence
python fullstudy.py triple                   # PP3 + PM1 + PP2 stacking
python fullstudy.py credit_dn                # PP3 credit beside PM1 by score band
```

## B6. Helper analyses (expert panels, labs, assays, tables) `[CPU]`

Several of these write to `results_fullstudy_backup/`.

```bash
python make_a12_5.py work_erepo/erepo_2026-10-07.tsv                  # A12.5 eRepo classifications: uncapped / capped / conditional
python make_stacking_audit.py work_erepo/erepo_2026-10-07.tsv         # stacking audit, all classifications
python make_stacking_audit.py work_erepo/erepo_2026-10-07.tsv missense
python make_stacking_cv.py work_erepo/erepo_2026-10-07.tsv            # 200 gene-split cross-validations
python make_erepo_eligibility.py work_erepo/erepo_2026-10-07.tsv      # selection into PP3 on expert truth
python make_parser_validation.py work_erepo/erepo_2026-10-07.tsv      # free-text ACMG parser precision/recall
python make_lab_stacking.py                                           # lab-practice replication (work/cvt 2026-10 files)
python make_review20_lab.py                                           # within-submitter credit + held-out PP2 rule test
python make_review21_lab.py                                           # tightened held-out test, frequency-handled credit
python make_mave_all.py                                               # PP3 x PM1 across all ProteinGym assays
python make_mave_curated.py                                           # PP3 x curated PM1 residues (needs work/erepo/ copy)
python make_bridge.py                                                 # estimand bridge to the ClinVar scale
python phase0_mave_coverage.py                                        # MaveDB coverage counts (API)
python make_benchmark.py                                              # supp_table7_benchmark.tsv
python make_calib_table.py results_fullstudy/a10_overlap.tsv          # supp_table4_conditional_pp3.tsv
python make_resource.py results_fullstudy                             # supp_table1_predictor_resource.tsv
python make_threshold_sensitivity.py results_fullstudy                # retention at top 1/2/10%
```

## B7. Figures `[CPU]`

```bash
python fig1_design.py                        # Figure 1 study design
python fig_benchmark.py                      # fig2_benchmark.png
python fig_credit.py                         # fig2_credit.png
python fig4_overlap.py                       # Figure 4 (reads a10_overlap_main.tsv in the repo root)
python fig_manuscript.py results_fullstudy   # manuscript figures 2-3 (default path is the author's \\wsl$ share)
python results_fullstudy/fig_aims.py         # fig_aims_1_2.png
```

## B8. Manuscript builders (optional)

These **overwrite** committed drafts; each reads the previous draft.

```bash
python build_v6.py && python build_v8.py && python build_v9.py && python build_v11.py
python build_supp_v12.py && python build_supp_v13.py
```

`sync_commit.sh` is the author's own Windows-to-WSL copy-and-commit helper. It isn't part of the analysis.
