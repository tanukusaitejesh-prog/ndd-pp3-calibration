# NDD de novo PP3 calibration: feasibility gate

This gate decides whether the full study is worth doing. It checks two things:

- **Null check.** Does the Roulette mutation model reproduce the synonymous de novo count across mutability classes and protein contexts?
- **Precision.** Is the de novo excess big enough to give AlphaMissense likelihood ratios precise enough for ACMG PP3/BP4 calibration?

The criteria are fixed in `GATE_PLAN.md` before any result is seen.

## What you need

| | |
|---|---|
| OS | Linux or macOS. On Windows, use WSL2 (Ubuntu), see step 1. |
| RAM | 16 GB recommended (8 GB probably works) |
| Disk | ~25 GB free. Roulette is processed one chromosome at a time and each VCF is deleted after extraction. |
| Download | ~75 GB in total: Roulette ≈ 73 GB, everything else ≈ 0.8 GB. Downloads resume if interrupted. |
| Time | The download is the slow part (hours, so run it overnight). Compute is under an hour. |

## Steps

**1. Windows only: install WSL2.** Open PowerShell as Administrator and run `wsl --install`, then reboot. Open the **Ubuntu** app and do every remaining step inside it. Keep the project in your Linux home (`~/`), not under `/mnt/c`, because files there are much slower.

**2. Install Miniforge (conda)** if you don't have conda yet:
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b
~/miniforge3/bin/conda init bash && exec bash      # on macOS (zsh): ~/miniforge3/bin/conda init zsh && exec zsh
```

**3. Unpack the bundle and create the environment.** This takes about 5 minutes.
```bash
cd ~
# WSL: the zip is in your Windows Downloads folder; copy it into Linux first:
#   cp /mnt/c/Users/<your-windows-name>/Downloads/ndd_pp3_gate.zip ~/   (and: sudo apt-get install -y unzip)
unzip ndd_pp3_gate.zip && cd ndd_pp3_gate
conda env create -f environment.yml
conda activate ndd-gate          # do this in every new terminal
```

**4. Self-test.** This runs the whole pipeline, bcftools included, on synthetic data with known answers. It takes about 4 minutes.
```bash
bash selftest.sh
```
It covers three scenarios:

- **Clean data.** The pipeline must recover the true coordinates, counts and likelihood ratios, and say GO.
- **Detection bias.** De novos are under-called in high-AlphaMissense codons, and the gate must say KILL.
- **CpG miscalibration.** The pooled scaling must fail and the pre-registered per-class scaling must rescue it.

It must end with `ALL SELF-TESTS PASSED`. If it doesn't, stop and send me the output.

**5. Freeze the plan before touching real data.** This is your pre-registration timestamp. Read `GATE_PLAN.md` and edit section 6 if you want a different held-out centre. Then commit the folder to a **private** GitHub repo:
```bash
git init && git add gate.py *.sh *.md environment.yml tests && git commit -m "Pre-registered feasibility gate"
# create an empty private repo on github.com, then:
git remote add origin https://github.com/<you>/ndd-pp3-gate.git && git push -u origin HEAD
```
Uploading `GATE_PLAN.md` to OSF also works.

**6. Download the small inputs** (≈ 0.8 GB):
```bash
bash get_data.sh
```

**7. Prepare.** This takes about 10 minutes.
```bash
python gate.py check-inputs     # every line must say OK
python gate.py prep-dnv         # lift the 45,221 de novos from GRCh37 to GRCh38
python gate.py prep-am          # split AlphaMissense (~5 min)
python gate.py prep-gtf         # coding-sequence models + region files for Roulette
```

**8. Smoke test on two small chromosomes** (≈ 1.9 GB download). This shows pipeline QC only and no calibration numbers, so nothing is peeked at:
```bash
CHROMS="21 22" bash get_roulette.sh
python gate.py sites --chroms 21-22
python gate.py run --allow-partial --qc-only
```
Check that the `G0 pipeline QC` line says PASS. If it says FAIL, send me `results/gate_report.txt` before downloading more.

**9. Full Roulette download.** Leave this running overnight. If it stops, re-run the same command and it resumes:
```bash
PARALLEL=2 bash get_roulette.sh       # 2 chromosomes at a time; use PARALLEL=1 if disk < 25 GB
```

**10. Build all sites and run the gate.**
```bash
python gate.py sites                  # ~20-40 min, skips chromosomes already done
python gate.py run                    # ~5 min; prints the verdict
```

**11. Send me these files** from `results/`:
- `gate_report.txt`
- `lr_bins.tsv`
- `null_checks.tsv`
- `thresholds.tsv`
- `cumulative_lr.tsv`
- `qc.tsv`
- `qc_by_chrom.tsv`
- `run_info.json`
- `lr_curves.png`
- `null_checks.png`

They are aggregate counts only.

`bash run_all.sh` runs steps 6–10 unattended, if you prefer that after the smoke test.

## How to read the verdict

- **GO**: the null calibrates and the LRs are precise. Next come the full study (about 30 predictors, local-window thresholds, an ML model with DDD held out) and the preprint.
- **GO WITH CAUTION**: a WARN item needs explaining first, usually a mutability or context group at 0.85–0.90 or 1.10–1.15.
- **KILL**: Roulette does not reproduce synonymous de novos in some stratum, so the expected counts can't be trusted. Changing the model after seeing this would break the pre-registration.
- **UNDERPOWERED**: the signal is too imprecise for ACMG-scale claims.

The single most interesting number is the **within-gene LR** at AlphaMissense ≥ 0.99. It measures what the score adds beyond knowing the gene.

## Troubleshooting

- **`conda env create` cannot find bcftools (Apple Silicon Mac):** run `brew install bcftools htslib`, remove bcftools/htslib from `environment.yml`, and create the environment again.
- **A download stops or fails:** re-run the same command. Finished files are marked with `.ok` and skipped.
- **Download is very slow:** the Roulette server is a single university host. Use `PARALLEL=3` if you have disk space, or leave it overnight.
- **`MemoryError` in `prep-am`:** run `python gate.py prep-am --force --chunksize 500000`.
- **Anything else:** send me the full error text and the last 30 lines of output.

## What to cite

- **De novo mutations:** Kaplanis J. et al. *Evidence for 28 genetic disorders discovered by combining healthcare and research data.* Nature 586, 757–762 (2020), doi:10.1038/s41586-020-2832-5, Supplementary Table 1. The file used is the re-normalised version (`…_2025_09_18`) in the DeNovoWEST repository (github.com/HurlesGroupSanger/DeNovoWEST, formerly queenjobo/DeNovoWEST; MIT licence), at the commit given in `run_info.json`.
- **AlphaMissense:** Cheng J. et al. *Accurate proteome-wide missense variant effect prediction with AlphaMissense.* Science 381, eadg7492 (2023). Predictions: Zenodo doi:10.5281/zenodo.10813168 (CC BY 4.0).
- **Roulette:** Seplyarskiy V. et al. *A mutation rate model at the basepair resolution identifies the mutagenic effect of polymerase III transcription.* Nat Genet 55, 2235–2242 (2023), doi:10.1038/s41588-023-01562-0. Per-site rate files `N_rate_v5.2_TFBS_correction_all.vcf.bgz` from genetics.bwh.harvard.edu/downloads/Vova/Roulette. These are a public download, but no licence is stated, so cite the paper.
- **Calibration framework:** Pejaver V. et al. *Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.* Am J Hum Genet 109, 2163–2177 (2022).
- **Annotation:** GENCODE release 32 (Frankish A. et al. *GENCODE reference annotation for the human and mouse genomes.* Nucleic Acids Res 47(D1), D766–D773 (2019), doi:10.1093/nar/gky955). Liftover: UCSC hg19ToHg38 chain.
- **s_het:** Cassa C.A. et al. *Estimating the selective effects of heterozygous protein-truncating variants from human exome data.* Nat Genet 49, 806–810 (2017), doi:10.1038/ng.3831, Supplementary Table 1, via DeNovoWEST. This was checked on 2026-10-04: DeNovoWEST's `s_het` equals Cassa's all-population s_het for every shared gene, and Kaplanis 2020 cites Cassa for s_het. It is not Weghorn et al. 2019.
- **DDG2P:** Thormann A. et al. *Flexible and scalable diagnostic filtering of genomic variants using G2P with Ensembl VEP.* Nat Commun 10, 2373 (2019).
