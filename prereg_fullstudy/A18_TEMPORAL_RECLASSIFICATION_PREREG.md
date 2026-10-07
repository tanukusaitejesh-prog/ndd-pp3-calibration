# Pre-registration A18: does stacked computational + location evidence predict later downgrade? (ClinVar 2019 → 2026)

*Written 7 October 2026, before any outcome below was tabulated. Fixed by the push of this file to branch
`claude/a18-temporal-reclassification`; not deposited in an independent registry (as for A16/A17).
Only file schemas were inspected before writing (column names and class labels), not any outcome by exposure.*

## Question
If PP3 and PM1 are credited as independent when they are not (R < 1, main text), variants called
Likely Pathogenic with *both* kinds of evidence carry less true evidence than the points credit them with.
Over time such calls should be revised downwards more often than LP calls that rest on only one of the two.
This is an outcome-based test that uses neither de novo data nor the classifications' own codes as truth.

## Data (all public)
- **Baseline:** ClinVar GRCh38 VCF 2019-09-23 (`clinvar_20190923.vcf.gz`, gnomAD public bucket
  `gs://gcp-public-data--gnomad/resources/grch38/clinvar/`; md5 recorded at run time).
- **Outcome:** the October 2026 class in `results_fullstudy/supp_data1_clinvar_pp3_pm1_points.tsv.gz`
  (ClinVar 2026-10-04; ≥ 1 star; P/LP, B/LB or VUS), joined on GRCh38 position + alt allele.
- **Drop-out bound:** ClinVar GRCh38 VCF 2025-05-04 (same bucket), for 2019 variants absent from the 2026 table.
- **Exposure flags** (already in the 2026 table, defined as in the main text):
  - PP3: REVEL ≥ Supporting threshold (0.644) — **primary**, because REVEL (2016) was the predictor in use
    before 2019 and AlphaMissense did not exist; AlphaMissense ≥ 0.792 and VEST4 ≥ 0.764 are secondary.
  - PM1-like: `pm1_like` (regional missense O/E < 0.4 or CCR percentile ≥ 95).

## Cohort
Autosomal missense SNVs with CLNSIG = `Likely_pathogenic` in 2019, ≥ 1 star (criteria provided, single or
multiple submitters without conflicts, expert panel, practice guideline), present in the 2026 table.

## Outcome
**Downgrade** = 2026 class VUS or B/LB (vs P/LP).

## Exposure groups
- **Both:** PP3 and PM1-like.
- **One:** exactly one of PP3, PM1-like.
- **Neither.**

## Primary endpoint and prediction
Odds ratio of downgrade, **Both vs One**, Mantel–Haenszel stratified by gene (genes with both groups
represented), 95% interval from a gene-clustered bootstrap (2,000 replicates, seed 20261018).
**Prediction:** OR > 1.
**Main-text rule:** reported in the main text if the 95% interval excludes 1, in either direction
(a credible OR < 1 would argue against over-crediting in practice). Otherwise reported in the supplement
as a null, with its minimum detectable OR.

## Secondary (descriptive, no rule)
1. Unstratified downgrade rates for Both / One / Neither, with gene-clustered intervals.
2. Same primary analysis with AlphaMissense and VEST4 as the PP3 proxy.
3. Baseline 2019 *Pathogenic* (not LP) cohort: expected to show a smaller effect, because P calls sit further
   from the class boundary.
4. Logistic regression of downgrade on Both / One / Neither with 2019 review-status (stars) as a covariate,
   gene-clustered errors.
5. Drop-out bound: among 2019 LP variants absent from the 2026 table, the 2025-05 status (conflicting,
   VUS, B/LB, P/LP, not found). Exposure cannot be assigned to them, so the bound is reported as the
   range of the primary OR obtained by assigning all drop-outs that are downgraded in 2025 to Both or to One.

## Known limitations, stated in advance
- Exposure flags are predictor/constraint proxies, not the codes curators applied in 2019; misclassification
  of exposure biases the OR towards 1.
- PM1-like uses regional constraint published after 2019; curators then used domains/hotspots.
- Variants with both kinds of evidence may truly be more pathogenic; conditioning on the same baseline class
  (LP) is what makes the comparison informative, but residual confounding by other evidence is possible.
