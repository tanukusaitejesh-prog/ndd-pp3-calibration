# SFARI Base data request: draft research plan (for the applicant to complete and submit)

> **Fill in before submitting:** applicant name, position, institution, signing official, IRB/ethics determination, and the exact SPARK release names offered on SFARI Base at the time of application. Fields marked [ ] are yours. Everything else can be pasted, edited as you prefer.

**Project title**
Label-free calibration of missense variant evidence within genes using de novo variants in autism.

**Applicant** [Name, position, institution, email]
**Institutional signing official** [Name, title, email]
**Ethics** [IRB/ethics committee determination for secondary analysis of de-identified data; reference number]

## Data requested

- **Primary:** de novo variant calls (SNVs, with annotations) for SPARK probands and unaffected siblings from the most recent exome or genome release, restricted to **individuals not included in Zhou et al. 2022 (Nat Genet)**. [Insert release names, e.g. "SPARK iWES vX / WGS vY".]
- **Minimal phenotype fields:** affected status, sex, and an NDD/intellectual-disability flag if available. No identifiable data are requested.
- **If individual-level de novo calls are not distributed:** trio genotype data for the same individuals, sufficient to call de novo SNVs. We would use the consortium's recommended pipeline.

## Background and rationale

Clinical classification of missense variants (ACMG/AMP) combines computational predictions (criterion PP3) with gene- and region-level evidence (PP2, PM1).

**Prior work by the applicant.** Using 31,058 developmental-disorder trios and 12,316 autism probands, we showed:

- Most of the evidence missense predictors provide for de novo variants is information about the gene and the protein position.
- Separately, a model trained directly on the de novo excess, with a within-gene objective, performs comparably to AlphaMissense within genes. Its point estimates are higher in two independent cohorts, but not significantly so.

**What is needed.** Determining whether such a model is better requires one large, independent autism/NDD cohort not used in its development. SPARK probands sequenced after the Zhou et al. 2022 data freeze are the largest such resource.

## Aims

1. **Primary.** Test, in a single pre-registered analysis, whether a de novo-trained model (design frozen and timestamped before data access) provides more within-gene evidence than AlphaMissense for de novo missense variants in autism probands.
   - **Estimand:** the within-gene likelihood ratio of the top 2% of scores, from the excess of observed de novo variants over a per-site mutation-rate expectation.
   - **Inference:** paired Poisson bootstrap over individuals.
2. **Secondary.** Replicate, in an independent cohort, how much evidence existing predictors retain after conditioning on gene, protein position and regional constraint.
3. **Negative control.** The same analyses in unaffected siblings.

## Methods summary

- **Variants and sites.** De novo SNVs are lifted to GRCh38 if needed and annotated against a fixed site universe: canonical GENCODE v32 coding sites.
- **Duplicates.** Individuals overlapping previously analysed cohorts are removed by shared-variant matching.
- **Expected counts.** From the Roulette mutation-rate model, scaled by the cohort's synonymous de novos.
- **Governance.** All analysis decisions are frozen and timestamped (git + OpenTimestamps) before data access. The plan is deposited [OSF link], and code is public on GitHub at publication.

## Data security and use

- Data are stored on an encrypted, access-controlled workstation. They are not redistributed.
- Only aggregate results are published (counts and likelihood ratios by score group and gene class). No individual-level data are reported.
- The work uses de novo variant calls only, for the stated research purpose. There will be no attempt to re-identify participants.
- On project completion, or as required by the data use agreement, data are deleted.

## Expected outputs

- A peer-reviewed publication and a preprint.
- Public release of model scores for all possible missense variants.
- Acknowledgement of SPARK participants and the SFARI as required.

## Timeline

| Months | Work |
|---|---|
| 1 | Data access and the frozen analysis |
| 2–3 | Manuscript |
