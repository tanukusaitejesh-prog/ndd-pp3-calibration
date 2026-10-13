# snRNA-VEP v2: selection-supervised, clinically calibrated variant-effect framework for spliceosomal snRNA genes in NDD

13 October 2026 · @Tejesh

This is the study plan as written on 13 October 2026. The frozen, executable version of its aims, endpoints, thresholds and splits is `PREREG.md`; where the two differ, `PREREG.md` wins and the difference is listed in its deviation log.

## Executive summary

The plan builds the first clinically calibrated variant-effect framework for spliceosomal snRNA genes, trained on population selection instead of scarce pathogenic labels. An August 2026 expert guidance says no in silico tool is designed or calibrated for non-coding RNA genes and asks the community to build one (guidance).

**Why it matters.** CADD thresholds that capture every reported pathogenic SNV also flag 62% (RNU4ATAC) and 55% (RNU4-2) of population SNVs. Variants in one 18-nucleotide RNU4-2 region are estimated to explain 0.4% of individuals with NDD (Chen 2024).

**Core design.** Five moves turn a thin method paper into a calibrated, generalizing, discovery-capable framework:

1. Selection-supervised training on gnomAD v4.1 with an snRNA-specific mutation model.
2. Triangulated validation: saturation genome editing (SGE), patient variants, population data.
3. Held-out-gene, leave-family-out and blind-hotspot tests.
4. ClinGen-style PP3/BP4 calibration with bootstrap confidence intervals.
5. Co-development with the RNUdb authors, who already host the field's variant resource.

**Honest ratings (my judgment, 1 to 5).** Novelty 4, Q1 fit 4 (about 4.5 for any Q1 journal, 3.5 if solo), feasibility 4. All three are conditional on the kill gates; see the ratings section.

**Kill gates.**

- **G1, week 1:** zero-shot Evo 2 reaches Spearman 0.6 or more and AUROC 0.85 or more on SGE, so drop the new-model framing and keep calibration.
- **G2, by week 4:** the known RNU4-2 depletion and two blind hotspots are not recovered, so drop selection-supervision.
- **G3, week 6:** the model fails to beat CADD by 0.15 AUROC on held-out genes, so pivot to calibration-only or stop.

## Verified gap and evidence

The gap is real and quantified: on RNU4-2, CADD reaches AUC 0.65 and ViennaRNA stability 0.72, against 0.93 for SGE function scores, and no tool is calibrated for non-coding RNA genes. Searches through 13 October 2026 found no snRNA-specific predictor.

| Finding | Evidence | Source (read status) |
|---|---|---|
| No snRNA-specific predictor; the community is asked to build one | None designed or calibrated for non-coding RNA genes; tools need calibration before PP3/BP4 use (Recommendation 6) | Guidance, Aug 2026 (read in full) |
| Genome-wide scores fail on snRNAs | CADD AUC 0.65. A threshold that catches all 12 ReNU SNVs (18.89 or more) flags 56.4% (195 of 346) of population SNVs | SGE, Nature 2026 (read in full) |
| Structure energetics are not enough | ViennaRNA U4/U6 ΔΔG: AUC 0.72, rho −0.27 across the gene, strong only in Stem II (rho −0.79) | SGE paper |
| Functional truth exists for one gene | RNU4-2 SGE: 539 variants, replicate r 0.83 to 0.86, 151 depleted (score below −0.302) | SGE paper |
| SGE cannot separate dominant from recessive | Haploid growth readout; haploid-diploid r = 0.75, and Stem III variants turn neutral in diploid cells | SGE paper |
| Population signal exists but is confounded | gnomAD v4.1: 76,215 genomes, 0.73 SNVs per nucleotide (1.45 in disease-linked snRNAs), about 50-fold de novo rate; only 87 snRNA genes pass coverage | Guidance |
| Pol III genes mutate faster | About 10-fold higher rate; Roulette misses a subset of high-rate sites in RNU and tRNA genes | Seplyarskiy 2023 (abstract and legends only) |
| Generic foundation models are unvalidated on snRNAs | Evo 2 claims state of the art on ncRNA fitness (preprint); a 2026 benchmark found sensitivity on mitochondrial tRNA variants fell from 52/79 to 4/79 when flanking context was permuted | Abstract level |
| Constraint-supervised learning works elsewhere | FLARE maps predicted regulatory effects to evolutionary constraint in non-coding DNA | Marderstein 2026 (abstract level) |

The SGE paper's population labels came from UK Biobank (490,640 genomes) and All of Us (414,840), both restricted. This project can use only gnomAD v4.1, roughly 12 times fewer genomes, which is the main statistical risk (see risk register).

## Brutal self-critique

Seven things can sink this project; the first four are the most likely, and the gates in this plan exist to find out cheaply whether they will.

| Threat | Why it is real | Likelihood (my estimate) | Defence |
|---|---|---|---|
| Population selection is underpowered | gnomAD v4.1 has about 12 times fewer genomes than the UK Biobank plus All of Us data behind the SGE labels. Pol III hypermutation can mask purifying selection, and Roulette misses some high-rate sites | High | Pool sites by structural element and by U4/U4atac analogs; gate G2; fall back to calibration-only |
| Foundation models already suffice | Evo 2 and EVEE are open and take an afternoon to run on snRNA variants | Medium | Gate G1 runs them first; if they work, the paper becomes calibration plus failure analysis |
| Scoop | The guidance authors host RNUdb and are aligned with ClinGen; the SGE team can extend to other snRNAs, as its paper says | Medium to high | Contact the RNUdb authors in week 1; pre-register; post a preprint early |
| Labels are noisy by design | 63 of 75 SGE-depleted variants outside the RNU4-2 critical region are seen in population cohorts, because they cause recessive disease. Presence in gnomAD does not mean benign | High | Benign-like labels only from high allele frequency, never from presence; region-specific analysis |
| Calibration data are tiny | The SGE team capped PS3 at strong because only 17 gold-standard pathogenic variants were available | High (certain) | Report supporting-strength thresholds with wide intervals; do not promise moderate or strong |
| Mapping and annotation errors | Paralogs are 99 to 100% identical, RNU4-2 differs from RNU4-1 at 4 of 145 nt, and transcript definitions changed (RNU4-2 141 to 145 nt; RNU4ATAC shifted by 1 nt) | Medium | Work in hg38 coordinates only; exclude low-VAF and problematic-region genes; QC against known RNU4-2 results |
| Not top-tier as a solo method paper | Method novelty is thin; the novelty is calibration and generalization | High | Frame as calibration plus held-out-gene generalization; collaborate |

Sources: guidance for transcript changes and paralogs, SGE paper for the labels and the 17-variant cap.

## Aims, hypotheses and pre-registered success criteria

Four aims, each with one primary endpoint fixed before any model is trained, so that a null result is reportable and a positive one is credible. Thresholds below are my proposals and will be frozen in an OSF pre-registration in week 0.

| Aim | Hypothesis | Primary endpoint | Pass criterion |
|---|---|---|---|
| 1. Selection map | After Pol III-aware mutation adjustment, depletion at the RNU4-2 critical region (n.62-70, n.75-78) is recovered, and two blind hotspots are flagged | Observed-to-expected rare-variant ratio per structural element; permutation test for hotspot overlap | Critical region below the gene-wide median with p < 0.01; both blind hotspots in the top-decile windows (permutation p < 0.01) |
| 2. Selection-supervised model | Training on population selection alone predicts SGE and patient pathogenicity without pathogenic labels | Spearman against SGE (521 in-transcript variants); AUROC for 12 ReNU SNVs against 346 population SNVs | Beats the best baseline (CADD, ViennaRNA, Evo 2) by at least 0.10 Spearman, or reaches AUROC 0.85 or more |
| 3. Generalization | The model transfers to genes it never saw | AUROC, patient P/LP against frequent population variants, for RNU4ATAC, RNU2-2, RNU5B-1 (leave-one-gene-out) and U4 to U4atac (leave-family-out) | At least 0.15 AUROC above CADD, 95% bootstrap interval excluding zero |
| 4. Clinical calibration | A calibrated score reaches PP3 and BP4 supporting strength | Likelihood ratios with bootstrap intervals, ClinGen method | Lower interval bound above the supporting threshold in at least one snRNA gene; otherwise report as not yet usable |

Secondary endpoints: a dominant-versus-recessive region head (AUROC against Table 2 regions), calibration drift across gnomAD ancestry groups, and agreement between the learned selection map and the SGE depletion map.

Out of scope: wet-lab work, restricted cohorts (UK Biobank, All of Us, Genomics England), and any claim of clinical readiness.

## Data inventory and access status

Every input is public, and the one high-value resource still to tally is the set of published patient variants. Two things are not used because they are restricted: UK Biobank, All of Us and Genomics England data, and the EGA RNA-seq (EGAS50000000889).

| Resource | What we take | Access | Status |
|---|---|---|---|
| RNU4-2 SGE | 539 function scores, depletion classes | Supplementary Table 1; reads at ENA PRJEB87505; code | Verified in the opened paper |
| gnomAD v4.1 genomes | Allele counts and coverage at snRNA loci (76,215 genomes) | gnomad.broadinstitute.org/data | Genome count verified; file layout to confirm in week 0 |
| GENCODE v49 | 111 snRNA genes and 1,630 pseudogenes (spliceosomal) | Open | Counts from the guidance |
| Roulette | Per-site mutation-rate model | GitHub | Code link verified; precomputed hg38 rates to confirm |
| Rfam | U4 (RF00015) and U4atac (RF00618) share a clan; also U6, U6atac, U1, U11 | U4 clan, U6 clan, U1 clan | Verified: 32,178 U4 and 1,648 U4atac regions |
| RNU4ATAC variants | 238 ClinVar entries, only 48 annotated to RNU4ATAC | Supplementary Table 2 (ClinVar as of 4 March 2025) | Counts from the guidance; re-download current ClinVar |
| Patient variants, other genes | RNU4-2 (Chen 2024, Greene 2024, Nava 2025, Rius 2026), RNU2-2 (Greene, Jackson, Leitão 2025-26), RNU5B-1, RNU6 and RNU6ATAC | Published supplements | To tally; my estimate is 150 to 250 distinct P/LP variants overall |
| RNUdb | Hotspot regions, patient and population overlays, equivalent variants | Web tool, downloadable | Beta; fully working for RNU4-2 only |
| Evo 2 7B | Zero-shot baseline | GitHub | 7B runs without FP8 hardware per the README |
| EVEE | Baseline scores | evee.goodfire.ai | ClinVar variants only |
| CADD v1.7, ViennaRNA 2.7.0 | Baselines | Open | Versions as used in the SGE paper |

Threshold note: the SGE preprint used −0.39 for depletion; the published paper uses −0.302 (q < 0.01) and −0.9 for strong depletion. Use the published values throughout.

## Methods 1: harmonization, QC and label sets

Every variant is keyed by hg38 position and alleles, never by transcript number, because snRNA transcript definitions changed and numbering breaks across resources. The guidance reports that RNU4-2 was extended from 141 to 145 nt and RNU4ATAC from 127 to 130 nt and shifted by 1 nt.

1. **Gene set.** Take GENCODE v49 spliceosomal snRNA genes (111). Keep genes whose gnomAD v4.1 genome coverage is median depth 20 or more and minimum depth 10 or more along the whole gene; the guidance's method yields 87 genes.
2. **Paralog QC.** Drop genes with mostly low-VAF calls or more than 50% overlap with GIAB problematic regions, as the RNU2-2 systematic screen did (43 genes excluded). Treat RNU4-2 against RNU4-1 (4 of 145 nt differ) as a high-risk locus.
3. **Coordinate check.** Convert every published variant to hg38, then confirm the reference allele. Reproduce known facts as a test: the SGE paper's 5 exact RNU4-2/RNU4ATAC equivalents and its 6 unassignable RNU4ATAC variants.
4. **Splits.** Sealed blind genes (RNU2-2, RNU5B-1): their patient labels stay in a locked file, checksum recorded in the pre-registration, opened once at final evaluation. Other splits are leave-one-gene-out, leave-family-out (U4 to U4atac) and leave-region-out inside RNU4-2.

| Label set | Definition | Used for |
|---|---|---|
| L1 SGE | Continuous function score; depleted if below −0.302, strong below −0.9 | Evaluation only |
| L2 Patient P/LP | Published pathogenic or likely pathogenic variants, tagged by gene, region and inheritance (dominant or recessive) | Evaluation, calibration |
| L3 Benign-like | Allele count of 9 or more in gnomAD v4.1 (about the SGE paper's AF threshold of 5.5e-5); sensitivity at 5, 20 and 50 | Calibration |
| L4 Unlabeled | All other observed or possible variants; never called benign | Selection training |

L3 is a scaled version of the SGE paper's rule (combined allele count above 100 in UK Biobank plus All of Us, about 905,000 genomes). Presence alone never makes a variant benign: SGE found that 63 of 75 depleted variants outside the critical region appear in population cohorts.

## Methods 2: snRNA mutation-rate and selection model

A count model converts gnomAD allele counts into per-variant depletion once mutation rate is accounted for, and every later step depends on getting that rate right. For each possible SNV v, the observed allele count A is modeled as a negative binomial around an expected count that shrinks with selection:

$$A_v \sim \mathrm{NegBin}(\lambda_v, \phi), \qquad \lambda_v = \mu_v \, T \, e^{-f_\theta(v)}, \qquad \mu_v = r^{\mathrm{Roulette}}_v \, c_{g,k(v)}$$

Here r is the Roulette rate, c is an snRNA-specific correction for gene g and structural element k, T converts rate to expected count at the gnomAD sample size, φ absorbs unmodeled heterogeneity, and f is the learned non-negative depletion score.

1. **Baseline rate.** Start from Roulette per-site rates, which already capture the roughly tenfold Pol III effect on average.
2. **snRNA correction.** Fit c on the observed site-frequency spectrum of well-covered genes, because Roulette misses a subset of high-rate sites in RNU genes. Fit it on elements expected to be neutral (loops with no known disease variants) and test on held-out genes.
3. **Site-rate covariate.** Add the number of other alternative alleles seen at the same position, excluding the variant being scored. Multiple alleles at a site are a direct readout of local hypermutability.
4. **Pooling.** Share element-level effects across analogous elements (U4 and U4atac by structural alignment) with partial pooling, because single sites have too few expected counts to estimate alone.
5. **Mutation-model QC.** Predicted against observed variant counts per gene and per trinucleotide context; residuals must show no trend against position outside known constrained regions.
6. **Element-level map (gate G2).** Observed-to-expected ratios per structural element, with permutation tests. This map is a result in its own right.

The learned score f is what the deep model in the next section predicts from sequence, structure and evolution.

## Methods 3: deep model, baselines and ablations

The model is deliberately small: one encoder over sequence, structure and evolution, three heads, trained on a single GPU and judged against seven baselines it must beat. The selection head trains only on unlabeled population data; SGE, patient and benign-like labels are touched only for evaluation, calibration and a post hoc mechanism probe.

Only the selection head is trained end to end, and only on population counts; the mechanism probe is fitted afterwards on non-blind genes, and everything else is used only to score predictions.

- **Inputs per variant.** The reference and alternative sequence of the whole gene plus 100 nt of flank; per-position profiles and pairwise covariation from the Rfam family alignment; a base-pairing graph from the consensus secondary structure, including U4/U6 intermolecular pairs; and an element label (stem, loop, k-turn, Sm site).
- **Encoder.** A small graph-augmented transformer with pair bias from the base-pairing graph, starting from frozen embeddings of Evo 2 7B or an open RNA language model (chosen by ablation).
- **Heads.** Selection (the depletion score f, negative-binomial likelihood) is the only head trained end to end. Mechanism (dominant hotspot, recessive region or neutral) is a separate probe on frozen embeddings, fitted on RNU4-2 and RNU4ATAC region labels only so the blind genes stay sealed. Uncertainty is a five-seed ensemble.

| Baseline | What it tests |
|---|---|
| CADD v1.7 | The field's default; AUC 0.65 on RNU4-2 |
| ViennaRNA U4/U6 ΔΔG | Pure thermodynamics; AUC 0.72 on RNU4-2 |
| Conservation and covariation alone | Whether evolution is enough |
| Evo 2 7B zero-shot | Whether a generic foundation model is enough (gate G1) |
| EVEE | Same, on the ClinVar subset |
| Allele count alone | Whether the model adds anything over what gnomAD shows directly (SGE correlation with allele count: rho 0.29) |
| Position-only model | The floor: element class and distance to hotspot |

**Ablations.** Remove the structure graph; remove covariation; remove the site-rate covariate; swap negative binomial for Poisson; replace Roulette with a context-only rate; and train once on SGE labels directly. The last gives an upper bound on how much selection supervision costs.

Hyperparameters are tuned only on gene-held-out folds of the unlabeled data, and the blind genes stay sealed until the final run.

## Validation design

The central validity risk is circularity: the model trains on gnomAD allele counts, and the SGE paper's patient-versus-population test uses population variants as controls. Every test below therefore removes the evaluated gene or family from training first.

| Test | Data and metric | Leakage control |
|---|---|---|
| SGE concordance (strict) | Spearman against 521 in-transcript SGE scores; AUROC for 12 ReNU SNVs against 346 population SNVs; by region (T-loop, Stem III, Stem II, k-turn, Sm site) | Train with all U4 genes removed (RNU4-2, RNU4-1, RNU4ATAC). Loose version removes only RNU4-2 |
| Held-out-gene patient test | AUROC for patient P/LP against benign-like (L3) in RNU4ATAC, RNU2-2, RNU5B-1, plus RNU6ATAC, RNU12 and RNU6 if each has 10 or more distinct variants | Leave-one-gene-out; count distinct variants, not patients, to avoid recurrent-variant inflation |
| Blind hotspots | Rank and AUROC of 5-nt windows for hotspot membership: RNU2-2 n.1-67, RNU5B-1 n.37-44 (Table 2); permutation p | Both genes' patient labels sealed; model trained without them |
| Leave-family-out | Train without U4 and U4atac, test on RNU4ATAC ClinVar P/LP against B/LB | Independent of the SGE-based mapping, which left 6 variants unassignable |
| Mechanism head | AUROC for dominant versus recessive region class; SGE strong versus moderate depletion as a severity proxy | Probe fitted on RNU4-2 and RNU4ATAC regions only; RNU2-2 and RNU5B-1 excluded |
| Ancestry check | Calibration thresholds and likelihood ratios per gnomAD ancestry group | Reported as drift, not tuned away |
| Negative controls | Shuffled labels, random structure graph, model trained on pseudogenes only | Must collapse to chance |
| Pseudo-prospective | Freeze scores in a dated Zenodo release, score against variants published afterwards | Timestamp is the control |

All metrics carry 95% stratified bootstrap intervals. A result counts only if it holds in the strict version of its test.

## Clinical calibration and mechanism head

Calibration follows the ClinGen method used for missense predictors, and with roughly 100 to 250 positives the honest expected outcome is supporting-strength evidence with wide intervals. The method is that of Pejaver 2022, extended to more tools in Bergquist 2025, with strengths converted to points as in Tavtigian 2020.

1. **Score sets.** Calibrate CADD, ViennaRNA ΔΔG, Evo 2 and the new model on the same variants, so the comparison is fair. The new model uses cross-fitted scores: each gene is scored by a model trained without it.
2. **Local likelihood ratios.** Estimate the ratio of pathogenic to benign-like density along each score using L2 and L3, bootstrapped by distinct variant.
3. **Thresholds.** Set a PP3 or BP4 threshold only where the conservative end of the bootstrap interval clears that evidence strength; otherwise report the score as not yet usable.
4. **Priors.** Use the ClinGen default prior as the primary analysis and show sensitivity from 0.01 to 0.5, because the prior of pathogenicity is far higher inside hotspots than outside them.
5. **Context split.** Calibrate dominant and recessive regions separately wherever each has enough variants. The guidance warns that de novo findings must be read in both dominant and recessive contexts.
6. **Circularity audit.** Flag positives whose original classification used CADD or conservation, since they bias the comparison toward those tools, and report results with and without them.

**Mechanism head (secondary).** It predicts whether a variant sits in a dominant-hotspot, recessive-region or neutral class, fitted on the guidance's critical-region labels for RNU4-2 and RNU4ATAC only, and evaluated against SGE strong versus moderate depletion as a severity proxy. It is reported as a research aid, never calibrated for classification.

**Output format.** One row per possible SNV in each well-covered gene: score, interval, evidence points at each threshold, region class, and a coordinate-safe key, in a JSON layout compatible with RNUdb so the authors can adopt it without reformatting.

No clinical-validity claim is made. Laboratories would need to validate any use of the thresholds themselves.

## Kill gates and decision tree

Read top to bottom: each gate lets the plan continue or routes to a fallback that still yields a paper. A position-only model within 0.05 Spearman of the best baseline would also end the deep-model framing inside RNU4-2.

```
G1 (week 1)  Evo 2 zero-shot on RNU4-2 SGE: Spearman >= 0.6 AND AUROC >= 0.85?
             yes -> drop new-model framing; paper = calibration + failure analysis
             no  -> continue
G2 (week 4)  RNU4-2 critical-region depletion AND both blind hotspots recovered?
             no  -> drop selection supervision; paper = calibration-only
             yes -> continue
G3 (week 6)  Model beats CADD by >= 0.15 AUROC on held-out genes?
             no  -> pivot to calibration-only, or stop
             yes -> full framework paper
```

## 16-week timeline and deliverables

Gate G1 lands in week 1, G2 in week 4 and G3 in week 6, so the project knows by week 6 which paper it is writing. Deliverables: the OSF pre-registration, the element-level selection map, the blind-hotspot result, calibrated thresholds, a dated Zenodo score release and the preprint.

The original timeline and decision tree were figures; the table above and the one below are reconstructed from the gate weeks stated in the text.

| Weeks | Work | Gate or deliverable |
|---|---|---|
| 0 | Pre-registration; data inventory; patient-variant tally (non-blind genes); seal blind labels | OSF record, sealed-label checksum |
| 1 | Week-1 protocol (Appendix A); outreach note | G1 |
| 2–4 | Mutation model, snRNA correction, element-level map | G2, selection map |
| 5–6 | Encoder and selection head; held-out-gene tests | G3 |
| 7–8 | Blind hotspots, leave-family-out, negative controls, ablations | Blind-hotspot result |
| 9–11 | Calibration, ancestry drift, mechanism probe | Calibrated thresholds |
| 12 | Score release | Dated Zenodo record |
| 13–16 | Manuscript; preprint by week 14 | Preprint |

## Operational risks and contingencies

The scientific threats are in the self-critique above; these are the practical ones, and each has a pre-agreed fallback that still yields a publishable paper.

| Risk | Trigger | Contingency |
|---|---|---|
| RNUdb authors do not reply | No answer by week 4 | Proceed solo, keep outputs RNUdb-compatible, cite the guidance, preprint by week 14 |
| gnomAD v4.1 extraction is slow | More than 2 days to pull 87 loci | Use regional downloads from the browser for the 87 loci only |
| Roulette rates not precomputed for hg38 | Not found in week 0 | Recompute from the repository, or fit a context-only rate with a Pol III factor learned from the spectrum |
| No 24 GB GPU for Evo 2 7B | Out of memory or no access | Use a hosted endpoint (NVIDIA NIM serves Evo 2; confirm quotas) or the 1B model |
| Fewer than 150 patient variants after tally | Count below 100 distinct P/LP | Drop per-gene calibration; report the selection map, SGE concordance and one pooled calibration |
| Blind labels leak | Any look at RNU2-2 or RNU5B-1 patient labels before final run | Log it, treat that test as exploratory, and say so in the paper |
| Pre-registration drift | Any threshold changed after data are seen | Record in a deviations log; changed thresholds are reported as post hoc |
| Licensing and attribution | Redistributing third-party tables | Release only derived scores and code; cite and link sources (SGE paper is CC BY 4.0) |

## Scoop mitigation and collaboration

The scoop risk is highest from the people who wrote the guidance, so the plan treats them as collaborators to approach in week 1 with a concrete result in hand. The guidance team hosts RNUdb and presented its text to ClinGen's variant classification group. The SGE paper says its strategy can be replicated for other snRNAs, and both papers share a corresponding author.

1. **Week 0.** Post the pre-registration on OSF. It timestamps the design and the sealed-label checksum.
2. **Week 1, after gate G1.** Send the outreach note below with the plan and the G1 result. Offer RNUdb-compatible scores for every well-covered snRNA gene and an independent blind-hotspot test.
3. **Ask for three things.** Willingness to co-develop or at least avoid duplication; access to not-yet-published variants for an external test; authorship terms agreed in writing early.
4. **Week 4.** No reply means proceed solo, keep the output format compatible, and cite the guidance.
5. **Week 12.** Release the scores as a dated Zenodo record, and post the preprint by week 14 even if gates have forced the calibration-only framing.

Do not use unpublished data without agreement, and do not wait for permission to run the pre-registered analyses.

> Subject: Calibrated in silico scores for snRNA variants (PP3/BP4): offer to co-develop with RNUdb
>
> Dear Dr Whiffin,
>
> Your August 2026 guidance asks for in silico tools calibrated on snRNA genes. We are building one: a model trained on population selection in gnomAD v4.1 with a Pol III-aware mutation model, checked against the RNU4-2 SGE map and held-out patient variants, with ClinGen-style PP3/BP4 calibration.
>
> Our pre-registered plan seals the RNU2-2 and RNU5B-1 hotspots as a blind test. Early result: [one line on gate G1].
>
> We would release scores in an RNUdb-compatible format and would welcome co-development, or at least coordination to avoid duplicated work. Could we have a 20-minute call in the next two weeks? Plan: [link].
>
> Best regards,
> [name, affiliation]

## Manuscript plan, target journals and reviewer Q&A

Frame the paper as calibration plus held-out-gene generalization, with the selection model and the blind-hotspot result as the new science. Working title: *A selection-supervised, ClinGen-calibrated variant-effect framework for spliceosomal snRNA genes.* If gate G1 forces the pivot, retitle it around calibrating in silico evidence for non-coding RNA genes.

**Targets.** First AJHG, Genome Medicine or Genetics in Medicine; fallbacks Bioinformatics or Nucleic Acids Research. Precedent: Pejaver 2022 calibrated predictors in AJHG, and Ellingford 2022 set non-coding interpretation recommendations in Genome Medicine. Confirm each journal's current JCR or Scimago quartile for the submission year; quartiles shift.

| Figure | Content | Depends on |
|---|---|---|
| 1 | The gap and design: CADD 0.65, ViennaRNA 0.72, SGE 0.93 on RNU4-2; study schematic | Nothing |
| 2 | Mutation-model QC and the element-level selection map across snRNAs | Gate G2 |
| 3 | Model versus SGE, strict split, by region | Gate G1 |
| 4 | Held-out-gene and leave-family-out generalization | Gate G3 |
| 5 | Blind hotspots: RNU2-2 and RNU5B-1 | Gate G2 |
| 6 | Calibration: reliability curves, thresholds with intervals, ancestry drift | Weeks 9 to 11 |
| Extended | Resource format, RNUdb integration, negative controls | Week 12 |

| Likely objection | Prepared answer |
|---|---|
| Only RNU4-2 has functional truth | Three independent truths: SGE, patient variants in six or more genes, population selection; held-out-gene tests carry the claim |
| Presence in gnomAD does not mean benign | Benign-like labels use high allele count only; selection training uses counts, not labels; sensitivity at four thresholds |
| Why not simply use Evo 2? | Gate G1 reports it first; where it fails, the failure analysis is a result |
| Training and controls share gnomAD | Strict leave-family-out and gene-removed training; shuffled-label and pseudogene controls collapse to chance |
| Calibration rests on about 100 positives | Conservative interval rule, supporting-strength ceiling, pooled and per-gene results shown side by side |
| Paralog mapping errors | Coverage and problematic-region filters; reproduces the SGE paper's known counts |
| Is it clinically useful? | No validity claim; RNUdb-compatible output; pseudo-prospective test on later variants |
| Does it work equally across ancestries? | Drift reported per gnomAD ancestry group, not tuned away |

## Resources, reproducibility and release

The project needs one 24 GB GPU, no wet lab and about four months of one person's time, and everything it produces is released openly. Figures marked approximate are my estimates, not sourced.

| Item | Plan | Note |
|---|---|---|
| Compute | One GPU with 24 GB or more for Evo 2 7B scoring and encoder training; CPU for gnomAD extraction | A 24 GB GPU is reported to fit the 7B model (third-party note); tens of GPU-hours, approximate |
| Software | Python and PyTorch, tabix and bcftools, ViennaRNA 2.7.0, CADD v1.7, Evo 2 (pip install), Roulette | Versions as in the SGE paper where it names them |
| Pre-registration | OSF record in week 0: hypotheses, thresholds, splits, sealed-label checksum | Changes go to a deviations log |
| Reproducibility | Git repository, container image, fixed seeds, workflow manager; every figure regenerates from one command | Archived with the paper |
| Release | Dated Zenodo record with scores, weights and code; RNUdb-compatible JSON | Derived scores under CC BY 4.0 |
| Input licenses | The SGE paper is CC BY 4.0; gnomAD and Rfam are open | Redistribute derived scores only, with attribution |
| Ethics | Only public, aggregate or published data; no individual-level records | Confirm with the host institution |
| Cost | Cloud GPU time of a few hundred to low thousands of US dollars, approximate | No wet-lab or licensing costs |
| People | One computational lead; a clinical variant scientist to review the label tally; a statistician to review calibration | Reviews are a few hours each |

## Honest ratings and what moves them

The redesign earns novelty 4, Q1 fit 4 and feasibility 4 on my judgment, but only if gates G2 and G3 pass; until then treat novelty as 3.5 to 4. Only the any-Q1 score sits above 4 today; novelty and feasibility reach 4.5 only after the gates are passed.

| Criterion (1 to 5) | v1: solo method paper | v2: this plan | What moves it |
|---|---|---|---|
| Novelty | 3.5 (problem 4.5, method 2.5) | 4 | To 4.5 if both blind hotspots and the held-out-gene test pass; to 3 if Evo 2 passes gate G1 |
| Q1 fit, any Q1 journal | 3.5 | 4.5 | Clean pre-registered results and an early preprint |
| Q1 fit, top tier (AJHG, Genome Medicine) | 2.5 | 4 with RNUdb collaboration, 3.5 solo | Co-authorship with the guidance team; calibration in several genes |
| Feasibility | 3.5 | 4 | To 4.5 if the patient tally exceeds 200 and extraction is smooth; to 3.5 if gate G2 fails |
| Scoop risk | High | Medium | Collaboration agreed by week 4 |

My probability estimates, not derived from data: any Q1 journal about 85%; a top-tier journal about 40% with collaboration and about 25% solo.

Three events would pull a rating below 4: Evo 2 passing gate G1 (novelty falls to about 3); no selection signal at gate G2 (feasibility falls to 3.5 and the framework reduces to calibration); or the RNUdb team releasing a calibrated predictor first (novelty falls to about 2.5, and the project pivots to an independent calibration and failure analysis).

## Appendix A: Week-1 protocol

Week 1 is a five-day test of whether a foundation model already solves RNU4-2, and its output decides gate G1. Run it after the pre-registration is posted in week 0.

1. **Day 1, data and coordinates.** Download SGE Supplementary Table 1 and clone the Crick repository. The library spans GRCh38 chr12:120291753-120291903; map every variant to hg38 and confirm reference alleles. Labels: depleted below −0.302, strong below −0.9; haploid scores primary, diploid for sensitivity; noise ceiling is replicate r of 0.83 to 0.86. Keep the RNU4-2 context, since RNU4-1 differs at only 4 of 145 nt.
2. **Day 2, floors.** Score all variants with CADD v1.7 and ViennaRNA U4/U6 ΔΔG, and build the position-only and allele-count-alone floors. Published reference points on RNU4-2: CADD AUC 0.65, ViennaRNA 0.72, SGE 0.93.
3. **Day 3, foundation models.** Run Evo 2 7B zero-shot, taking the alternative-minus-reference log-likelihood over a window around the locus and averaging forward and reverse-complement scores. Look up EVEE for the ClinVar-overlap subset only.
4. **Day 4, evaluate.** Spearman against the 521 in-transcript SGE scores; AUROC for depleted versus normal; AUROC for the 12 ReNU SNVs against population SNVs taken from gnomAD v4.1, so the control list will be shorter than the paper's 346 and CADD must be re-scored on it; per-region results; bootstrap intervals.
5. **Day 5, decision.** Apply gate G1 (see Kill gates and decision tree). Then send the outreach note and record the result in the pre-registration log.

Do not use the preprint's −0.39 threshold; the published −0.302 supersedes it.

## References and links

Two sources were opened in full and carry most of the numbers; the rest were read at abstract or snippet level and are marked so. Checked 13 October 2026.

**Opened in full**

- Guidance for clinical variant classification in genes for spliceosomal snRNAs (medRxiv, 4 Aug 2026; preprint, not yet peer reviewed; doi 10.64898/2026.08.03.26359558).
- Saturation editing of RNU4-2 reveals distinct dominant and recessive disorders (Nature 654:429-436, 2026; doi 10.1038/s41586-026-10334-9), with analysis code (github.com/FrancisCrickInstitute/RNU4-2_Saturation_Genome_Editing).

**Abstract or snippet level**

- Chen 2024, RNU4-2 and ReNU syndrome (doi 10.1038/s41586-024-07773-7); Nava 2025, U4 and U5 snRNAs.
- Systematic snRNA screen and RNU2-2 variants (preprint version; Nature Genetics 2026, doi 10.1038/s41588-026-02547-5).
- Seplyarskiy 2023, Roulette and Pol III mutagenesis.
- Evo 2 (Nature 2026, doi 10.1038/s41586-026-10176-5), EVEE, and the Evo 2 blind-spot benchmark.
- FLARE, Marderstein 2026.
- Calibration: Pejaver 2022, Bergquist 2025, Tavtigian 2020, Brnich 2019, Ellingford 2022.
- Data: gnomAD downloads, Rfam clans U4, U6, U1, Roulette code, EVEE portal.

**Patient-variant sources to tally (not yet opened)**

- RNU4-2 recessive: Rius 2026, doi 10.1038/s41588-026-02554-6. RNU2-2: Greene 2025 (10.1038/s41588-025-02159-5), Jackson 2025 (10.1038/s41588-025-02209-y), Jackson 2026 (10.1038/s41588-026-02551-9), Greene 2026 (10.1038/s41588-026-02539-5).
- RNU4 and RNU6 in retinitis pigmentosa: Quinodoz 2026 (10.1038/s41588-025-02451-4). RNU6ATAC: Mendez 2026 (10.1016/j.xhgg.2026.100588), Johnson 2026 (10.1016/j.ajhg.2026.02.017). RNU4ATAC: Benoit-Pilven 2020 (10.1371/journal.pone.0235655).
