# Novelty check (13 October 2026)

**Question.** Before we build the model: is snRNA-VEP v2 new, and which claims can the paper make?

**Searches run.**
- Web, extended mode: snRNA predictors, snRNA constraint, non-coding PP3/BP4 calibration, selection-supervised models, RNA language models for variant effects.
- Europe PMC: every paper and preprint since 2024 with snRNA, RNU4-2, RNU2-2, RNU5B-1, RNU4ATAC, RNU6ATAC, RNU12, ReNU or RNUopathy in the title. That is 270 records, about 60 of them relevant to disease.
- Full text: the closest overlaps (Jackson 2025; the region-aware CADD paper).

## Verdict by claim

| Claim | Status | Closest prior work | How to frame it |
|---|---|---|---|
| **An snRNA-specific variant-effect predictor** | **Novel.** None found | The guidance (Aug 2026) says none exists and asks for one. Leitão 2026 prioritised snRNA variants with MobiDeep, a general non-coding meta-score (CADD, ReMM, GPN-MSA, phyloP, Cactus), neither snRNA-specific nor calibrated | Main claim |
| **ClinGen PP3/BP4 calibration for snRNA (or any ncRNA-gene) in silico evidence** | **Novel.** None found | Pejaver 2022 and Bergquist 2025 are missense only. The SGE paper calibrated its *functional assay* (PS3/BS3), not an in silico score. Region-aware CADD thresholds (NAR Genom Bioinform 2025) are for prioritisation, with no likelihood ratios | Main claim |
| **Training on population selection (counts plus a mutation model) instead of pathogenic labels** | **Not new as a method class** | LINSIGHT (2017), ExtRaINSIGHT (2022), JARVIS (2021), GeneBayes (2024, gene level), Gnocchi (2024, regional) | Claim only the *application*: per-variant selection learned at hypermutable Pol III/snRNA loci with an snRNA-aware mutation model and structure-aware inputs. Cite all five |
| **Depletion at the RNU4-2 critical region** | **Known** | Chen 2024 (UK Biobank: 26% of possible SNVs observed vs 78% elsewhere) | Positive control only |
| **Depletion at the RNU2-2 and RNU5B-1 hotspots ("blind" G2)** | **Known** | Jackson 2025 already mapped gnomAD v4 depletion in RNU2-2 (n.1–60) and RNU5B-1 (n.30–50) with a sliding-window method; Nava 2025 reported RNU5B-1/RNU5A-1 variants as depleted in UK Biobank | Replication under a mutation-rate-adjusted, pre-registered model. **Not** a discovery; drop "blind hotspot discovery" from the title and abstract |
| **Systematic mutation-adjusted selection map across 55 snRNA genes, and the dominant-depleted vs recessive-tolerated contrast** | **Partly new** | Leitão 2026 proposes a "gradient-of-impact" model; the SGE paper shows recessive-region variants in population cohorts within RNU4-2 | Systematic and quantitative across genes and Table 2 regions; secondary result |
| **Evo 2 fails on snRNA variants (ρ 0.30; T-loop ρ 0.00)** | **New result.** No prior Evo 2 benchmark on RNU4-2 SGE found | A 2026 benchmark reported an Evo 2 blind spot on mitochondrial tRNAs (abstract level) | Part of the failure analysis |

## What this changes

1. **Lead with the predictor and calibration.** The headline is the first snRNA-specific in silico tool calibrated to PP3/BP4, plus held-out-gene generalisation. That is unchanged and still open.
2. **Demote the blind-hotspot result.** It is a replication that validates the mutation model, not a discovery. `PLAN.md` rated novelty 4.5 if both blind hotspots passed; that upgrade no longer applies either way.
3. **Add baselines reviewers will ask for.**
   - MobiDeep, because it was used for snRNA screening in Leitão 2026.
   - GPN-MSA and ReMM (precomputed scores are public).
   - One RNA language model (RiNALMo or RNA-FM).

   Each is logged in `PREREG.md` before its results are seen, and none changes a gate.
4. **Scoop watch.** Three groups are the most likely to build an snRNA predictor:
   - the guidance and RNUdb team (Whiffin);
   - the French and Montpellier groups behind Nava 2025, Leitão 2026 and MobiDeep;
   - the SGE team (Findlay).

   Outreach in week 1, as planned, still makes sense.

## Updated rating (my judgment)

| | `PLAN.md` v2 | Now |
|---|---|---|
| Novelty | 4 | **4** for the predictor and calibration, which stay unscooped; the selection map drops to about 3. The overall story holds if G3 passes |
| Main risk to novelty | Scoop | Unchanged; plus a reviewer saying "selection-supervised = LINSIGHT/JARVIS", answered by framing and citing them |

## Sources

- Guidance for clinical variant classification in spliceosomal snRNA genes, medRxiv 2026, doi 10.64898/2026.08.03.26359558
- Chen et al. 2024, Nature, doi 10.1038/s41586-024-07773-7 (RNU4-2 depletion in UK Biobank)
- Jackson et al. 2025, Nature Genetics, doi 10.1038/s41588-025-02209-y (gnomAD v4 constrained regions in RNU2-2 and RNU5B-1)
- Nava et al. 2025, Nature Genetics, doi 10.1038/s41588-025-02184-4
- Leitão et al. 2026, Nature Genetics, doi 10.1038/s41588-026-02547-5 (MobiDeep prioritisation; gradient-of-impact model)
- Saturation editing of RNU4-2, Nature 2026, doi 10.1038/s41586-026-10334-9
- Small nuclear RNA genes in Mendelian disorders (review), Nature Genetics 2025, doi 10.1038/s41588-025-02440-7
- Genome region aware CADD thresholds for noncoding variant prioritization, NAR Genomics and Bioinformatics 2025 (PMC12629842)
- LINSIGHT: Huang et al. 2017, Nature Genetics, doi 10.1038/ng.3810
- ExtRaINSIGHT: Gulko and Siepel 2022, Nature Communications, doi 10.1038/s41467-022-31872-6
- JARVIS: Vitsios et al. 2021, Nature Communications, doi 10.1038/s41467-021-21790-4
- GeneBayes: Zeng et al. 2024 (bioRxiv 10.1101/2023.05.19.541520; Nature Genetics 2024)
- Gnocchi: Chen et al. 2024, Nature 625:92–100, doi 10.1038/s41586-023-06045-0
- Pejaver et al. 2022, AJHG (PMC9748256); Bergquist et al. 2025 (calibration of additional tools)
- Whole-genome discovery of pathogenic snRNA variants, iScience 2026, doi 10.1016/j.isci.2026.116814
- Further delineation of snRNA NDDs, Research Square preprint 2026, doi 10.21203/rs.3.rs-9609995/v1
