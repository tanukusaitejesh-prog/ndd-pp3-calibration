# Novelty audit (14 October 2026, after the single unsealing)

Supersedes the claim-level verdicts in `NOVELTY.md` (13 October). Four parallel literature searches were run:
- snRNA predictors;
- PP3/BP4 calibration and clinical guidance;
- constraint, gene discovery and inheritance;
- functional screens.

Each search covered PubMed, Europe PMC, bioRxiv/medRxiv, Research Square, MaveDB and the IGVF portal, and each was told to look for prior work. Every reference below was opened unless marked UNVERIFIED.

## Verdict by claim

| # | Claim | Verdict | Closest prior work | What remains new |
|---|---|---|---|---|
| 1 | An snRNA-specific variant-effect predictor | **NEW** | D'Souza…Whiffin, medRxiv 2026 (10.64898/2026.08.03.26359558): no tool "designed or calibrated specifically" for ncRNA genes. MobiDeep (Bouazzaoui…Baux, Research Square 2026, 10.21203/rs.3.rs-8823759/v1) is a general noncoding meta-score; it "prioritised" ClinVar RNU4-2 variants, with no held-out gene and no AUROC. Leitão 2026 used it for annotation only | The whole claim |
| 2 | Frozen v4 beats CADD blind on RNU2-2 (0.89 vs 0.50) | **NEW (number); qualitatively expected** | De Jonghe…Findlay, Nature 2026 (10.1038/s41586-026-10334-9): CADD AUC 0.65 on RNU4-2 ReNU vs population. D'Souza 2026: RNU2-2 n.4G>A / n.35A>G have CADD ~11.5, which missense cutoffs read as benign. Nava 2025: no CADD difference | First quantitative blind benchmark on RNU2-2, and the first to include phyloP447, GPN-MSA, Evo 2 and ReMM on an snRNA |
| 3 | Rfam family profiles carry the signal (alone ~0.85 blind) where genome-wide conservation fails | **NEW for nuclear ncRNA**; concept precedented in mt-tRNA | MitoTIP (Sonney 2017, PLoS Comput Biol, 10.1371/journal.pcbi.1005867): cross-tRNA alignment. nAPOGEE (Bianco…Mazza, bioRxiv 2025, 10.64898/2025.12.20.695671) | Covariance-model column features as a replacement for phyloP at multicopy loci. **The most generalisable finding** |
| 4 | Selection-supervised per-variant model at snRNA loci | **NEW as a negative result** | LINSIGHT, ExtRaINSIGHT, JARVIS, Gnocchi (genome-wide). Chen 2024 Nature warns that mutation models are poorly calibrated for snRNAs. Seplyarskiy 2023 Nat Genet (10.1038/s41588-023-01562-0): Roulette misses the RNU excess | v2's failure on SGE (ρ 0.14) as an informative negative; cite Chen and Seplyarskiy as the prior expectation |
| 5 | ClinGen local-LR PP3/BP4 calibration, per snRNA gene | **NEW for nuclear snRNA genes.** Not "first ncRNA" | Villani…Spurdle AJHG 2024 (10.1016/j.ajhg.2024.05.002): 5′ regulatory, Tavtigian LR. Walker 2023 AJHG: SpliceAI. nAPOGEE: mt-tRNA/rRNA posteriors. Tejura…Fowler AJHG 2024: per-gene heterogeneity. The guidance (Rec. 6) allows PP3/BP4 only after snRNA calibration and computes no LR | It answers the guidance's own stated gap. A PM1 region LR would also be new; the guidance says one "can" be computed but computes none |
| 6 | ClinVar snRNA VUS reclassifiable by calibrated thresholds | **NEW** | De Jonghe (SGE on a handful of VUS) | First systematic count; reviewers will compare it against SGE for RNU4-2 |
| 7 | Pre-registration (OpenTimestamps) + AES-sealed labels opened once | **NEW in VEP evaluation** | CAGI: blind challenges. Zhang bioRxiv 2026 (10.64898/2026.09.14.751496): hash-frozen data, no pre-registration | A design contribution; keep it as a methods strength, not the headline |
| 8 | Scan of all snRNA loci for depleted windows; five NDD genes ranked top | **PARTLY DONE, high overlap** | Chen 2024 Nature, Fig. 4 / ED Table 4: 18-bp windows across 28 snRNAs in UKB; RNU4-1 central depletion already shown. **Messaoud…O'Donnell-Luria, medRxiv Aug 2026 (10.64898/2026.08.06.26358951)**: 10-bp windows over 41 genes + 1,613 pseudogenes in gnomAD v4.1, 71 regions, new candidates (RNU5A-1, RNU5F-1, RNU2-63P, RNU5E-6P). Jackson 2025: RNU2-2 and RNU5B-1 only | Mutation-model (Roulette) calibrated scan; quantitative disease-gene benchmark (AUROC 1.0, permutation p); top-window = disease-region overlap. **RNU4-1 is not new.** **Messaoud reports depletion in the adRP RNU6 genes, contradicting "RNU6 not enriched"; reconcile** |
| 9 | Whole-gene O/E fails because snRNAs are mutation hotspots | **PARTLY DONE** | D'Souza 2026: ~50× de novo rate, gnomAD SNVs/nt 0.73 vs intergenic 0.23. Seplyarskiy 2023: 7× SNV / 32× de novo excess, attributed to Pol III | The explicit demonstration that gene-level O/E ranks disease genes near the bottom. Do not claim transcription-associated mutagenesis for Pol II snRNAs; it is not established |
| 10 | Dominant regions depleted, recessive not | **PARTLY DONE** | Chen 2024 (RNU4ATAC weakly depleted). Rius…Whiffin Nat Genet 2026 (10.1038/s41588-026-02554-6) and Greene 2026 (10.1038/s41588-026-02539-5) define recessive regions | O/E for recessive vs dominant regions within the same gene |
| 11 | Inheritance mode predicted from local constraint | **NEW as a classifier; qualitative pattern in Chen 2024** | MOI-Pred (Petrazzini…Do, Cell Rep Methods 2024): missense only | Variant-level classifier. The non-blind AUROC is close to circular (the classes occupy different regions); it needs a region- or gene-held-out test |
| 12 | Healthy-carrier minor-intron readout (Geuvadis) | **NEW, null** | Merico 2015 Nat Commun (carrier parents); Arriaga 2025 AJHG (GREGoR/UDN outliers all biallelic) | Population RNA-seq carrier test; a null is expected |

## A third validation set: none is untouched
- RNU4-2 SGE is still the only saturation screen of any snRNA. MaveDB and IGVF hold no snRNA score sets, and the guidance says no functional assays are established for the other snRNAs.
- The best other panel, Benoit-Pilven 2020 (PLoS One, 10.1371/journal.pone.0235655; RNU4ATAC minigene, 24 variants, public S3 Table), is **not untouched**:
  - v4 was trained on RNU4ATAC L2 labels;
  - the assay was scored in the earlier MuST-VEP work;
  - its 16/17 deleterious patient variants vs 7/7 neutral gnomAD variants make it nearly a label set.
- Other candidates:
  - Norppa & Frilander 2021 NAR (RNU12 stem-loop III, ~8 designed variants, stability only): a weak secondary check.
  - Greene 2026 RNU2-2 allelic-imbalance data overlap the blind labels.
  - Episignatures are group-level and controlled-access.

## Must-cite (additions to NOVELTY.md)
- D'Souza 2026 (guidance); De Jonghe 2026 (SGE)
- Messaoud 2026; Chen 2024 Fig. 4; Seplyarskiy 2023
- MobiDeep 2026; MitoTIP 2017; nAPOGEE 2025
- Villani 2024; Walker 2023; Tejura 2024; Tenywa 2025 NARGAB (10.1093/nargab/lqaf157; ncRNA-exonic CADD threshold 11.44, AUC 0.872, the contrast to RNU2-2 CADD 0.50)
- RNU2-2 cohort papers: Greene 2025, Greene 2026, Jackson 2026 (10.1038/s41588-026-02551-9), Chiu 2026 Ann Neurol (10.1002/ana.78071); Rius 2026
- Quinodoz 2026 Nat Genet (10.1038/s41588-025-02451-4; RNU6 adRP)
- Merico 2015; Arriaga 2025; Petrazzini 2024

## Scoop risk
- **High: the Whiffin/Findlay group.** The guidance's Recommendation 6 calls for a calibrated snRNA tool, and the group holds the SGE data and RNUdb.
- **Medium:** MobiDeep's published version (Baux lab).
- **Gene discovery is already scooped** (Messaoud 2026).

## Reviewer risks found in the audit
- L3 = gnomAD AC ≥ 9 at loci with ~50× mutation rate is a weak benign proxy; add a sensitivity analysis at a higher AC.
- Positives overlap the guidance authors' own curation (RNUdb). Label-side circularity is mitigated: no PP3/BP4 was applied to any blind RNU2-2 positive (PREREG s14).
- Supporting-only thresholds give a small expected yield.

## Honest strongest claim
In a pre-registered study with labels sealed before modelling, an snRNA-specific model frozen before unsealing separated RNU2-2 patient variants from population variants (AUROC 0.89) where CADD performed at chance (0.50), and the signal came from RNA-family covariance profiles rather than genome-wide conservation. Its scores give the first ClinGen-calibrated PP3/BP4 thresholds for nuclear snRNA genes, which the 2026 classification guidance identifies as missing. Gene-level constraint, hypermutability and the dominant/recessive depletion contrast replicate and quantify prior observations (Chen 2024; Messaoud 2026; D'Souza 2026) rather than discover them.

**Rating:** novelty ~3.5/5 as is. It could reach ~4.5 if the family-profile result is shown zero-shot and pre-registered across other ncRNA disease genes (RMRP, SNORD118, TERC, RNU7-1, mt-tRNAs).
