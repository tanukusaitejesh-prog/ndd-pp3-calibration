# Novelty audit of the rarity-inflation claims (4 Nov 2026)

This is a web literature search. Several bioRxiv full texts could not be opened (HTTP 429), so some entries rest on abstracts or snippets.

**Verdict:** the central idea is a known caveat. What the manuscript adds is a new quantification and a new correction.

| Claim | Verdict | Closest prior work | What remains new |
|---|---|---|---|
| R1 Common negatives inflate rarity-tracking scores | ALREADY DONE as a caveat | Niroula & Vihinen 2019 PLoS Comput Biol 10.1371/journal.pcbi.1006481; Rastogi...Ioannidis 2025 Hum Genet 10.1007/s00439-025-02732-2 (AF-matched subsets, AF bins; closest); Grimm 2015 Hum Mutat 10.1002/humu.22768; Gunning 2021 J Med Genet PMID 32843488 | Stating it for non-coding conservation and DNA-model scores |
| R2 Rarity-only null for each score; excess AUROC | PARTLY DONE | Rastogi 2025 (AF matching); Caron 2019 Genome Biol NCBoost 10.1186/s13059-019-1634-2 (rare-negative sensitivity analysis); Lu, Lin & Brandes 2025 bioRxiv 10.1101/2025.09.05.674459 (context-prevalence null) | Drawing random rare unlabelled variants as a null for each score, and reporting the excess: no exact precedent found |
| R3 ClinVar non-coding: rarity = 40-58% of above-chance AUROC | NEW (as a number) | - | The numbers |
| R4 AlphaMissense 0.92 -> 0.69; lead over CADD vanishes | PARTLY DONE; CONFLICT | Tabet 2024 Genome Biol 10.1186/s13059-024-03314-7: rare-only UKB trait benchmark where AlphaMissense ranks first | Must be reconciled with Tabet (different endpoint) |
| R5 NCBoost and TraitGym re-analysis; rankings change | NEW with caveats | TraitGym Mendelian controls are MAF > 5% and not matched (Benegas 2025 bioRxiv 10.1101/2025.02.11.637758); Caron 2019 had a partial rare-negative check | TraitGym Mendelian inflation and ranking changes |
| R6 Simulation; deattenuated rarity-matched AUROC | NEW (method) | - | The estimator and the simulation |
| R7 PP3/BP4 biased by common benign variants | MIS-PREMISED as written | Pejaver 2022 AJHG (PMC9748256) already used AF < 0.01 benign sets | Only residual bias within AF < 0.01 (e.g. singletons vs 0.1-1%) could be new; reframe |

**Must cite:** Rastogi 2025; Niroula & Vihinen 2019; Grimm 2015; Gunning 2021; Caron 2019; Benegas 2025; Lu/Lin/Brandes 2025; Pejaver 2022; Tejura 2024 (AF content unverified); Tabet 2024; Livesey & Marsh 2023 Mol Syst Biol 10.15252/msb.202211474; Cheng 2023 Science 10.1126/science.adg7492.

**Not found:** "Zhang 2026 hash-frozen benchmark". Do not cite it unless it is verified.

**Actions for the manuscript:**
- Reframe the rarity result as a correction that quantifies a known bias, not a discovery; cite Rastogi 2025 as the closest prior work.
- Rewrite R7 as residual bias inside AF < 0.01.
- Reconcile R4 with Tabet 2024.
