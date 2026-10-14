You are checking the novelty of a research paper's claims against the published and preprint literature. Use web
search (PubMed, Europe PMC, bioRxiv/medRxiv, Google Scholar). Be skeptical: your job is to find prior work that
already did each thing, not to praise the paper. Today is October 2026; include 2024-2026 work and preprints.

## Context
The paper is about spliceosomal small nuclear RNA (snRNA) genes in neurodevelopmental disorders (NDD): RNU4-2
(ReNU syndrome, 2024), RNU2-2, RNU5B-1, RNU4ATAC (MOPD1/Roifman/Lowry-Wood), RNU12, RNU6ATAC, and RNU6-1/2/8/9
(retinitis pigmentosa). It uses gnomAD v4 and UK Biobank allele counts, Evo 2, CADD, ViennaRNA, RNU4-2 saturation
genome editing (SGE) data, and ClinGen-style (PP3/BP4) calibration.

## Claims to check (for each: has it been done? by whom? how close?)
C1. Pre-registered, timestamped (OpenTimestamps) evaluation of a variant-effect predictor for snRNA genes, with a
    blind test on sealed patient labels (RNU2-2, RNU5B-1).
C2. Transfer of RNU4-2 SGE functional data to its minor-spliceosome homolog RNU4ATAC (via structural alignment),
    used to score RNU4ATAC variants; this beat CADD (AUROC 0.85 vs 0.77) for recessive disease.
C3. A framework score that uses population selection for dominant genes and lab-derived function scores for
    recessive genes (inheritance-aware scoring), beating CADD on blind genes (delta AUROC +0.27).
C4. A systematic scan of all ~55 snRNA loci and paralogs/pseudogenes for local (within-gene, 20-nt window)
    depletion of variation in gnomAD + UK Biobank, recovering all five snRNA NDD genes at the top and the known
    disease regions (e.g. RNU4-2 n.62-79); plus the observation that whole-gene comparisons fail because expressed
    snRNA genes are germline mutation hotspots (transcription-associated mutagenesis).
C5. Predicting inheritance mode (dominant vs recessive action) of pathogenic noncoding/snRNA variants from local
    population constraint (dominant variants depleted in heterozygotes; recessive variants are not).
C6. A stacked, regularised classifier (CADD, GPN-MSA, phyloP, ReMM, Evo 2, ViennaRNA ddG, Rfam profile, local
    constraint, homology transfer) for snRNA variants, validated leave-one-paralog-cluster-out and blind.
C7. ClinGen-calibrated PP3/BP4 thresholds specific to snRNA genes.
C8. Healthy-carrier readout (1000 Genomes + Geuvadis junctions) of minor-intron splicing for snRNA variants.
C9. A released score table for every possible SNV in ~55 snRNA loci.

## Search hints
snRNA constraint gnomAD; RNU4-2 depletion population; ReNU syndrome; RNU2-2 neurodevelopmental; RNU5B-1;
RNU4ATAC variant interpretation; snRNA variant effect predictor; RIBOSCOPE; noncoding constraint Gnocchi;
transcription-associated mutagenesis snRNA hotspot; U snRNA germline mutation rate; saturation genome editing
RNU4-2; minor intron retention RNU4ATAC carriers; predicting dominant recessive inheritance variants; ClinGen
snRNA guidance; PP3 calibration noncoding.

## Output
1. A table: claim | verdict (NEW / PARTLY DONE / ALREADY DONE) | closest prior work (authors, year, journal,
   PMID/DOI) | what exactly overlaps | what remains new.
2. A short list of papers we MUST cite.
3. One paragraph: the honest strongest novelty claim the paper can make.
Cite only papers you actually found; mark anything uncertain as UNVERIFIED. Do not invent references.
