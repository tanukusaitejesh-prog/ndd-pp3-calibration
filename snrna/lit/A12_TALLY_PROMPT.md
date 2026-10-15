# Prompt for the A12 literature tally (paste into a research AI; save the TSV as lit/a12_<GENE-or-batch>_<YYYY-MM-DD>.tsv)

TASK: Literature tally of patient variants in human spliceosomal snRNA genes, for a pre-registered prospective test
(amendment A12). The model being tested is already frozen; your job is only to record what papers report.
GENES, in priority order:
1. RNU6-1, RNU6-2, RNU6-7, RNU6-8, RNU6-9 (autosomal dominant retinitis pigmentosa; e.g. Quinodoz et al. 2026 Nat Genet)
2. RNU5A-1, RNU5B-1 (neurodevelopmental disorder)
3. RNU2-2 (only variants first reported on or after 2026-10-14; earlier ones are already in the study)
4. RNU11, RNU4-1, RNU5D-1, RNU5E-1, RNU5F-1, RNVU1-8, RNVU1-32, RNVU1-33, RNVU1-34
5. RNU4-2, RNU4ATAC, RNU12, RNU6ATAC (only variants first reported on or after 2026-10-14)
RULES:
- Record every distinct variant reported (peer-reviewed or preprint, incl. supplementary tables, OMIM, LOVD, RNUdb) as
  pathogenic, likely pathogenic or disease-causing in affected individuals. Never invent or infer a variant.
- Copy the classification wording exactly as stated (classification_as_stated). If a source calls it a VUS, write that.
- Do NOT look up or record any computational score (CADD, phyloP, conservation, predictors) for any variant.
OUTPUT: one TSV, one row per variant per source, header exactly:
gene	transcript	hgvs_n	genomic_grch38	genomic_grch37	variant_type	classification_as_stated	acmg_classification_given	inheritance	zygosity_observed	phenotype	n_families	source	first_publication_year	first_publication_date	insilico_or_conservation_used_for_classification	notes
(hgvs_n like n.51G>A on the gene's own transcript; source = PMID or DOI of the EARLIEST report, others in notes;
inheritance dominant/recessive; insilico_or_conservation_used_for_classification yes/no/unknown.)
Then report counts per gene, every source checked (with date), uncertain items, and whether any gene had no report.
NEXT (analyst): python a12_new_variants.py map lit/a12_*.tsv && python a12_new_variants.py build && commit.
