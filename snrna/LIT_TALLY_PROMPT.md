# Prompt for the literature tally (paste into a research AI; output goes to curate_literature.py)

TASK: Literature tally of patient variants in human spliceosomal snRNA genes, for a pre-registered genetics study.
GENES (only these): RNU4ATAC, RNU12, RNU6ATAC, RNU11, RNU4-1, RNU4-2, RNU5A-1, RNU5D-1, RNU5E-1, RNU5F-1, RNU6-1, RNU6-2,
RNU6-7, RNU6-8, RNU6-9, RNVU1-8, RNVU1-32, RNVU1-33, RNVU1-34. Priority: RNU4ATAC, RNU12, RNU6 paralogs, then the rest.
STRICT BLINDING RULE: never search for, open, read, record or mention any variant in RNU2-2 or RNU5B-1.
WHAT COUNTS: every distinct variant reported (peer-reviewed or preprint, incl. supplementary tables, OMIM, LOVD) as
pathogenic, likely pathogenic or disease-causing in affected individuals. Never invent a variant.
OUTPUT: one TSV, one row per variant, header exactly:
gene	transcript	hgvs_n	genomic_grch38	genomic_grch37	variant_type	classification_as_stated	acmg_classification_given	inheritance	zygosity_observed	phenotype	n_families	source	first_publication_year	first_publication_date	insilico_or_conservation_used_for_classification	notes
(hgvs_n like n.51G>A; source = PMID of the EARLIEST report, others in notes; inheritance dominant/recessive.)
Then report counts per gene, sources checked, uncertain items, and confirm nothing about RNU2-2/RNU5B-1 was recorded.
