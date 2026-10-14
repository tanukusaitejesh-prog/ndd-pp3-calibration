# UK Biobank Allele Frequency Browser exports (v3 input)

Exported by the authors from the UK Biobank WGS Allele Frequency Browser (about 490,640 genomes), one CSV per region in
`../../ukb_regions.txt`: all VEP impact classes, no column filters. Columns: Variant ID (GRCh38 chrom-pos-ref-alt),
rsID, Gene, VEP Annotation, HGVS Consequence, Allele Count, Allele Number, Allele Frequency, Number of Homozygotes.
Export date: 2026-10-14. Use is subject to the browser's terms; check them before reuse or redistribution.

Build the v3 inputs: `python v3_counts.py <(...)` with all CSVs combined, see README "v3".
