# Zenodo data audit

Record: 10.5281/zenodo.18221260; licence cc-by-4.0. Version relation is saved in manifest.json.
Inventory: 15 files, 170,568,661,658 compressed bytes.
Downloaded only the two prespecified small archives; both provided checksums verified.
Full SHA256 values, URLs and file sizes are recorded in manifest.json.

Counts: {'gene': 403011, 'susie': 4193682}. All 56 table members were streamed through to EOF.
Member names, columns, missing-value denominators, first rows and PIP counts are in archive_profiles.json.
Delimiter is detected from contents; .tsv members use commas.

Gene-level schema: gene, ACAT_p, top_MarkerID, top_pval. No beta, SE, effect allele, frequency or donor count.
SuSiE schema: SNP, Index, PIP, Credible_Set, gene, celltype. No beta, SE or frequency.
Variant IDs parse as chromosome:one-based-position:allele:allele; REF/ALT meanings still need reference validation.
No oriented observed betas can be constructed from these two archives.
Missing and convergence-failure rows are counted rather than silently discarded.

Observed cells: ASDC, B_intermediate, B_memory, B_naive, CD14_Mono, CD16_Mono, CD4_CTL, CD4_Naive, CD4_Proliferating, CD4_TCM, CD4_TEM, CD8_Naive, CD8_Proliferating, CD8_TCM, CD8_TEM, HSPC, ILC, MAIT, NK, NK_CD56bright, NK_Proliferating, Plasmablast, Treg, cDC1, cDC2, dnT, gdT, pDC.
Explicit proposed cell mapping: {'B naive': 'B_naive', 'CD4 naive T': 'CD4_Naive', 'CD14 classical monocyte': 'CD14_Mono', 'CD16 nonclassical monocyte': 'CD16_Mono'}; labels validated during selection.
SuSiE credible sets are 95% sets per Zenodo metadata, not certain causal labels.
No donor-level data are requested or inferred. Sample counts and normalization are not recoverable from these schemas.

Official project source inspection finds GENCODE v44 basic annotation paths in
Scanpy/combine_files_add_gene_info.py; this supports release choice but does not establish one unique canonical TSS per gene.
GRCh38, exact reference alleles, beta allele convention and assay normalization require further verification against full summaries and project code before efficacy scoring.
Full 16.6-GB summary archive deferred because inference credentials are missing; download_large_files=false.
The preprint linked by the record is https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v2 .
Its sample counts have not been independently verified here.
