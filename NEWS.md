# ivolcano 0.0.6

+ `ivolcano()` and `ivolcano_point()` gain `tooltip_cols` and `tooltip_labels` parameters,
  allowing users to display arbitrary columns from the input data frame in the interactive
  tooltip (e.g., gene biotype, GO terms, functional descriptions joined from external databases). (2026-06-04, Thu)


# ivolcano 0.0.5

+ `pathway_volcano()` now highlights pathway genes in the volcano plot after selection,
  while keeping all points and graying out non-pathway genes.
+ After pathway selection, recompute and render labels within the pathway subset
  following `geom_ivolcano_gene()` rules (top_n/label_mode/label_sig_only).
+ Fix an issue where volcano points/labels could be unintentionally hidden after selection
  by caching and restoring original SVG styles.
+ `idotplot()` removes deprecated ggplot2 usage (replace `aes_string()` with tidy-eval).
+ `ivolcano()` adds `data-pval`/`data-logfc`/`data-sig` attributes to interactive points
  for front-end linkage.
+ `geom_ivolcano_gene()` writes label settings into plot_env for `pathway_volcano()`.

# ivolcano 0.0.4

+ explicitly call `knitr::knit_print` when running in knitr (2025-12-02, Tue)
+ Add more `onclick` functions (2025-12-01, Mon)
    - onclick_ncbi to open NCBI Gene database.
    - onclick_ensembl to open Ensembl database.
    - onclick_hgnc to open HGNC symbol report.
    - onclick_uniprot to search UniProt.
    - onclick_pubmed to search PubMed for the gene.

# ivolcano 0.0.3

+ geom_ivolcano_line, geom_ivolcano_gene, and scale_color_figuareya/scale_fill_figureya (2025-11-21, Fri)

# ivolcano 0.0.2

+ fixed typo

# ivolcano 0.0.1

+ ivolcano
+ onclick_fanyi
+ onclick_genecards
