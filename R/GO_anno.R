#' Get GO annotaion for each gene
#' 
#' @param genes genes to annotation.
#' @param orgdb_pkg package for GO annotation.
#' @param keytype type for gene name.
#' @param n_each top n item to show.
get_gene_go_summary_batch <- function(genes,
                                      orgdb_pkg = "org.Hs.eg.db",
                                      keytype = "SYMBOL",
                                      n_each = 2) {
  rlang::check_installed(
    "AnnotationDbi",
    reason = "For retrieving gene-to-GO mappings."
  )
  rlang::check_installed(
    "GO.db",
    reason = "For retrieving GO term annotation."
  )
  rlang::check_installed(
    orgdb_pkg,
    reason = "For retrieving organism-specific gene annotation."
  )

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0) {
    return(data.frame(
      gene = character(),
      GO_id = character(),
      GO_term_detail = character(),
      ONTOLOGY = character(),
      stringsAsFactors = FALSE
    ))
  }

  orgdb <- NULL
  go_db <- NULL

  if (requireNamespace(orgdb_pkg, quietly = TRUE)) {
    orgdb <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))
  }

  if (requireNamespace("GO.db", quietly = TRUE)) {
    go_db <- get("GO.db", envir = asNamespace("GO.db"))
  }

  gene_go <- suppressMessages(
    AnnotationDbi::select(
      orgdb,
      keys = genes,
      keytype = keytype,
      columns = c(keytype, "GO", "ONTOLOGY")
    )
  )

  gene_go <- dplyr::distinct(gene_go)
  gene_go <- dplyr::filter(gene_go, !is.na(.data$GO), !is.na(.data$ONTOLOGY))

  if (nrow(gene_go) == 0) {
    return(data.frame(
      gene = character(),
      GO_id = character(),
      GO_term_detail = character(),
      ONTOLOGY = character(),
      stringsAsFactors = FALSE
    ))
  }

  key_col <- keytype

  go_detail <- suppressMessages(
    AnnotationDbi::select(
      go_db,
      keys = unique(gene_go$GO),
      keytype = "GOID",
      columns = c("GOID", "TERM", "ONTOLOGY")
    )
  )

  go_detail <- dplyr::distinct(go_detail, .data$GOID, .data$TERM, .data$ONTOLOGY)

  result <- dplyr::inner_join(
    gene_go,
    go_detail,
    by = c("GO" = "GOID", "ONTOLOGY" = "ONTOLOGY")
  )

  result <- dplyr::transmute(
    result,
    gene = .data[[key_col]],
    GO_id = .data$GO,
    GO_term_detail = .data$TERM,
    ONTOLOGY = .data$ONTOLOGY
  )

  result <- dplyr::distinct(result)

  result_top <- result |>
    dplyr::group_by(.data$gene, .data$ONTOLOGY) |>
    dplyr::slice_head(n = n_each) |>
    dplyr::ungroup()

  result_top
}


#' Format GO annotaion for each gene
#' 
#' @param genes genes to annotation.
#' @param orgdb_pkg package for GO annotation.
#' @param keytype type for gene name.
#' @param n_each top n item to show.
format_gene_go_tooltip_batch <- function(genes,
                                         orgdb_pkg = "org.Hs.eg.db",
                                         keytype = "SYMBOL",
                                         n_each = 2) {
  go_df <- get_gene_go_summary_batch(
    genes = genes,
    orgdb_pkg = orgdb_pkg,
    keytype = keytype,
    n_each = n_each
  )

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  if (length(genes) == 0) {
    return(stats::setNames(character(), character()))
  }

  if (nrow(go_df) == 0) {
    return(stats::setNames(rep("", length(genes)), genes))
  }

  go_txt_df <- go_df |>
    dplyr::group_by(.data$gene, .data$ONTOLOGY) |>
    dplyr::summarise(
      txt = paste(.data$GO_term_detail, collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      txt = paste0("GO ", .data$ONTOLOGY, ": ", .data$txt)
    ) |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarise(
      tooltip_go = paste(.data$txt, collapse = "\n"),
      .groups = "drop"
    )

  out <- stats::setNames(rep("", length(genes)), genes)
  idx <- match(go_txt_df$gene, names(out))
  out[idx] <- go_txt_df$tooltip_go
  out
}