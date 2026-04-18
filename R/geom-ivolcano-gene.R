#' Add gene labels to ivolcano plot
#'
#' @title geom_ivolcano_gene
#' @param top_n top N genes to display the labels (gene ID)
#' @param label_mode one of 'all' or 'separate' (default).
#'    If label_mode = 'all', top_n genes with minimal p values will be displayed,
#'    otherwise, top_n up-regulated and top_n down-regulated genes will be displayed.
#' @param fontface one of 'plain', 'bold', 'italic' (default) and their combination, e.g. 'bold.italic'
#' @param label_sig_only whether filter significant genes before subset 'top_n' genes
#' @param filter custom filter expression to select genes for labeling
#' @param orgdb_pkg package for GO annotation.
#' @param orgdb_keytype type for gene name.
#' @param go_n_each top n GO item to show.
#' @rdname ivolcano
#' @export
geom_ivolcano_gene <- function(
  top_n = 10,
  label_mode = "separate",
  fontface = "italic",
  label_sig_only = TRUE,
  filter = NULL,
  orgdb_pkg = NULL,
  orgdb_keytype = "SYMBOL",
  go_n_each = 2
) {
  structure(
    list(
      top_n = top_n,
      label_mode = label_mode,
      fontface = fontface,
      label_sig_only = label_sig_only,
      filter = filter,
      orgdb_pkg = orgdb_pkg,
      orgdb_keytype = orgdb_keytype,
      go_n_each = go_n_each
    ),
    class = "ivolcano_gene"
  )
}

#' @importFrom ggplot2 ggplot_add
#' @method ggplot_add ivolcano_gene
#' @export
ggplot_add.ivolcano_gene <- function(object, plot, ...) {
  df <- plot@data
  gene_col <- plot@plot_env$gene_col
  logFC_col <- plot@plot_env$logFC_col
  pval_col <- plot@plot_env$pval_col

  top_n <- object$top_n
  label_mode <- object$label_mode
  fontface <- object$fontface
  label_sig_only <- object$label_sig_only
  filter <- object$filter
  orgdb_pkg <- object$orgdb_pkg
  orgdb_keytype <- object$orgdb_keytype
  go_n_each <- object$go_n_each

  plot@plot_env$label_top_n <- top_n
  plot@plot_env$label_mode <- label_mode
  plot@plot_env$label_sig_only <- label_sig_only
  plot@plot_env$label_filter <- filter
  plot@plot_env$label_fontface <- fontface
  plot@plot_env$orgdb_pkg <- orgdb_pkg
  plot@plot_env$orgdb_keytype <- orgdb_keytype
  plot@plot_env$go_n_each <- go_n_each

  if (!is.null(orgdb_pkg) && isTRUE(plot@plot_env$interactive)) {
    genes <- as.character(df[[gene_col]])

    go_txt_map <- tryCatch(
      format_gene_go_tooltip_batch(
        genes = genes,
        orgdb_pkg = orgdb_pkg,
        keytype = orgdb_keytype,
        n_each = go_n_each
      ),
      error = function(e) {
        stats::setNames(rep("", length(unique(genes))), unique(genes))
      }
    )

    go_txt <- unname(go_txt_map[genes])
    go_txt[is.na(go_txt)] <- ""

    if (!"tooltip" %in% colnames(df)) {
      df$tooltip <- as.character(df[[gene_col]])
    }

    df$tooltip <- ifelse(
      nzchar(go_txt),
      paste0(df$tooltip, "\n", go_txt),
      df$tooltip
    )

    plot@data <- df
  }

  df_label <- df
  if (label_sig_only) {
    df_label <- dplyr::filter(df_label, .data$sig != "Not_Significant")
  }

  if (!is.null(filter)) {
    filter_expr <- rlang::parse_expr(filter)
    df_label <- dplyr::filter(df_label, !!filter_expr)
  } else if (top_n > 0) {
    if (label_mode == "separate") {
      df_label <- dplyr::group_by(df_label, sign(!!rlang::sym(logFC_col)))
    }

    df_label <- dplyr::slice_min(
      df_label,
      order_by = !!rlang::sym(pval_col),
      n = top_n
    )
  }

  df_label$label <- df_label[[gene_col]]

  obj <- ggrepel::geom_text_repel(
    data = df_label,
    ggplot2::aes(label = !!rlang::sym("label")),
    fontface = fontface,
    box.padding = 0.3,
    max.overlaps = 20,
    bg.colour = "white",
    bg.r = 0.15
  )

  ggplot2::ggplot_add(obj, plot, ...)
}
