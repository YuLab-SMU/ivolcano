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
#' @rdname ivolcano
#' @export
geom_ivolcano_gene <- function(
  top_n = 10,
  label_mode = "separate",
  fontface = "italic",
  label_sig_only = TRUE,
  filter = NULL
) {
  structure(
    list(
      top_n = top_n,
      label_mode = label_mode,
      fontface = fontface,
      label_sig_only = label_sig_only,
      filter = filter
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

  df_label <- df
  if (label_sig_only) {
    df_label <- dplyr::filter(df, .data$sig != "Not_Significant")
  }

  if (!is.null(filter)) {
    filter_expr <- rlang::parse_expr(filter)
    df_label <- dplyr::filter(df_label, !!filter_expr)
  } else if (top_n > 0) {
    # label topN genes

    if (label_mode == "separate") {
      df_label <- dplyr::group_by(df_label, sign(!!sym(logFC_col)))
    }

    df_label <- dplyr::slice_min(
      df_label,
      order_by = !!sym(pval_col),
      n = top_n
    )
  }

  df_label$label <- df_label[[gene_col]]

  obj <- geom_text_repel(
    data = df_label,
    aes(label = !!sym("label"), fontface = fontface),
    box.padding = 0.3,
    max.overlaps = 20,
    bg.colour = "white",
    bg.r = 0.15
  )

  ggplot_add(obj, plot, ...)
}
