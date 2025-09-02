#' interactive volcano plot
#' 
#' @title ivolcano
#' @param data A data frame that contains minimal information with gene id, logFC and adjusted P values
#' @param logFC_col column name in 'data' that stored the logFC values
#' @param pval_col column name in 'data' that stored the adjusted P values
#' @param gene_col column name in 'data' that stored the gene IDs
#' @param title plot title
#' @param interactive whether plot the graph in interactive mode
#' @param onclick_fun effects when click on the dot (gene), default is NULL
#' @param pval_cutoff cutoff of the adjusted P values
#' @param logFC_cutoff cutoff of the logFC values
#' @param top_n top N genes to display the labels (gene ID)
#' @param label_mode one of 'all' or 'separate' (default). 
#'    If label_mode = 'all', top_n genes with minimal p values will be displayed, 
#'    otherwise, top_n up-regulated and top_n down-regulated genes will be displayed.
#' @param fontface one of 'plain', 'bold', 'italic' (default) and their combination, e.g. 'bold.italic'
#' @param label_sig_only whether filter significant genes before subset 'top_n' genes
#' @param threshold_line customize threshold line style (e.g., line color, type, and width)
#' @param sig_colors customize colors for up- and down-regulated, and non-significant genes
#' @param size_by one of "none" (default), "negLogP", or "absLogFC" to scale dot sizes.
#' @return volcano plot
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggiraph girafe
#' @importFrom ggiraph geom_point_interactive
#' @importFrom ggiraph opts_hover
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang sym
#' @importFrom rlang .data
#' @export
#' @examples
#' # example data
#' f <- system.file('extdata/airway.rds', package='ivolcano')
#' df <- readRDS(f)
#' # plot
#' ivolcano(df,
#'         logFC_col = "log2FoldChange",
#'         pval_col = "padj",
#'         gene_col = "symbol",
#'         top_n = 5,
#'         onclick_fun=onclick_genecards)
#' @author Guangchuang Yu
ivolcano <- function(
  data,
  logFC_col = "logFC",
  pval_col = "adj.P.Val",
  gene_col = "gene",
  title = "Volcano Plot",
  interactive = TRUE,
  onclick_fun = NULL,
  pval_cutoff = 0.05,
  logFC_cutoff = 1,
  top_n = 10,
  label_mode = "separate",
  fontface = "italic",
  label_sig_only = TRUE,
  threshold_line = list(color = "black", linetype = "dashed", linewidth = 0.5),
  sig_colors = c(Up = "red", Down = "blue", Not_Significant = "grey70"),
  size_by = "none"
) {
  stopifnot(all(c(logFC_col, pval_col, gene_col) %in% colnames(df)))
  label_mode <- match.arg(label_mode, c("all", "separate"))
  size_by <- match.arg(size_by, c("none", "negLogP", "absLogFC"))

  # prepare data
  df <- data[!is.na(data[[logFC_col]]) & !is.na(data[[pval_col]]), ]

  df <- df |>
    dplyr::mutate(
      negLogP = -log10(!!sym(pval_col)),
      sig = dplyr::case_when(
        !!sym(pval_col) < pval_cutoff & !!sym(logFC_col) > logFC_cutoff ~ "Up",
        !!sym(pval_col) < pval_cutoff & !!sym(logFC_col) < -logFC_cutoff ~
          "Down",
        TRUE ~ "Not_Significant"
      )
    )

  df$tooltip <- sprintf(
    "Gene: %s\nlogFC: %s\nadj.P.val: %s",
    df[[gene_col]],
    signif(df[[logFC_col]], 3),
    signif(df[[pval_col]], 3)
  )

  # onclick
  if (!is.null(onclick_fun)) {
    if (is.function(onclick_fun)) {
      df$onclick <- vapply(
        df[[gene_col]],
        function(g) {
          x <- onclick_fun(g)
          if (is.na(x) || !nzchar(x)) "" else as.character(x)
        },
        character(1)
      )
    } else {
      stop("onclick_fun must be a function or NULL")
    }
  } else {
    df$onclick <- ""
  }
  aes_args <- aes(
    x = !!sym(logFC_col),
    y = !!sym("negLogP"),
    color = !!sym("sig"),
    tooltip = paste0(
      gene_col,
      ": ",
      !!sym(gene_col),
      "\nlogFC: ",
      round(!!sym(logFC_col), 2),
      "\nFDR: ",
      signif(!!sym(pval_col), 3)
    ),
    data_id = !!sym(gene_col),
    onclick = !!sym("onclick")
  )
  if (!is.null(size_by) && size_by %in% names(df)) {
    aes_args$size <- df[[size_by]]
  }

  p <- ggplot(df, aes_args) +
    {
      if (interactive) {
        ggiraph::geom_point_interactive(alpha = 0.7)
      } else {
        geom_point(alpha = 0.7)
      }
    } +
    scale_color_manual(values = sig_colors) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(FDR)") +
    theme_minimal()

  # threshold line
  p <- p +
    geom_hline(
      yintercept = -log10(pval_cutoff),
      color = threshold_line$color,
      linetype = threshold_line$linetype,
      linewidth = threshold_line$linewidth
    ) +
    geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      color = threshold_line$color,
      linetype = threshold_line$linetype,
      linewidth = threshold_line$linewidth
    )

  # label topN genes
  if (top_n > 0) {
    df_label <- df
    if (label_sig_only) {
      df_label <- dplyr::filter(df, .data$sig != "Not_Significant")
    }
    if (label_mode == "separate") {
      df_label <- dplyr::group_by(df_label, sign(!!sym(logFC_col)))
    }
    df_label <- dplyr::slice_min(
      df_label,
      order_by = !!sym(pval_col),
      n = top_n
    )

    df_label$label <- df_label[[gene_col]]

    p <- p +
      geom_text_repel(
        data = df_label,
        aes(label = !!sym("label"), fontface = fontface),
        box.padding = 0.3,
        max.overlaps = 20
      )
  }

  if (interactive) {
    girafe(ggobj = p, options = list(opts_hover(css = "fill:black;r:6")))
  } else {
    p
  }
}

