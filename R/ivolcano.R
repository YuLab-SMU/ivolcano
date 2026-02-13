#' interactive volcano plot
#'
#' @title ivolcano
#' @inheritParams ivolcano_point
#' @inheritParams geom_ivolcano_gene
#' @param threshold_line customize threshold line style (e.g., line color, type, and width)
#' @return volcano plot
#' @export
#' @examples
#' # example data
#' f <- system.file("extdata/airway.rds", package = "ivolcano")
#' df <- readRDS(f)
#' # plot
#' ivolcano(df,
#'   logFC_col = "log2FoldChange",
#'   pval_col = "padj",
#'   gene_col = "symbol",
#'   onclick_fun = onclick_genecards
#' )
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
  pval_cutoff2 = NULL,
  logFC_cutoff2 = NULL,
  size_by = "none",
  point_size = list(base = 2, medium = 4, large = 6),
  threshold_line = list(color = "black", linetype = "dashed", linewidth = 0.5),
  top_n = 10,
  label_mode = "separate",
  fontface = "italic",
  label_sig_only = TRUE,
  filter = NULL
) {
  p <- ivolcano_point(
    data = data,
    logFC_col = logFC_col,
    pval_col = pval_col,
    gene_col = gene_col,
    title = title,
    interactive = interactive,
    onclick_fun = onclick_fun,
    pval_cutoff = pval_cutoff,
    logFC_cutoff = logFC_cutoff,
    pval_cutoff2 = pval_cutoff2,
    logFC_cutoff2 = logFC_cutoff2,
    size_by = size_by,
    point_size = point_size
  )

  p <- p +
    geom_ivolcano_line(
      color = threshold_line$color %||% "black",
      linetype = threshold_line$linetype %||% "dashed",
      linewidth = threshold_line$linewidth %||% 0.5
    )

  p <- p +
    geom_ivolcano_gene(
      top_n = top_n,
      label_mode = label_mode,
      fontface = fontface,
      label_sig_only = label_sig_only,
      filter = filter
    )
  return(p)
}


#' Visualize points in volcano plot
#' @title ivolcano_point
#' @param data A data frame that contains minimal information with gene id, logFC and adjusted P values
#' @param logFC_col column name in 'data' that stored the logFC values
#' @param pval_col column name in 'data' that stored the adjusted P values
#' @param gene_col column name in 'data' that stored the gene IDs
#' @param title plot title
#' @param interactive whether plot the graph in interactive mode
#' @param onclick_fun effects when click on the dot (gene), default is NULL
#' @param pval_cutoff cutoff of the adjusted P values
#' @param logFC_cutoff cutoff of the logFC values
#' @param pval_cutoff2 second cutoff of the adjusted P values for advanced mode
#' @param logFC_cutoff2 second cutoff of the logFC values for advanced mode
#' @param size_by one of "none" (default), "manual" (set by `point_size`), "negLogP", "absLogFC", or other variable in the input data to scale dot sizes.
#' @param point_size set point size when `size_by` is "manual", a list with three elements: base, medium, large.
#' @return base plot of a volcano plot
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 scale_size_identity
#' @importFrom ggiraph girafe
#' @importFrom ggiraph geom_point_interactive
#' @importFrom ggiraph opts_hover
#' @importFrom ggrepel geom_text_repel
#' @importFrom rlang sym
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom utils modifyList
#' @export
ivolcano_point <- function(
  data,
  logFC_col = "logFC",
  pval_col = "adj.P.Val",
  gene_col = "gene",
  title = "Volcano Plot",
  interactive = TRUE,
  onclick_fun = NULL,
  pval_cutoff = 0.05,
  logFC_cutoff = 1,
  pval_cutoff2 = NULL,
  logFC_cutoff2 = NULL,
  size_by = "none",
  point_size = list(base = 2, medium = 4, large = 6)
) {
  stopifnot(all(c(logFC_col, pval_col, gene_col) %in% colnames(data)))
  size_by <- match.arg(
    size_by,
    c("none", "negLogP", "absLogFC", "manual", colnames(data))
  )

  if (!is.null(pval_cutoff2) & !is.null(logFC_cutoff2)) {
    mode <- "advanced"
  } else {
    mode <- "classic"
  }

  # prepare data
  df <- data[!is.na(data[[logFC_col]]) & !is.na(data[[pval_col]]), ]

  x_range <- range(df[[logFC_col]], na.rm = TRUE)
  xmin_val <- x_range[1]
  xmax_val <- x_range[2]
  # margin, 5% of the range, or at least 1
  margin <- max(1, (xmax_val - xmin_val) * 0.05)
  xlim <- c(xmin_val - margin, xmax_val + margin)

  df <- df |>
    dplyr::mutate(
      negLogP = -log10(!!sym(pval_col)),
      sig = dplyr::case_when(
        !!sym(pval_col) < pval_cutoff & !!sym(logFC_col) > logFC_cutoff ~ "Up",
        !!sym(pval_col) < pval_cutoff &
          !!sym(logFC_col) < -logFC_cutoff ~ "Down",
        TRUE ~ "Not_Significant"
      )
    )

  if (mode == "advanced") {
    df <- df |>
      dplyr::mutate(
        sig = dplyr::case_when(
          !!sym(pval_col) < pval_cutoff2 &
            !!sym(logFC_col) > logFC_cutoff2 ~ "Up2",
          !!sym(pval_col) < pval_cutoff2 &
            !!sym(logFC_col) < -logFC_cutoff2 ~ "Down2",
          .default = !!sym("sig")
        )
      )
  }

  ylab_pval <- "P-value"
  fmt <- "%s: %s\nlogFC: %.3f\nP.val: %.3e"
  if (grepl("adj", tolower(pval_col))) {
    ylab_pval <- "Adjusted P-value"
    fmt <- "%s: %s\nlogFC: %.3f\nadj.P.val: %.3e"
  }

  # Remove single and double quotes from gene names
  if (any(grepl('[\'"]', df[[gene_col]]))) {
    df[[gene_col]] <- gsub("['\"]+", "", df[[gene_col]])
  }

  if (interactive) {
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

    # tooltip
    df$tooltip <- sprintf(
      fmt,
      gene_col,
      df[[gene_col]],
      df[[logFC_col]],
      df[[pval_col]]
    )
  }

  aes_args <- aes(
    x = !!sym(logFC_col),
    y = !!sym("negLogP"),
    color = !!sym("sig")
  )

  if (interactive) {
    aes_args <- modifyList(
      aes_args,
      aes(
        tooltip = .data$tooltip,
        data_id = !!sym(gene_col),
        onclick = !!sym("onclick"),
        `data-plot` = I("volcano"),
        `data-pval` = !!sym(pval_col),
        `data-logfc` = !!sym(logFC_col),
        `data-sig` = !!sym("sig")
      )
    )
  }

  point_params <- list(alpha = 0.6)
  if (interactive) {
    point_params$extra_interactive_params <- c(
      "data-plot",
      "data-pval",
      "data-logfc",
      "data-sig"
    )
  }
  if (size_by != "none") {
    if (size_by == "manual") {
      point_size_cat <- setNames(
        c("large", "large", "medium", "medium", "base"),
        c("Up2", "Down2", "Up", "Down", "Not_Significant")
      )
      point_size <- setNames(point_size[point_size_cat], names(point_size_cat))
      df$point_size <- unlist(point_size[df$sig])
      size_by <- "point_size"
    }

    if (size_by == "absLogFC") {
      df$absLogFC <- abs(df[[logFC_col]])
    }
    point_params$mapping <- aes(size = !!sym(size_by))
  }

  p <- ggplot(df, aes_args) +
    {
      if (interactive) {
        do.call(ggiraph::geom_point_interactive, point_params)
      } else {
        do.call(geom_point, point_params)
      }
    } +
    scale_color_figureya(mode) +
    # labs(title = title, x = "log2 Fold Change", y = "-log10(FDR)") +
    labs(
      x = bquote(~ Log[2] ~ "(Fold Change)"),
      y = bquote(~ -Log[10] ~ "(" ~ italic(.(ylab_pval)) ~ ")"),
      title = title
    ) +
    coord_cartesian(xlim = xlim) +
    theme_minimal()

  if (size_by == "point_size") {
    p <- p +
      scale_size_identity()
  }

  p@plot_env$interactive <- interactive

  class(p) <- c("ivolcano", class(p))

  #ivolcano_toplot(p)
  return(p)
}
