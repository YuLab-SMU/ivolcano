#' Add threshold lines to an iVolcano plot
#'
#' @title geom_ivolcano_line
#' @param linetype line type for the threshold lines
#' @param color line color for the threshold lines
#' @param linewidth line width for the threshold lines
#' @return A ggplot2 layer object
#' @export
geom_ivolcano_line <- function(
  linetype = "longdash",
  color = 'grey40',
  linewidth = 0.5
) {
  structure(
    list(
      linetype = linetype,
      color = color,
      linewidth = linewidth
    ),
    class = "ivolcano_line"
  )
}

#' @importFrom ggplot2 ggplot_add
#' @importFrom ggplot2 theme
#' @method ggplot_add ivolcano_line
#' @export
ggplot_add.ivolcano_line <- function(object, plot, ...) {
  xintercept <- c(-plot@plot_env$logFC_cutoff, plot@plot_env$logFC_cutoff)
  yintercept <- -log10(plot@plot_env$pval_cutoff)

  obj <- add_ivocano_line(
    xintercept = xintercept,
    yintercept = yintercept,
    linetype = object$linetype,
    color = object$color,
    linewidth = object$linewidth
  )

  if (plot@plot_env$mode == "advanced") {
    xintercept2 <- c(-plot@plot_env$logFC_cutoff2, plot@plot_env$logFC_cutoff2)
    yintercept2 <- -log10(plot@plot_env$pval_cutoff2)
    obj <- c(
      obj,
      add_ivocano_line(
        xintercept = xintercept2,
        yintercept = yintercept2,
        linetype = object$linetype,
        color = object$color,
        linewidth = object$linewidth
      )
    )
  }

  obj <- list(obj, theme(panel.grid = ggplot2::element_blank()))
  ggplot_add(obj, plot, ...)
}


add_ivocano_line <- function(
  xintercept,
  yintercept,
  linetype = "longdash",
  color = "grey40",
  linewidth = 0.5
) {
  list(
    geom_hline(
      yintercept = yintercept,
      color = color,
      linetype = linetype,
      linewidth = linewidth
    ),
    geom_vline(
      xintercept = xintercept,
      color = color,
      linetype = linetype,
      linewidth = linewidth
    )
  )
}
