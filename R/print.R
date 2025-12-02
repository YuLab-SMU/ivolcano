#' @importFrom ggiraph girafe
#' @importFrom ggiraph opts_hover
#' @importFrom knitr knit_print
#' @export
print.ivolcano <- function(x, ...) {
  x <- ivolcano_toplot(x)

  # If running in knitr (Quarto), delegate to knit_print
  if (isTRUE(getOption('knitr.in.progress'))) {
    return(knitr::knit_print(x, ...))
  }

  print(x, ...)
}



ivolcano_toplot <- function(x) {
  # use NextMethod() or strip class to avoid infinite recursion
  class(x) <- setdiff(class(x), "ivolcano")

  interactive <- x@plot_env$interactive
  if (is.null(interactive)) {
    interactive <- FALSE
  }

  if (interactive) {
    opts_list <- list(opts_hover(css = "fill:black;r:6"))

    #if (x@plot_env$toolbar) {
    #  opts_list <- c(opts_list, list(
    #    opts_toolbar(saveaspng = TRUE)
    #  ))
    #}

    x <- girafe(ggobj = x, options = opts_list)
  } 

  return(x)
}
