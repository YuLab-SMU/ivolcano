#' @method print ivolcano
#' @export
print.ivolcano <- function(x, ...) {
  class(x) <- class(x)[-1]
  if (x@plot_env$interactive) {
    girafe(ggobj = x, options = list(opts_hover(css = "fill:black;r:6")))
  } else {
    print(x)
  }
}
