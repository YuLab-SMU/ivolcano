#' scale_color_figureya
#'
#' @title scale_color_figureya
#' @param mode one of "classic" or "advanced"
#' @param na.value color for NA values
#' @param guide guide parameter passed to ggplot2::scale_color_manual
#' @param ... additional parameters passed to ggplot2::scale_color_manual
#' @return A ggplot2 scale object
#' @importFrom ggplot2 scale_color_manual
#' @export
scale_color_figureya <- function(
  mode = c("classic", "advanced"),
  na.value = "grey80",
  guide = "none",
  ...
) {
  mode <- match.arg(mode)
  cols <- figureya_pal(mode)
  ggplot2::scale_color_manual(
    values = cols,
    na.value = na.value,
    guide = guide,
    ...
  )
}


figureya_pal <- function(mode = c("classic", "advanced")) {
  mode <- match.arg(mode)
  default_palettes <- list(
    classic = c(Up = "red", Down = "blue", Not_Significant = "grey"),
    advanced = c(
      Up = "#FB9A99",
      Up2 = "#ED4F4F",
      Down = "#B2DF8A",
      Down2 = "#329E3F",
      Not_Significant = "grey"
    )
  )
  base <- default_palettes[[mode]]
  return(base)
}
