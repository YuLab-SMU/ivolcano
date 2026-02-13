#' Interactive Dot Plot
#'
#' @title idotplot
#' @param object enrichment result object (e.g. enrichResult) or a data frame
#' @param x x-axis variable, one of 'GeneRatio', 'Count', etc.
#' @param color variable for color, e.g., 'p.adjust'
#' @param showCategory number of categories to show
#' @param size variable for size, e.g., 'Count'
#' @param label_format a numeric value sets wrap length, alternatively a custom function to format axis labels.
#' @param trigger interaction trigger, one of "click" or "hover"
#' @param title plot title
#' @param ... additional parameters passed to fortify
#' @return ggplot object with interactive layers
#' @importFrom ggplot2 ggplot aes scale_color_continuous theme_minimal ylab labs fortify
#' @importFrom ggiraph geom_point_interactive
#' @importFrom rlang .data sym
#' @importFrom stats reorder
#' @importFrom methods is
#' @export
#' @author Guangchuang Yu
idotplot <- function(object, x = "GeneRatio", color = "p.adjust", showCategory = 10, size = "Count", label_format = 30, trigger = c("click", "hover"), title = "", ...) {
  trigger <- match.arg(trigger)
  
  # Try to convert object to data.frame
  # We rely on fortify method being available (e.g. from enrichplot)
  if (inherits(object, "data.frame")) {
    df <- object
  } else {
    # We use ggplot2::fortify which should dispatch if method is available
    df <- tryCatch(
      {
        ggplot2::fortify(object, showCategory = showCategory, ...)
      },
      error = function(e) {
        if (requireNamespace("enrichplot", quietly = TRUE)) {
          return(tryCatch(
            ggplot2::fortify(object, showCategory = showCategory, ...),
            error = function(e2) {
              stop("Cannot convert object to data.frame via ggplot2::fortify(). Please ensure the object has a fortify method.")
            }
          ))
        }
        stop("Cannot convert object to data.frame via ggplot2::fortify(). Please install and load 'enrichplot' if using enrichResult objects.")
      }
    )
  }

  # Handle GeneRatio sorting if x is GeneRatio
  if (x == "GeneRatio" && "GeneRatio" %in% colnames(df) && is.character(df$GeneRatio)) {
    df$GeneRatioVal <- sapply(df$GeneRatio, function(val) eval(parse(text=val)))
    df <- df[order(df$GeneRatioVal), ]
    # Ensure Description is a factor ordered by GeneRatioVal
    if ("Description" %in% colnames(df)) {
        df$Description <- factor(df$Description, levels = unique(df$Description))
    }
  } else if ("Description" %in% colnames(df)) {
      # Default sort by whatever x is if possible, otherwise keep order
      # If x is in df, order by it
      if (x %in% colnames(df)) {
          df <- df[order(df[[x]]), ]
          df$Description <- factor(df$Description, levels = unique(df$Description))
      }
  }

  escape_js_string <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", x, fixed = FALSE)
    x <- gsub("\"", "\\\\\"", x, fixed = TRUE)
    x <- gsub("\n", " ", x, fixed = TRUE)
    x
  }

  if ("Description" %in% colnames(df)) {
    df$data_id <- as.character(df$Description)
  } else {
    df$data_id <- rownames(df)
  }

  if (!"Description" %in% colnames(df)) {
    df$Description <- factor(df$data_id, levels = unique(df$data_id))
  }

  if ("geneID" %in% colnames(df)) {
    df$geneID <- as.character(df$geneID)
    df$geneID[is.na(df$geneID)] <- ""
    df$geneID <- gsub("['\"]+", "", df$geneID)
  } else {
    df$geneID <- ""
    warning("Column 'geneID' not found in data. Interactive highlighting will not work.")
  }

  genes_preview <- vapply(
    strsplit(df$geneID, "/", fixed = TRUE),
    function(x) {
      x <- x[nzchar(x)]
      if (!length(x)) return("")
      paste(utils::head(x, 10), collapse = ", ")
    },
    character(1)
  )

  tooltip_fields <- list(
    Description = if ("Description" %in% colnames(df)) as.character(df$Description) else df$data_id,
    p.adjust = if ("p.adjust" %in% colnames(df)) df[["p.adjust"]] else NA_real_,
    Count = if ("Count" %in% colnames(df)) df[["Count"]] else NA_real_,
    Genes = genes_preview
  )

  df$tooltip <- sprintf(
    "Description: %s<br>p.adjust: %s<br>Count: %s<br>Genes: %s",
    tooltip_fields$Description,
    ifelse(is.na(tooltip_fields$p.adjust), "NA", format(tooltip_fields$p.adjust, digits = 3, scientific = TRUE)),
    ifelse(is.na(tooltip_fields$Count), "NA", as.character(tooltip_fields$Count)),
    tooltip_fields$Genes
  )

  if (trigger == "click") {
    df$onclick <- ifelse(
      nzchar(df$geneID),
      paste0(
        "highlight_genes(\"",
        escape_js_string(df$geneID),
        "\",\"",
        escape_js_string(df$data_id),
        "\")"
      ),
      ""
    )
  } else {
    df$onclick <- ""
  }

  # Handle label formatting
  label_func <- function(x) x
  if (is.numeric(label_format)) {
      label_func <- function(x) {
          vapply(x, function(s) {
              if (nchar(s) > label_format) {
                  paste0(substr(s, 1, label_format), "...")
              } else {
                  s
              }
          }, character(1))
      }
  } else if (is.function(label_format)) {
      label_func <- label_format
  }
  
  if ("Description" %in% colnames(df)) {
      levels(df$Description) <- label_func(levels(df$Description))
  }

  # Plotting
  x_sym <- rlang::sym(x)
  size_sym <- rlang::sym(size)
  color_sym <- rlang::sym(color)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = !!x_sym,
      y = .data$Description,
      size = !!size_sym,
      color = !!color_sym
    )
  ) +
    ggiraph::geom_point_interactive(
      mapping = ggplot2::aes(
        tooltip = .data$tooltip,
        data_id = .data$data_id,
        onclick = .data$onclick,
        `data-genes` = .data$geneID,
        `data-plot` = I("dotplot")
      ),
      extra_interactive_params = c("data-genes", "data-plot")
    ) +
    scale_color_continuous(low = "red", high = "blue") +
    theme_minimal() +
    ylab(NULL) +
    labs(title = title)
    
  return(p)
}
