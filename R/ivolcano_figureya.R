#' FigureYa volcano plot theme
#' 
#' @title ivolcano_figureya
#' @param data A data frame that contains minimal information with gene id, logFC and P values
#' @param logFC_col column name in 'data' that stored the logFC values
#' @param pval_col column name in 'data' that stored the P values (P.Value or adj.P.Val)
#' @param gene_col column name in 'data' that stored the gene IDs
#' @param selected_genes optional data frame with genes to highlight (first column should be gene IDs, second column can be pathway)
#' @param plot_mode one of 'classic' or 'advanced' (default: 'advanced')
#' @param logFCcut log2 fold change cutoff for basic threshold (default: 1.5)
#' @param logFCcut2 log2 fold change cutoff for second threshold in advanced mode (default: 2.5)
#' @param logFCcut3 log2 fold change cutoff for third threshold in advanced mode (default: 5)
#' @param pvalCut P value cutoff for basic threshold (default: 0.05)
#' @param pvalCut2 P value cutoff for second threshold in advanced mode (default: 0.0001)
#' @param pvalCut3 P value cutoff for third threshold in advanced mode (default: 0.00001)
#' @param xlim custom x-axis limits (default: NULL, will use symmetric range)
#' @param ylim custom y-axis limits (default: NULL, will auto-calculate)
#' @param label_extreme logical, whether to label genes with extreme logFC values (default: FALSE)
#' @param extreme_logFC_threshold threshold for extreme logFC labeling (default: 9)
#' @param show_pathway logical, whether to show pathway legend for selected genes (default: FALSE)
#' @param pathway_colors custom colors for pathways (default: NULL, will use predefined colors)
#' @param title plot title (default: "")
#' @param interactive whether plot the graph in interactive mode (default: FALSE)
#' @param onclick_fun effects when click on the dot (gene), default is NULL
#' @param point_colors custom colors for points. For classic mode: named vector/list with 'up', 'down', 'ns'. For advanced mode: named vector/list with 'up1', 'up2', 'down1', 'down2', 'ns' (default: NULL, uses default colors)
#' @param point_alpha transparency of points (default: 0.6 for main points, 0.5 for advanced mode)
#' @param point_size_base base size for non-significant points (default: NULL, uses mode-specific defaults)
#' @param point_size_sig size for significant points in classic mode (default: NULL, uses 4)
#' @param point_size_small size for first threshold in advanced mode (default: NULL, uses 2)
#' @param point_size_medium size for second threshold in advanced mode (default: NULL, uses 4)
#' @param point_size_large size for third threshold in advanced mode (default: NULL, uses 6)
#' @param label_color color for gene labels (default: NULL, uses pathway colors if available, otherwise black)
#' @param label_color_extreme color for extreme logFC labels (default: "darkred")
#' @param label_size size for gene labels (default: 5)
#' @return volcano plot with FigureYa theme
#' @importFrom dplyr mutate filter case_when
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline labs scale_x_continuous coord_cartesian theme_bw theme scale_color_manual guides guide_legend
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggiraph geom_point_interactive girafe opts_hover
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom grid unit
#' @importFrom grDevices adjustcolor
#' @importFrom rlang sym
#' @importFrom rlang .data
#' @export
#' @examples
#' # Load example data
#' f1 <- system.file('extdata/easy_input_limma.rds', package='ivolcano')
#' df <- readRDS(f1)
#' 
#' # Basic plot
#' ivolcano_figureya(df, 
#'                   logFC_col = "logFC",
#'                   pval_col = "P.Value",
#'                   gene_col = "X")
#' 
#' # With selected genes
#' f2 <- system.file('extdata/easy_input_selected.rds', package='ivolcano')
#' selected <- readRDS(f2)
#' 
#' ivolcano_figureya(df,
#'                   logFC_col = "logFC",
#'                   pval_col = "P.Value",
#'                   gene_col = "X",
#'                   selected_genes = selected,
#'                   show_pathway = TRUE)
#' 
#' # Interactive mode with onclick function
#' ivolcano_figureya(df,
#'                   logFC_col = "logFC",
#'                   pval_col = "P.Value",
#'                   gene_col = "X",
#'                   interactive = TRUE,
#'                   onclick_fun = onclick_genecards)
#' @author Based on FigureYa59volcanoV2.Rmd by Haitao Wang
ivolcano_figureya <- function(
  data,
  logFC_col = "logFC",
  pval_col = "P.Value",
  gene_col = "X",
  selected_genes = NULL,
  plot_mode = "advanced",
  logFCcut = 1.5,
  logFCcut2 = 2.5,
  logFCcut3 = 5,
  pvalCut = 0.05,
  pvalCut2 = 0.0001,
  pvalCut3 = 0.00001,
  xlim = NULL,
  ylim = NULL,
  label_extreme = FALSE,
  extreme_logFC_threshold = 9,
  show_pathway = FALSE,
  pathway_colors = NULL,
  title = "",
  interactive = FALSE,
  onclick_fun = NULL,
  point_colors = NULL,
  point_alpha = NULL,
  point_size_base = NULL,
  point_size_sig = NULL,
  point_size_small = NULL,
  point_size_medium = NULL,
  point_size_large = NULL,
  label_color = NULL,
  label_color_extreme = "darkred",
  label_size = 5
) {
  # Validate inputs
  stopifnot(all(c(logFC_col, pval_col, gene_col) %in% colnames(data)))
  plot_mode <- match.arg(plot_mode, c("classic", "advanced"))
  
  # Prepare data
  x <- data[!is.na(data[[logFC_col]]) & !is.na(data[[pval_col]]), ]
  x$label <- x[[gene_col]]
  
  # Calculate y-axis range if not provided
  if (is.null(ylim)) {
    ymin <- 0
    ymax <- max(-log10(x[[pval_col]])) * 1.1
  } else {
    ymin <- ylim[1]
    ymax <- ylim[2]
  }
  
  # Calculate x-axis range if not provided
  if (is.null(xlim)) {
    x_range <- range(x[[logFC_col]], na.rm = TRUE)
    # 确保范围包含所有数据，并添加适当的边距
    xmin_val <- x_range[1]
    xmax_val <- x_range[2]
    # 计算边距（数据范围的 5%，但至少为 1）
    margin <- max(1, (xmax_val - xmin_val) * 0.05)
    # 如果数据不对称，确保范围是对称的（可选，或者保持原始范围）
    # 这里我们保持原始范围，但添加边距
    xlim <- c(xmin_val - margin, xmax_val + margin)
  }
  
  # Set default colors and sizes based on mode
  if (plot_mode == "classic") {
    # Default colors
    default_colors <- list(up = "red", down = "blue", ns = "grey")
    if (!is.null(point_colors)) {
      if (is.list(point_colors) || is.vector(point_colors)) {
        default_colors <- modifyList(default_colors, as.list(point_colors))
      }
    }
    
    # Simple color setting
    x$color_transparent <- ifelse(
      (x[[pval_col]] < pvalCut & x[[logFC_col]] > logFCcut), default_colors$up,
      ifelse(
        (x[[pval_col]] < pvalCut & x[[logFC_col]] < -logFCcut), default_colors$down,
        default_colors$ns
      )
    )
    
    # Simple size setting
    size_base <- ifelse(is.null(point_size_base), 2, point_size_base)
    size_sig <- ifelse(is.null(point_size_sig), 4, point_size_sig)
    size <- ifelse(
      (x[[pval_col]] < pvalCut & abs(x[[logFC_col]]) > logFCcut), size_sig, size_base
    )
  } else if (plot_mode == "advanced") {
    # Default colors
    default_colors <- list(
      up1 = "#FB9A99", up2 = "#ED4F4F",
      down1 = "#B2DF8A", down2 = "#329E3F",
      ns = "grey"
    )
    if (!is.null(point_colors)) {
      if (is.list(point_colors) || is.vector(point_colors)) {
        default_colors <- modifyList(default_colors, as.list(point_colors))
      }
    }
    
    # Complex color setting
    n1 <- nrow(x)
    cols <- rep(default_colors$ns, n1)
    
    # Different threshold colors
    cols[x[[pval_col]] < pvalCut & x[[logFC_col]] > logFCcut] <- default_colors$up1
    cols[x[[pval_col]] < pvalCut2 & x[[logFC_col]] > logFCcut2] <- default_colors$up2
    cols[x[[pval_col]] < pvalCut & x[[logFC_col]] < -logFCcut] <- default_colors$down1
    cols[x[[pval_col]] < pvalCut2 & x[[logFC_col]] < -logFCcut2] <- default_colors$down2
    
    # Set alpha
    alpha_val <- ifelse(is.null(point_alpha), 0.5, point_alpha)
    color_transparent <- adjustcolor(cols, alpha.f = alpha_val)
    x$color_transparent <- color_transparent
    
    # Complex size setting
    size_base <- ifelse(is.null(point_size_base), 1, point_size_base)
    size_small <- ifelse(is.null(point_size_small), 2, point_size_small)
    size_medium <- ifelse(is.null(point_size_medium), 4, point_size_medium)
    size_large <- ifelse(is.null(point_size_large), 6, point_size_large)
    
    size <- rep(size_base, n1)
    size[x[[pval_col]] < pvalCut & x[[logFC_col]] > logFCcut] <- size_small
    size[x[[pval_col]] < pvalCut2 & x[[logFC_col]] > logFCcut2] <- size_medium
    size[x[[pval_col]] < pvalCut3 & x[[logFC_col]] > logFCcut3] <- size_large
    size[x[[pval_col]] < pvalCut & x[[logFC_col]] < -logFCcut] <- size_small
    size[x[[pval_col]] < pvalCut2 & x[[logFC_col]] < -logFCcut2] <- size_medium
    size[x[[pval_col]] < pvalCut3 & x[[logFC_col]] < -logFCcut3] <- size_large
  }
  
  # Set alpha for classic mode if provided
  if (plot_mode == "classic") {
    alpha_val <- ifelse(is.null(point_alpha), 0.6, point_alpha)
    # Apply alpha to colors
    x$color_transparent <- adjustcolor(x$color_transparent, alpha.f = alpha_val)
  }
  
  # Prepare selected genes data if provided
  selectgenes <- NULL
  if (!is.null(selected_genes)) {
    # Assume first column is gene ID, second column (if exists) is pathway
    gene_col_sel <- colnames(selected_genes)[1]
    x$gsym <- x[[gene_col]]
    selectgenes <- merge(selected_genes, x, by.x = gene_col_sel, by.y = "gsym", all.x = FALSE)
    if (nrow(selectgenes) == 0) {
      warning("No matching genes found in selected_genes. Ignoring selected_genes.")
      selectgenes <- NULL
    }
  }
  
  # Default pathway colors
  if (is.null(pathway_colors)) {
    mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", 
               "#088247", "#58CDD9", "#7A142C", "#5D90BA", "#431A3D", 
               "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D", "#7CC767")
  } else {
    mycol <- pathway_colors
  }
  
  # Prepare onclick and tooltip for interactive mode
  if (interactive) {
    # onclick
    if (!is.null(onclick_fun)) {
      if (is.function(onclick_fun)) {
        x$onclick <- vapply(
          x[[gene_col]],
          function(g) {
            result <- onclick_fun(g)
            if (is.na(result) || !nzchar(result)) "" else as.character(result)
          },
          character(1)
        )
      } else {
        stop("onclick_fun must be a function or NULL")
      }
    } else {
      x$onclick <- ""
    }
    
    # Build tooltip - use sprintf and clean quotes
    x$tooltip <- sprintf(
      "%s: %s\nlogFC: %.3f\nP-value: %.3e",
      gene_col,
      x[[gene_col]],
      x[[logFC_col]],
      x[[pval_col]]
    )
    # Remove single and double quotes from tooltip
    x$tooltip <- gsub("['\"]", "", x$tooltip)
  }
  
  # Build plot
  if (interactive) {
    aes_args <- ggplot2::aes(
      x = !!rlang::sym(logFC_col),
      y = -log10(!!rlang::sym(pval_col)),
      label = .data$label,
      tooltip = .data$tooltip,
      data_id = !!rlang::sym(gene_col),
      onclick = .data$onclick
    )
  } else {
    aes_args <- ggplot2::aes(
      x = !!rlang::sym(logFC_col),
      y = -log10(!!rlang::sym(pval_col)),
      label = .data$label
    )
  }
  
  # Determine alpha for points (already applied in color_transparent for classic mode)
  point_alpha_val <- if (plot_mode == "classic") {
    1.0  # Alpha already in color
  } else {
    ifelse(is.null(point_alpha), 0.6, point_alpha)
  }
  
  p1 <- ggplot2::ggplot(data = x, aes_args) +
    {
      if (interactive) {
        ggiraph::geom_point_interactive(alpha = point_alpha_val, size = size, colour = x$color_transparent)
      } else {
        ggplot2::geom_point(alpha = point_alpha_val, size = size, colour = x$color_transparent)
      }
    } +
    ggplot2::labs(x = bquote(~Log[2]~"(fold change)"), 
         y = bquote(~-Log[10]~italic("P-value")), 
         title = title) +
    # Use coord_cartesian for both axes to avoid data clipping
    ggplot2::coord_cartesian(ylim = c(ymin, ymax), xlim = xlim) +
    {
      # Dynamically calculate x-axis breaks based on data range
      x_range_full <- range(x[[logFC_col]], na.rm = TRUE)
      x_abs_max <- max(abs(x_range_full))
      
      if (x_abs_max > 10) {
        # Large range: use larger step size
        step <- ifelse(x_abs_max > 20, 5, 2)
        x_breaks <- seq(
          floor(min(xlim) / step) * step,
          ceiling(max(xlim) / step) * step,
          by = step
        )
        # Ensure key threshold points are included
        key_points <- c(-logFCcut, 0, logFCcut)
        x_breaks <- sort(unique(c(x_breaks, key_points)))
        # Filter breaks within the xlim range
        x_breaks <- x_breaks[x_breaks >= min(xlim) & x_breaks <= max(xlim)]
      } else {
        # Small range: use finer breaks
        x_breaks <- c(-10, -5, -logFCcut, 0, logFCcut, 5, 10)
        # Filter breaks within the xlim range
        x_breaks <- x_breaks[x_breaks >= min(xlim) & x_breaks <= max(xlim)]
      }
      
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        labels = x_breaks
      )
    } +
    # Draw threshold boundary lines
    ggplot2::geom_vline(xintercept = c(-logFCcut, logFCcut), 
               color = "grey40", linetype = "longdash", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(pvalCut), 
               color = "grey40", linetype = "longdash", linewidth = 0.5) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  
  # Add additional threshold lines for advanced mode
  if (plot_mode == "advanced") {
    p1 <- p1 +
      ggplot2::geom_vline(xintercept = c(-logFCcut2, logFCcut2), 
                 color = "grey40", linetype = "longdash", linewidth = 0.5) +
      ggplot2::geom_hline(yintercept = -log10(pvalCut2), 
                 color = "grey40", linetype = "longdash", linewidth = 0.5)
  }
  
  # Label extreme genes if requested
  if (label_extreme) {
    p1 <- p1 +
      ggrepel::geom_text_repel(
        ggplot2::aes(x = !!rlang::sym(logFC_col), y = -log10(!!rlang::sym(pval_col)),
            label = ifelse(!!rlang::sym(logFC_col) > extreme_logFC_threshold, .data$label, "")),
        colour = label_color_extreme, size = label_size, 
        box.padding = grid::unit(0.35, "lines"),
        point.padding = grid::unit(0.3, "lines")
      )
  }
  
  # Highlight selected genes if provided
  if (!is.null(selectgenes) && nrow(selectgenes) > 0) {
    # Get pathway column name if exists
    pathway_col <- if (ncol(selected_genes) > 1) colnames(selected_genes)[2] else NULL
    
    # Add pathway colors if pathway column exists
    if (!is.null(pathway_col) && pathway_col %in% colnames(selectgenes)) {
      # Map pathways to colors
      unique_pathways <- unique(selectgenes[[pathway_col]])
      np <- length(unique_pathways)
      pathway_color_map <- setNames(mycol[1:np], unique_pathways)
      
      # For selected genes, use non-interactive layers to avoid aes conflicts
      p2 <- p1 +
        # Draw black circle around selected genes (non-interactive)
        ggplot2::geom_point(
          data = selectgenes,
          ggplot2::aes(x = !!rlang::sym(logFC_col), y = -log10(!!rlang::sym(pval_col))),
          alpha = 1, size = 4.6, shape = 1,
          stroke = 1, color = "black", inherit.aes = FALSE
        ) +
        # Display gene names with pathway colors
        ggrepel::geom_text_repel(
          data = selectgenes,
          ggplot2::aes(
            x = !!rlang::sym(logFC_col),
            y = -log10(!!rlang::sym(pval_col)),
            color = !!rlang::sym(pathway_col),
            label = .data$label
          ),
          show.legend = show_pathway,
          size = label_size,
          box.padding = grid::unit(0.35, "lines"),
          point.padding = grid::unit(0.3, "lines"),
          inherit.aes = FALSE
        ) +
        ggplot2::scale_color_manual(values = pathway_color_map) +
        ggplot2::guides(color = ggplot2::guide_legend(title = NULL))
      
      # Show pathway legend as table if requested
      if (show_pathway) {
        labelsInfo <- data.frame(
          pathway = unique_pathways,
          col = mycol[1:np]
        )
        # Create a simple text table for pathway legend
        pathway_table <- gridExtra::tableGrob(
          matrix(labelsInfo$pathway, ncol = 1),
          rows = rep("", np),
          cols = "",
          theme = gridExtra::ttheme_minimal()
        )
        # Apply colors to text (this is a simplified approach)
        p2 <- p2 +
          ggplot2::annotation_custom(
            pathway_table,
            ymin = ymax - 2, ymax = ymax,
            xmin = xlim[1] - 1.5, xmax = xlim[1]
          )
      }
      
      p_final <- p2
    } else {
      # Determine label color for selected genes
      label_color_sel <- ifelse(is.null(label_color), "black", label_color)
      
      # For selected genes, use non-interactive layers to avoid aes conflicts
      p2 <- p1 +
        # Draw black circle around selected genes (non-interactive)
        ggplot2::geom_point(
          data = selectgenes,
          ggplot2::aes(x = !!rlang::sym(logFC_col), y = -log10(!!rlang::sym(pval_col))),
          alpha = 1, size = 4.6, shape = 1,
          stroke = 1, color = "black", inherit.aes = FALSE
        ) +
        # Display gene names
        ggrepel::geom_text_repel(
          data = selectgenes,
          ggplot2::aes(
            x = !!rlang::sym(logFC_col),
            y = -log10(!!rlang::sym(pval_col)),
            label = .data$label
          ),
          show.legend = FALSE,
          colour = label_color_sel,
          size = label_size,
          box.padding = grid::unit(0.35, "lines"),
          point.padding = grid::unit(0.3, "lines"),
          inherit.aes = FALSE
        )
      
      p_final <- p2
    }
  } else {
    p_final <- p1
  }
  
  # Return interactive or static plot
  if (interactive) {
    return(ggiraph::girafe(
      ggobj = p_final,
      options = list(ggiraph::opts_hover(css = "fill:black;r:6"))
    ))
  } else {
    return(p_final)
  }
}

