#' Launch a pathway-selection Shiny explorer
#'
#' The application presents pathway enrichment results in a searchable DT
#' table. Selecting one or more rows reactively highlights the union of their
#' matched genes in a ggplot2/ggiraph volcano plot and lists the selected genes
#' below the plot.
#'
#' @param enrichment enrichment result object or data frame. Objects supported
#'   by `ggplot2::fortify()` (for example, `clusterProfiler::enrichResult`)
#'   are converted automatically.
#' @param gene_data data frame containing one row per gene and its effect size.
#' @param pathway_col column containing pathway descriptions.
#' @param pathway_gene_col column containing pathway member IDs. IDs may be
#'   separated by `/`, `;`, `,`, or `|`.
#' @param pathway_p_col pathway adjusted P-value column. If absent, the table
#'   displays `NA` for FDR.
#' @param pathway_ratio_col optional GeneRatio column.
#' @param gene_col gene identifier column in `gene_data`.
#' @param gene_effect_col gene effect column, usually logFC.
#' @param gene_p_col gene P-value column, usually adjusted P-value.
#' @param show_pathways optional number of pathways to show. `NULL` keeps all
#'   pathways.
#' @param gene_id_case whether IDs should be matched as-is, upper-cased, or
#'   lower-cased.
#' @param pval_cutoff gene P-value cutoff used for volcano classification.
#' @param logFC_cutoff absolute effect cutoff used for volcano classification.
#' @param launch whether to immediately launch the application with
#'   `shiny::runApp()`. The default returns a `shiny.appobj` for embedding or
#'   testing.
#' @param app_title title shown above the application.
#' @return A `shiny.appobj` invisibly when `launch = FALSE`; when `launch = TRUE`
#'   the function runs the application and returns after it is closed.
#' @importFrom ggiraph girafe girafeOutput geom_point_interactive opts_hover renderGirafe
#' @importFrom ggplot2 fortify ggplot aes geom_hline geom_vline labs scale_colour_manual
#'   scale_alpha_identity scale_size_identity theme_minimal theme element_text
#'   element_blank margin
#' @export
#' @author Guangchuang Yu
pathway_app <- function(
  enrichment,
  gene_data,
  pathway_col = "Description",
  pathway_gene_col = "geneID",
  pathway_p_col = "p.adjust",
  pathway_ratio_col = "GeneRatio",
  gene_col = "gene",
  gene_effect_col = "logFC",
  gene_p_col = "adj.P.Val",
  show_pathways = NULL,
  gene_id_case = c("asis", "upper", "lower"),
  pval_cutoff = 0.05,
  logFC_cutoff = 1,
  launch = FALSE,
  app_title = "Pathway-gene explorer"
) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The pathway_app() function requires the 'shiny' package.", call. = FALSE)
  }
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("The pathway_app() function requires the 'DT' package.", call. = FALSE)
  }
  if (!requireNamespace("ggiraph", quietly = TRUE)) {
    stop("The pathway_app() function requires the 'ggiraph' package.", call. = FALSE)
  }
  gene_id_case <- match.arg(gene_id_case)
  prepared <- .pathway_app_prepare(
    enrichment = enrichment,
    gene_data = gene_data,
    pathway_col = pathway_col,
    pathway_gene_col = pathway_gene_col,
    pathway_p_col = pathway_p_col,
    pathway_ratio_col = pathway_ratio_col,
    gene_col = gene_col,
    gene_effect_col = gene_effect_col,
    gene_p_col = gene_p_col,
    show_pathways = show_pathways,
    gene_id_case = gene_id_case,
    pval_cutoff = pval_cutoff,
    logFC_cutoff = logFC_cutoff
  )

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML(
        paste0(
          ".pa-shell{max-width:1400px;margin:0 auto;padding:12px 20px;} ",
          ".pa-status{padding:8px 12px;background:#F8FAFC;border:1px solid #E2E8F0;border-radius:6px;margin:10px 0;color:#334155;} ",
          ".pa-help{color:#64748B;font-size:13px;margin:4px 0 10px;} ",
          ".pa-panel{margin-bottom:18px;} ",
          ".pa-panel h3{margin:8px 0 6px;color:#0F172A;}"
        )
      ))
    ),
    shiny::div(
      class = "pa-shell",
      shiny::h2(app_title),
      shiny::div(class = "pa-help", "Select one or more pathways in the table to highlight their matched genes in the volcano plot."),
      shiny::div(class = "pa-panel", DT::DTOutput("pathway_table")),
      shiny::div(class = "pa-status", shiny::textOutput("selection_status")),
      shiny::div(class = "pa-panel", ggiraph::girafeOutput("volcano", width = "100%", height = "620px")),
      shiny::div(class = "pa-panel", shiny::h3("Selected genes"), DT::DTOutput("selected_genes")),
      shiny::div(class = "pa-panel", shiny::h3("Matching diagnostics"), shiny::verbatimTextOutput("diagnostics"))
    )
  )

  server <- function(input, output, session) {
    selected_rows <- shiny::reactive({
      rows <- input$pathway_table_rows_selected
      if (is.null(rows)) integer() else as.integer(rows)
    })

    selected_gene_keys <- shiny::reactive({
      rows <- selected_rows()
      if (!length(rows)) return(character())
      unique(unlist(prepared$pathway_table$gene_keys[rows], use.names = FALSE))
    })

    selected_pathway_names <- shiny::reactive({
      rows <- selected_rows()
      if (!length(rows)) character() else prepared$pathway_table$pathway[rows]
    })

    output$pathway_table <- DT::renderDT({
      DT::datatable(
        prepared$pathway_display,
        rownames = FALSE,
        filter = "top",
        selection = list(mode = "multiple", target = "row"),
        options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
        class = "compact stripe hover"
      )
    }, server = FALSE)

    output$selection_status <- shiny::renderText({
      keys <- selected_gene_keys()
      paths <- selected_pathway_names()
      sprintf(
        "Selected pathways: %d | Matched genes: %d | %s",
        length(paths),
        length(keys),
        if (length(paths)) paste(paths, collapse = "; ") else "No pathway selected"
      )
    })

    output$volcano <- ggiraph::renderGirafe({
      dat <- prepared$gene_table
      keys <- selected_gene_keys()
      dat$selected <- dat$gene_key %in% keys
      dat$point_alpha <- if (length(keys)) ifelse(dat$selected, 1, 0.16) else 0.7
      dat$point_size <- ifelse(dat$selected, 2.8, 1.7)
      dat$tooltip <- sprintf(
        "Gene: %s<br>logFC: %s<br>P-value: %s%s",
        dat$gene_label,
        format(dat$gene_effect, digits = 3),
        format(dat$gene_p, digits = 3, scientific = TRUE),
        ifelse(dat$selected, "<br>Selected pathway gene", "")
      )
      p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$gene_effect, y = .data$neg_log_p)) +
        ggiraph::geom_point_interactive(
          ggplot2::aes(
            colour = .data$sig,
            alpha = .data$point_alpha,
            size = .data$point_size,
            tooltip = .data$tooltip,
            data_id = .data$gene_key
          )
        ) +
        ggplot2::scale_colour_manual(
          values = c("Down" = "#0072B2", "NS" = "#CBD5E1", "Up" = "#D55E00"),
          drop = FALSE,
          name = NULL
        ) +
        ggplot2::scale_alpha_identity(guide = "none") +
        ggplot2::scale_size_identity(guide = "none") +
        ggplot2::geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", colour = "#64748B", linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", colour = "#64748B", linewidth = 0.4) +
        ggplot2::labs(x = gene_effect_col, y = expression(-log[10](P))) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          legend.position = "bottom",
          plot.margin = ggplot2::margin(8, 12, 8, 12)
        )
      selected_dat <- dat[dat$selected, , drop = FALSE]
      if (nrow(selected_dat)) {
        p <- p + ggiraph::geom_point_interactive(
          data = selected_dat,
          ggplot2::aes(
            x = .data$gene_effect,
            y = .data$neg_log_p,
            tooltip = .data$tooltip,
            data_id = .data$gene_key
          ),
          inherit.aes = FALSE,
          shape = 21,
          fill = NA,
          colour = "#111827",
          stroke = 1.1,
          size = 4.2
        )
      }
      ggiraph::girafe(
        ggobj = p,
        options = list(
          ggiraph::opts_hover(css = "stroke:#111827;stroke-width:2px;opacity:1;"),
          ggiraph::opts_toolbar(saveaspng = TRUE)
        )
      )
    })

    output$selected_genes <- DT::renderDT({
      keys <- selected_gene_keys()
      if (!length(keys)) {
        return(DT::datatable(
          data.frame(Message = "Select one or more pathways above."),
          rownames = FALSE,
          options = list(dom = "t")
        ))
      }
      rows <- prepared$gene_table$gene_key %in% keys
      selected <- prepared$gene_table[rows, c("gene_label", "gene_effect", "gene_p", "sig"), drop = FALSE]
      selected$pathways <- vapply(prepared$gene_table$gene_key[rows], function(key) {
        hits <- vapply(prepared$pathway_table$gene_keys, function(x) key %in% x, logical(1))
        paste(prepared$pathway_table$pathway[hits & seq_along(hits) %in% selected_rows()], collapse = "; ")
      }, character(1))
      names(selected) <- c("Gene", gene_effect_col, gene_p_col, "Direction", "Selected pathways")
      DT::datatable(
        selected,
        rownames = FALSE,
        filter = "top",
        options = list(pageLength = 10, scrollX = TRUE),
        class = "compact stripe hover"
      )
    }, server = FALSE)

    output$diagnostics <- shiny::renderText({
      pt <- prepared$pathway_table
      sprintf(
        paste0(
          "Pathways: %d\n",
          "Pathways with >=1 matched gene: %d\n",
          "Genes in gene_data: %d\n",
          "Genes referenced by pathways: %d\n",
          "Matched gene IDs: %d"
        ),
        nrow(pt),
        sum(pt$matched_count > 0),
        nrow(prepared$gene_table),
        length(unique(unlist(pt$gene_keys, use.names = FALSE))),
        length(prepared$matched_gene_keys)
      )
    })
  }

  app <- shiny::shinyApp(ui = ui, server = server)
  if (isTRUE(launch)) shiny::runApp(app)
  app
}

.pathway_app_prepare <- function(
  enrichment,
  gene_data,
  pathway_col,
  pathway_gene_col,
  pathway_p_col,
  pathway_ratio_col,
  gene_col,
  gene_effect_col,
  gene_p_col,
  show_pathways,
  gene_id_case,
  pval_cutoff,
  logFC_cutoff
) {
  if (!is.data.frame(gene_data)) stop("gene_data must be a data.frame.", call. = FALSE)
  pathway_df <- if (is.data.frame(enrichment)) {
    enrichment
  } else if (is.null(show_pathways)) {
    tryCatch(ggplot2::fortify(enrichment), error = function(e) {
      stop("Cannot convert enrichment to a data.frame via ggplot2::fortify().", call. = FALSE)
    })
  } else {
    tryCatch(ggplot2::fortify(enrichment, showCategory = show_pathways), error = function(e) {
      tryCatch(ggplot2::fortify(enrichment), error = function(e2) {
        stop("Cannot convert enrichment to a data.frame via ggplot2::fortify().", call. = FALSE)
      })
    })
  }
  missing_pathway <- setdiff(c(pathway_col, pathway_gene_col), names(pathway_df))
  if (length(missing_pathway)) stop("Missing enrichment columns: ", paste(missing_pathway, collapse = ", "), call. = FALSE)
  missing_gene <- setdiff(c(gene_col, gene_effect_col, gene_p_col), names(gene_data))
  if (length(missing_gene)) stop("Missing gene_data columns: ", paste(missing_gene, collapse = ", "), call. = FALSE)

  normalise <- function(x) {
    x <- trimws(as.character(x))
    if (gene_id_case == "upper") x <- toupper(x)
    if (gene_id_case == "lower") x <- tolower(x)
    x
  }
  split_ids <- function(x) {
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    if (!length(x) || all(is.na(x))) return(character())
    ids <- unlist(strsplit(gsub("['\"]+", "", paste(as.character(x), collapse = "/")), "[/;,|]"), use.names = FALSE)
    unique(normalise(ids[nzchar(trimws(ids)) & !is.na(ids)]))
  }
  ratio_value <- function(x) {
    if (!length(x) || is.na(x)) return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x)[1])
    bits <- strsplit(as.character(x)[1], "/", fixed = TRUE)[[1]]
    if (length(bits) != 2) return(suppressWarnings(as.numeric(x)[1]))
    as.numeric(bits[1]) / as.numeric(bits[2])
  }

  gene_table <- gene_data
  gene_table$gene_label <- as.character(gene_table[[gene_col]])
  gene_table$gene_key <- normalise(gene_table[[gene_col]])
  gene_table$gene_effect <- suppressWarnings(as.numeric(gene_table[[gene_effect_col]]))
  gene_table$gene_p <- suppressWarnings(as.numeric(gene_table[[gene_p_col]]))
  gene_table <- gene_table[!is.na(gene_table$gene_key) & nzchar(gene_table$gene_key), , drop = FALSE]
  gene_table <- gene_table[!duplicated(gene_table$gene_key), , drop = FALSE]
  gene_table$neg_log_p <- -log10(pmax(gene_table$gene_p, .Machine$double.xmin))
  gene_table$neg_log_p[!is.finite(gene_table$neg_log_p)] <- 0
  gene_table$sig <- "NS"
  gene_table$sig[gene_table$gene_p < pval_cutoff & gene_table$gene_effect > logFC_cutoff] <- "Up"
  gene_table$sig[gene_table$gene_p < pval_cutoff & gene_table$gene_effect < -logFC_cutoff] <- "Down"

  pathway_table <- pathway_df[!is.na(pathway_df[[pathway_col]]), , drop = FALSE]
  pathway_table$pathway <- as.character(pathway_table[[pathway_col]])
  pathway_table$pathway_id <- make.unique(pathway_table$pathway)
  pathway_table$gene_keys <- lapply(pathway_table[[pathway_gene_col]], split_ids)
  pathway_table$total_count <- lengths(pathway_table$gene_keys)
  pathway_table$matched_count <- vapply(pathway_table$gene_keys, function(x) sum(x %in% gene_table$gene_key), integer(1))
  pathway_table$fdr <- if (pathway_p_col %in% names(pathway_table)) suppressWarnings(as.numeric(pathway_table[[pathway_p_col]])) else NA_real_
  pathway_table$gene_ratio <- if (pathway_ratio_col %in% names(pathway_table)) vapply(pathway_table[[pathway_ratio_col]], ratio_value, numeric(1)) else NA_real_
  pathway_table <- pathway_table[order(pathway_table$fdr, na.last = TRUE), , drop = FALSE]
  if (!is.null(show_pathways)) pathway_table <- utils::head(pathway_table, show_pathways)
  pathway_table$gene_keys <- lapply(pathway_table$gene_keys, unique)

  pathway_display <- data.frame(
    Pathway = pathway_table$pathway,
    FDR = format(pathway_table$fdr, digits = 3, scientific = TRUE),
    GeneRatio = ifelse(is.na(pathway_table$gene_ratio), "NA", format(pathway_table$gene_ratio, digits = 3)),
    `Matched genes` = pathway_table$matched_count,
    `Total genes` = pathway_table$total_count,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    pathway_table = pathway_table,
    pathway_display = pathway_display,
    gene_table = gene_table,
    matched_gene_keys = intersect(unique(unlist(pathway_table$gene_keys, use.names = FALSE)), gene_table$gene_key)
  )
}
