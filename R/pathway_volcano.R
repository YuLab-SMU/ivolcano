#' Combine Interactive Pathway Dot Plot and Volcano Plot
#'
#' @title pathway_volcano
#' @param p1 interactive dot plot (from idotplot)
#' @param p2 interactive volcano plot (from ivolcano)
#' @param widths relative widths of the two plots, default is c(1, 1)
#' @param trigger interaction trigger, one of "click" or "hover"
#' @param hover_css css for hovered elements
#' @param hover_inv_css css for non-hovered elements
#' @param tooltip_css css for tooltip
#' @param non_selected_color color for non-selected volcano points
#' @param non_selected_opacity opacity for non-selected volcano points
#' @param label_follow whether volcano labels follow pathway selection
#' @param ... additional parameters passed to ggiraph::girafe
#' @return a girafe object
#' @importFrom patchwork wrap_plots
#' @importFrom ggiraph girafe
#' @importFrom htmlwidgets prependContent
#' @importFrom htmltools tags HTML
#' @export
#' @author Guangchuang Yu
pathway_volcano <- function(
  p1,
  p2,
  widths = c(1, 1),
  trigger = c("click", "hover"),
  hover_css = "stroke:red;stroke-width:2px;",
  hover_inv_css = "opacity:0.15;",
  tooltip_css = "background:#fff;color:#000;border:1px solid #999;padding:6px;border-radius:4px;",
  non_selected_color = "#BDBDBD",
  non_selected_opacity = 0.25,
  label_follow = TRUE,
  ...
) {
  trigger <- match.arg(trigger)

  p_combined <- patchwork::wrap_plots(p1, p2, widths = widths)
  g <- ggiraph::girafe(ggobj = p_combined, ...)
  g <- ggiraph::girafe_options(
    g,
    ggiraph::opts_hover(css = hover_css),
    ggiraph::opts_hover_inv(css = hover_inv_css),
    ggiraph::opts_tooltip(css = tooltip_css)
  )

  trigger_js <- paste0("\"", trigger, "\"")
  non_selected_color_js <- paste0("\"", non_selected_color, "\"")
  non_selected_opacity_js <- as.character(non_selected_opacity)
  label_follow_js <- if (isTRUE(label_follow)) "true" else "false"

  label_top_n <- 0L
  label_mode0 <- "separate"
  label_sig_only0 <- TRUE
  label_fontface0 <- "italic"
  if (!is.null(p2$plot_env$label_top_n)) {
    label_top_n <- as.integer(p2$plot_env$label_top_n)
  }
  if (!is.null(p2$plot_env$label_mode)) {
    label_mode0 <- as.character(p2$plot_env$label_mode)
  }
  if (!is.null(p2$plot_env$label_sig_only)) {
    label_sig_only0 <- isTRUE(p2$plot_env$label_sig_only)
  }
  if (!is.null(p2$plot_env$label_fontface)) {
    label_fontface0 <- as.character(p2$plot_env$label_fontface)
  }
  label_top_n_js <- as.character(label_top_n)
  label_mode_js <- paste0("\"", label_mode0, "\"")
  label_sig_only_js <- if (isTRUE(label_sig_only0)) "true" else "false"
  label_fontface_js <- paste0("\"", label_fontface0, "\"")

  js_code <- paste0(
    "(function(){\n",
    "  var TRIGGER = ", trigger_js, ";\n",
    "  var NON_SELECTED_COLOR = ", non_selected_color_js, ";\n",
    "  var NON_SELECTED_OPACITY = ", non_selected_opacity_js, ";\n",
    "  var LABEL_FOLLOW = ", label_follow_js, ";\n",
    "  var LABEL_TOP_N = ", label_top_n_js, ";\n",
    "  var LABEL_MODE = ", label_mode_js, ";\n",
    "  var LABEL_SIG_ONLY = ", label_sig_only_js, ";\n",
    "  var LABEL_FONTFACE = ", label_fontface_js, ";\n",
    "\n",
    "  function toGeneArray(genes){\n",
    "    if(genes === null || genes === undefined) return [];\n",
    "    if(Array.isArray(genes)) return genes.filter(function(x){return x !== null && x !== undefined && String(x).length > 0;}).map(function(x){return String(x).trim();}).filter(function(x){return x.length > 0;});\n",
    "    if(typeof genes === 'string'){\n",
    "      if(genes.length === 0) return [];\n",
    "      return genes.split('/').map(function(x){return x.trim();}).filter(function(x){return x.length > 0;});\n",
    "    }\n",
    "    return [String(genes).trim()].filter(function(x){return x.length > 0;});\n",
    "  }\n",
    "\n",
    "  function volcanoSelection(){\n",
    "    var sel = d3.selectAll('[data-plot=\\\"volcano\\\"]');\n",
    "    if(sel.size && sel.size() > 0) return sel;\n",
    "    return d3.selectAll('circle[data-id]').filter(function(){\n",
    "      return d3.select(this).attr('data-plot') !== 'dotplot';\n",
    "    });\n",
    "  }\n",
    "\n",
    "  function volcanoGeneUniverse(){\n",
    "    var s = new Set();\n",
    "    volcanoSelection().each(function(){\n",
    "      var id = d3.select(this).attr('data-id');\n",
    "      if(id) s.add(id);\n",
    "    });\n",
    "    return s;\n",
    "  }\n",
    "\n",
    "  function labelTexts(){\n",
    "    var universe = volcanoGeneUniverse();\n",
    "    return d3.selectAll('text').filter(function(){\n",
    "      if(this.classList && this.classList.contains('pv-dyn-label')) return false;\n",
    "      var t = (this.textContent || '').trim();\n",
    "      return universe.has(t);\n",
    "    });\n",
    "  }\n",
    "\n",
    "  function clearDynamicLabels(){\n",
    "    d3.selectAll('text.pv-dyn-label').remove();\n",
    "  }\n",
    "\n",
    "  function parseNum(x){\n",
    "    var v = parseFloat(x);\n",
    "    if(isNaN(v)) return null;\n",
    "    return v;\n",
    "  }\n",
    "\n",
    "  function collectLabelCandidates(geneSet){\n",
    "    var out = [];\n",
    "    volcanoSelection().each(function(){\n",
    "      var el = d3.select(this);\n",
    "      var id = el.attr('data-id');\n",
    "      if(!id || !geneSet.has(id)) return;\n",
    "      var p = parseNum(el.attr('data-pval'));\n",
    "      if(p === null) return;\n",
    "      var logfc = parseNum(el.attr('data-logfc'));\n",
    "      var sig = el.attr('data-sig') || '';\n",
    "      if(LABEL_SIG_ONLY && sig === 'Not_Significant') return;\n",
    "      out.push({ id: id, p: p, logfc: (logfc === null ? 0 : logfc) });\n",
    "    });\n",
    "    return out;\n",
    "  }\n",
    "\n",
    "  function pickLabelGenes(cands){\n",
    "    if(!LABEL_TOP_N || LABEL_TOP_N <= 0) return [];\n",
    "    cands.sort(function(a,b){ return a.p - b.p; });\n",
    "    if(LABEL_MODE === 'separate'){\n",
    "      var up = [];\n",
    "      var down = [];\n",
    "      for(var i=0;i<cands.length;i++){\n",
    "        if(cands[i].logfc >= 0) up.push(cands[i]); else down.push(cands[i]);\n",
    "      }\n",
    "      return up.slice(0, LABEL_TOP_N).concat(down.slice(0, LABEL_TOP_N)).map(function(x){ return x.id; });\n",
    "    }\n",
    "    return cands.slice(0, LABEL_TOP_N).map(function(x){ return x.id; });\n",
    "  }\n",
    "\n",
    "  function hideOriginalLabels(){\n",
    "    labelTexts().each(function(){\n",
    "      var t = d3.select(this);\n",
    "      if(t.attr('data-pv-orig-display') === null) t.attr('data-pv-orig-display', t.style('display') || '');\n",
    "      if(t.attr('data-pv-orig-opacity') === null) t.attr('data-pv-orig-opacity', t.style('opacity') || '');\n",
    "      t.style('display', 'none');\n",
    "    });\n",
    "  }\n",
    "\n",
    "  function renderDynamicLabels(ids){\n",
    "    clearDynamicLabels();\n",
    "    if(!ids || ids.length === 0) return;\n",
    "    hideOriginalLabels();\n",
    "    for(var i=0;i<ids.length;i++){\n",
    "      var gene = ids[i];\n",
    "      var node = volcanoSelection().filter(function(){ return d3.select(this).attr('data-id') === gene; }).node();\n",
    "      if(!node) continue;\n",
    "      var el = d3.select(node);\n",
    "      var cx = parseNum(el.attr('cx'));\n",
    "      var cy = parseNum(el.attr('cy'));\n",
    "      if(cx === null || cy === null) continue;\n",
    "      var t = d3.select(node.parentNode).append('text')\n",
    "        .attr('class', 'pv-dyn-label')\n",
    "        .attr('x', cx + 3)\n",
    "        .attr('y', cy - 3)\n",
    "        .style('font-size', '10px')\n",
    "        .style('pointer-events', 'none')\n",
    "        .text(gene);\n",
    "      if(LABEL_FONTFACE){\n",
    "        if(LABEL_FONTFACE.indexOf('bold') >= 0) t.style('font-weight', 'bold');\n",
    "        if(LABEL_FONTFACE.indexOf('italic') >= 0) t.style('font-style', 'italic');\n",
    "      }\n",
    "    }\n",
    "  }\n",
    "\n",
    "  function updateDynamicLabels(geneSet){\n",
    "    if(!LABEL_FOLLOW) return;\n",
    "    if(!LABEL_TOP_N || LABEL_TOP_N <= 0) { clearDynamicLabels(); return; }\n",
    "    var cands = collectLabelCandidates(geneSet);\n",
    "    var ids = pickLabelGenes(cands);\n",
    "    renderDynamicLabels(ids);\n",
    "  }\n",
    "\n",
    "  function cachePoint(el){\n",
    "    if(el.attr('data-pv-cached')) return;\n",
    "    el.attr('data-pv-cached', '1');\n",
    "    el.attr('data-pv-orig-r', el.attr('r'));\n",
    "    var fill = el.attr('fill');\n",
    "    if(fill !== null) el.attr('data-pv-orig-fill', fill);\n",
    "    var stroke = el.attr('stroke');\n",
    "    if(stroke !== null) el.attr('data-pv-orig-stroke', stroke);\n",
    "    var sw = el.attr('stroke-width');\n",
    "    if(sw !== null) el.attr('data-pv-orig-stroke-width', sw);\n",
    "    var op = el.style('opacity');\n",
    "    if(op !== null && op !== '') el.attr('data-pv-orig-opacity', op);\n",
    "  }\n",
    "\n",
    "  function restorePoint(el){\n",
    "    var r = el.attr('data-pv-orig-r');\n",
    "    if(r !== null) el.attr('r', r);\n",
    "    var fill = el.attr('data-pv-orig-fill');\n",
    "    if(fill !== null) el.attr('fill', fill); else el.attr('fill', null);\n",
    "    var stroke = el.attr('data-pv-orig-stroke');\n",
    "    if(stroke !== null) el.attr('stroke', stroke); else el.attr('stroke', null);\n",
    "    var sw = el.attr('data-pv-orig-stroke-width');\n",
    "    if(sw !== null) el.attr('stroke-width', sw); else el.attr('stroke-width', null);\n",
    "    var op = el.attr('data-pv-orig-opacity');\n",
    "    if(op !== null) el.style('opacity', op); else el.style('opacity', null);\n",
    "  }\n",
    "\n",
    "  function resetVolcano(){\n",
    "    clearDynamicLabels();\n",
    "    volcanoSelection().each(function(){\n",
    "      var el = d3.select(this);\n",
    "      cachePoint(el);\n",
    "      restorePoint(el);\n",
    "    });\n",
    "    if(LABEL_FOLLOW){\n",
    "      labelTexts().each(function(){\n",
    "        var t = d3.select(this);\n",
    "        var od = t.attr('data-pv-orig-display');\n",
    "        if(od === null || od === '') t.style('display', null); else t.style('display', od);\n",
    "        var oo = t.attr('data-pv-orig-opacity');\n",
    "        if(oo === null || oo === '') t.style('opacity', null); else t.style('opacity', oo);\n",
    "      });\n",
    "    }\n",
    "    window._pathwayVolcanoActive = null;\n",
    "  }\n",
    "\n",
    "  function applyHighlight(genesArr){\n",
    "    var geneSet = new Set(genesArr);\n",
    "    volcanoSelection().each(function(){\n",
    "      var el = d3.select(this);\n",
    "      cachePoint(el);\n",
    "      var id = el.attr('data-id');\n",
    "      if(geneSet.has(id)) {\n",
    "        restorePoint(el);\n",
    "        var r0 = parseFloat(el.attr('r')) || 0;\n",
    "        el.style('opacity', 1)\n",
    "          .attr('stroke', 'red')\n",
    "          .attr('stroke-width', '2px')\n",
    "          .attr('r', Math.max(r0, 6));\n",
    "        this.parentNode.appendChild(this);\n",
    "      } else {\n",
    "        restorePoint(el);\n",
    "        el.style('opacity', NON_SELECTED_OPACITY)\n",
    "          .attr('fill', NON_SELECTED_COLOR)\n",
    "          .attr('stroke', null)\n",
    "          .attr('stroke-width', null);\n",
    "      }\n",
    "    });\n",
    "\n",
    "    updateDynamicLabels(geneSet);\n",
    "  }\n",
    "\n",
    "  function highlight_genes(genes, pathway_id){\n",
    "    var genesArr = toGeneArray(genes);\n",
    "    if(genesArr.length === 0){\n",
    "      resetVolcano();\n",
    "      return;\n",
    "    }\n",
    "    if(TRIGGER === 'click') {\n",
    "      if(window._pathwayVolcanoActive && window._pathwayVolcanoActive === pathway_id) {\n",
    "        resetVolcano();\n",
    "        return;\n",
    "      }\n",
    "      window._pathwayVolcanoActive = pathway_id || null;\n",
    "      applyHighlight(genesArr);\n",
    "      return;\n",
    "    }\n",
    "    applyHighlight(genesArr);\n",
    "  }\n",
    "\n",
    "  window.highlight_genes = highlight_genes;\n",
    "  window.reset_pathway_volcano = resetVolcano;\n",
    "\n",
    "  function bindHover(el){\n",
    "    var root = d3.select(el);\n",
    "    var pts = root.selectAll('[data-plot=\\\"dotplot\\\"]');\n",
    "    if(!pts.size || pts.size() === 0) {\n",
    "      pts = root.selectAll('[data-genes]');\n",
    "    }\n",
    "    var hoverTimer = null;\n",
    "    pts.on('mouseover', function(){\n",
    "      if(TRIGGER !== 'hover') return;\n",
    "      if(hoverTimer) { clearTimeout(hoverTimer); hoverTimer = null; }\n",
    "      var genes = this.getAttribute('data-genes');\n",
    "      var pid = this.getAttribute('data-id');\n",
    "      highlight_genes(genes, pid);\n",
    "    }).on('mouseout', function(){\n",
    "      if(TRIGGER !== 'hover') return;\n",
    "      if(hoverTimer) clearTimeout(hoverTimer);\n",
    "      hoverTimer = setTimeout(function(){ resetVolcano(); }, 60);\n",
    "    });\n",
    "\n",
    "    root.on('click', function(event){\n",
    "      if(TRIGGER !== 'click') return;\n",
    "      var t = event && event.target ? event.target : null;\n",
    "      if(!t) return;\n",
    "      var tag = (t.tagName || '').toLowerCase();\n",
    "      if(tag === 'circle' || tag === 'path' || tag === 'text') return;\n",
    "      resetVolcano();\n",
    "    });\n",
    "  }\n",
    "\n",
    "  if(window.HTMLWidgets && HTMLWidgets.addPostRenderHandler){\n",
    "    HTMLWidgets.addPostRenderHandler(function(el){ bindHover(el); });\n",
    "  } else {\n",
    "    document.addEventListener('DOMContentLoaded', function(){\n",
    "      var nodes = document.querySelectorAll('.girafe');\n",
    "      for(var i=0;i<nodes.length;i++){ bindHover(nodes[i]); }\n",
    "    });\n",
    "  }\n",
    "})();"
  )

  g <- htmlwidgets::prependContent(
    g,
    htmltools::tags$script(htmltools::HTML(js_code))
  )

  g
}
