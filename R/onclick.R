#' generate onclick function from `fanyi::gene_summary()` output
#' 
#' @title onclick_fanyi
#' @param gene_summary output from `fanyi::gene_summary()`
#' @param cols selected columns from 'gene_summary'
#' @return onclick function
#' @export
#' @author Guangchuang Yu
onclick_fanyi <- function(gene_summary, cols) {
  function(gene) {
    if (!gene %in% rownames(gene_summary)) {
      return("")
    }

    # extract columns, support multiple columns
    info <- gene_summary[gene, cols, drop = FALSE]
    txt <- paste(apply(info, 1, paste, collapse = "\n"), collapse = "\n")

    # safe to use ggiraph
    txt <- gsub("['\"]", "", txt) # rm single and double quote
    txt <- gsub("\n", " | ", txt) # replace new line

    # JS alert
    paste0("alert(\"", gene, ": ", txt, "\")")
  }
}

#' onclick function to popup genecards webpage
#' 
#' @title onclick_genecards
#' @param gene query gene
#' @return onclick action
#' @export
#' @author Guangchuang Yu
onclick_genecards <- function(gene) {
  url <- paste0(
    "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
    utils::URLencode(gene, reserved = TRUE)
  )
  paste0('window.open("', url, '", "_blank");')
}
