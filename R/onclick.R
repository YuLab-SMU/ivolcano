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

#' onclick function to popup NCBI gene webpage
#'
#' @title onclick_ncbi
#' @param gene query gene
#' @return onclick action
#' @export
#' @author Guangchuang Yu
onclick_ncbi <- function(gene) {
  url <- paste0(
    "https://www.ncbi.nlm.nih.gov/gene/?term=",
    utils::URLencode(gene, reserved = TRUE)
  )
  paste0('window.open("', url, '", "_blank");')
}

#' onclick function to popup Ensembl gene webpage
#'
#' @title onclick_ensembl
#' @param gene query gene
#' @return onclick action
#' @export
#' @author Guangchuang Yu
onclick_ensembl <- function(gene) {
  url <- paste0(
    "https://www.ensembl.org/Multi/Search/Results?q=",
    utils::URLencode(gene, reserved = TRUE),
    ";site=ensembl"
  )
  paste0('window.open("', url, '", "_blank");')
}

#' onclick function to popup HGNC gene webpage
#'
#' @title onclick_hgnc
#' @param gene query gene
#' @return onclick action
#' @export
#' @author Guangchuang Yu
onclick_hgnc <- function(gene) {
  url <- paste0(
    "https://www.genenames.org/tools/search/#!/?query=",
    utils::URLencode(gene, reserved = TRUE)
  )
  paste0('window.open("', url, '", "_blank");')
}

#' onclick function to popup UniProt webpage
#'
#' @title onclick_uniprot
#' @param gene query gene
#' @return onclick action
#' @export
#' @author Guangchuang Yu
onclick_uniprot <- function(gene) {
  url <- paste0(
    "https://www.uniprot.org/uniprot/?query=",
    utils::URLencode(gene, reserved = TRUE),
    "&sort=score"
  )
  paste0('window.open("', url, '", "_blank");')
}

#' onclick function to popup PubMed webpage
#'
#' @title onclick_pubmed
#' @param gene query gene
#' @return onclick action
#' @export
#' @author Guangchuang Yu
onclick_pubmed <- function(gene) {
  url <- paste0(
    "https://pubmed.ncbi.nlm.nih.gov/?term=",
    utils::URLencode(gene, reserved = TRUE)
  )
  paste0('window.open("', url, '", "_blank");')
}
