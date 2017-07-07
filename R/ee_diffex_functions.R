################################################################################

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

################################################################################

#' get_diffexed_genes
#'
#' Extract genes from a limma-fit object (resulting from an eBayes call) that
#' are differentially expressed. Two different methods are used for calling
#' differential expression:
#'   \code{any_probe} - at least one probe for that gene must be differentially
#'   expressed.
#'   \code{all_probes} - every probe for that gene must be differentially
#'   expressed (and must be differentially expressed in the same direction).
#' This doesn't return information about the direction of change, and only
#' returns a vector of gene symbols or ids, as extracted from the \code{genes}
#' dataframe in the limma-fit object \code{fit}.
#'
#' @param        fit           A limma fit object (class: MArrayLM) as returned
#'   by eBayes. If the initial model/analysis contained more than one contrast
#'   the user should select a single contrast prior to passing \code{fit} into
#'   \code{get_diffexed_genes}, by doing fit[, some.column] - if info on more
#'   than one contrast is present in the \code{fit} object, the function throws
#'   an error.
#'
#' @param        gene_column   The name of the column within \code{fit$genes}
#'    that contains the gene-symbols or gene-ids that are to be returned by
#'    this function. For example, if you want to receive a subset of the
#'    "entrez.id" column of \code{fit$genes} then specify
#'    \code{gene_column = "entrez.id"}. The function dies if the specified
#'    column is absent from the \code{fit$genes} dataframe.
#'
#' @param        diffex_type   The user selects either \code{"any_probe"} or
#'   \code{"all_probes"} to indicate the type of analysis to be performed.
#'
#' @param        ...           Further arguments that are to be passed through
#'   to limma::decideTests. For example, if you want to set the \code{lfc}
#'   filter to 1 (in decideTests; corresponding to keeping probes only if they
#'   are at least 2-fold up/down-regulated), you can specify it as follows:
#'   \code{get_diffexed_genes(some.fit, some.column, some.type, lfc = 1)}.
#'
#' @return       A subvector of the inidicated 'gene_column' in fit fit$genes
#'   annotation dataframe.
#'
#' @importFrom   magrittr      %>%   extract2
#' @importFrom   limma         decideTests
#' @importFrom   tibble        data_frame
#' @importFrom   dplyr         ......
#'
#' @importClassesFrom   limma   MArrayLM
#'
#' @export
#'
get_diffexed_genes <- function(
    fit,
    gene_column = "symbol",
    diffex_type = c("any_probe", "all_probes"),
    ...
){
  diffex_type <- match.arg(diffex_type)


  stopifnot(
    is.character(gene_column) &&
      length(gene_column) == 1
    )
  stopifnot(
    is(fit, "MArrayLM")           &&
      "genes"     %in% names(fit) &&
      gene_column %in% colnames(fit[["genes"]])
    )
  if (ncol(fit) != 1) {
    stop("fit should contain details for one contrast; try rerunning with
         get_diffexed_genes(fit[, some.column], ...)")
  }

  genes         <- fit[["genes"]][, gene_column]
  dt            <- limma::decideTests(fit, ...)
  diffexed_rows <- which(dt != 0)

  diffexed_genes <- if (
      diffex_type == "any_probe"
    ) {
    genes[diffexed_rows]
  } else {
    tibble::data_frame(
      diffex = dt[, 1],
      gene   = genes
    ) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(diffex = all(diffex == 1) |
                                all(diffex == -1)
                       ) %>%
      dplyr::filter(diffex) %>%
      magrittr::extract2("gene")
  }

  diffexed_genes %>%
    unique %>%
    sort
}
