#' Optimising the Definition of Expressed Regions
#'
#' \code{ODER} is a wrapper function that runs the pipeline from generation of ERs
#' across different MCCs/MRGs to selecting the optimal set of ERs. Either a gtf
#' is required, to optain a set of non-overlapping exons to use as the
#' \code{opt_gr} OR users can provide their own set of ranges to optimise ER
#' definitions towards.
#'
#' @inheritParams calc_mean_cov
#' @inheritParams gen_ERs
#' @inheritParams get_no_overlap_exons
#' @inheritParams calc_ERs_delta
#' @inheritParams get_opt_ERs
#'
#' @return list containing optimised ERs, optimal pair of MCC/MRGs and
#'   \code{delta_df}
#'
#' @export
ODER <- function(bw_paths, aucs, target_auc, chr_to_filter = stringr::str_c("chr", c(1:22, "X", "Y", "M")),
                 MCCs, MRGs,
                 gtf = NULL, ucsc_chr, ignore.strand,
                 opt_gr = NULL){

  if(is.null(gtf) && is.null(opt_gr)) stop("One of gtf OR opt_gr must be provided.")

  chrs_mean_cov <-
    calc_mean_cov(bw_paths, aucs, target_auc, chr_to_filter)

  ERs_MCCs_MRGs <- gen_ERs(chrs_mean_cov, MCCs, MRGs)

  if(!is.null(gtf)){

    opt_gr <- get_no_overlap_exons(gtf, ucsc_chr, ignore.strand)

  }

  delta_df <- calc_ERs_delta(ERs_MCCs_MRGs, opt_gr)

  opt_ERs <- get_opt_ERs(ERs_MCCs_MRGs, delta_df)

  return(opt_ERs)

}
