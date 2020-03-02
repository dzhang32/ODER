#'Optain set of non-overlapping exons
#'
#'\code{get_no_overlap_exons} uses as input a gtf to optain set of
#'non-overlapping exons. These can be used in \code{\link{gen_ERs}} as the
#'\code{opt_gr} to optimise ER definitions.
#'
#'@param gtf Either a string if path to a .gtf file or a pre-imported gtf using
#'  \code{\link[rtracklayer]{import}}.
#'@param add_chr logical scalar, determining whether to add "chr" prefix to the
#'  seqnames of non-overlapping exons. Note, if set to TRUE and seqnames already
#'  have "chr", it will not add another.
#'
#'@return GRanges object containing non-overlapping exons.
#'@export
get_no_overlap_exons <- function(gtf, add_chr){

  if(is.character(gtf)){

    print(str_c(Sys.time(), " - Loading in GTF..."))

    gtf_gr <- rtracklayer::import(gtf)

  }else{

    gtf_gr <- gtf

  }

  print(str_c(Sys.time(), " - Obtaining non-overlapping exons"))

  exons_gr <- gtf_gr[gtf_gr$type == "exon"]
  exons_gr <- exons_gr[!duplicated(exons_gr$exon_id)] %>% unique()

  exons_hits <- GenomicRanges::findOverlaps(exons_gr, drop.self = T)

  # all(S4Vectors::queryHits(gtf_gr_exons_hits) %in% S4Vectors::subjectHits(gtf_gr_exons_hits))

  exons_no_overlap_gr <- exons_gr[-c(S4Vectors::queryHits(exons_hits) %>% unique())]

  # check - no overlaps
  # GenomicRanges::findOverlaps(exons_no_overlap_gr, drop.self = T)

  if(add_chr){

    GenomeInfoDb::seqlevels(exons_no_overlap_gr) <-
      GenomeInfoDb::seqlevels(exons_no_overlap_gr) %>%
      stringr::str_replace("chr", "") %>%
      stringr::str_c("chr", .)

  }

  return(exons_no_overlap_gr)

}

#' Calculates delta for sets of ERs
#'
#' \code{calc_ERs_delta} calculates the delta/difference between a set of ERs
#' and another given set of GRanges that are the optimal
#'
#' @param ERs_MCCs_MRGs Sets of ERs across various MCCs/MRGs - output of
#'   \code{\link{gen_ERs}}.
#' @param opt_gr GRanges object that contains the regions that ideally, you want
#'   to the ER definitions to match
#' @param delta_func Function that calculates the delta between ERs and
#'   \code{opt_gr}. Takes as input a set of ERs from \code{ERs_MCCs_MRGs} and
#'   \code{opt_gr}. Then outputs a tibble/dataframe containing the summarised
#'   delta scores for that set of one set of ERs.
#'
#' @return tibble/dataframe containing summarised delta values. One row per set
#'   of ERs.
#'
#' @export
calc_ERs_delta <- function(ERs_MCCs_MRGs, opt_gr, delta_func = .delta){

  print(str_c(Sys.time(), " - Calculating delta for ERs..."))

  MCC_labels <- names(ERs_MCCs_MRGs)

  delta_df <- dplyr::tibble()

  for(i in 1:length(MCC_labels)){

    MRG_labels <- names(ERs_MCCs_MRGs[[MCC_labels[i]]])

    for(j in 1:length(MRG_labels)){

      delta_summarised <-
        delta_func(query = ERs_MCCs_MRGs[[MCC_labels[i]]][[MRG_labels[j]]],
                   subject = opt_gr)

      delta_df <- delta_df %>%
        dplyr::bind_rows(delta_summarised %>%
                           dplyr::mutate(MCC = MCC_labels[i],
                                         MRG = MRG_labels[j]))

    }
  }

  delta_df <- delta_df %>%
    dplyr::select(MCC, MRG, dplyr::everything())

  return(delta_df)

}

#' Optains optimised set of ERs
#'
#' \code{get_opt_ERs} optains
#'
#' @inheritParams calc_ERs_delta
#'
#' @param delta_df Output of \code{\link{calc_ERs_delta}}.
#'
#' @return list containing optimised ERs, optimal pair of MCC/MRGs and
#'   \code{delta_df}
#'
#' @export
get_opt_ERs <- function(ERs_MCCs_MRGs, delta_df){

  print(str_c(Sys.time(), " - Obtaining optimal set of ERs..."))

  delta_opt <-
  delta_df %>%
    dplyr::filter(median == min(median)) %>% # with the lowest median ER delta
    dplyr::filter(n_eq_0 == max(n_eq_0)) # and highest num of delta equal to 0

  opt_ERs <-
      list(opt_ERs = ERs_MCCs_MRGs[[delta_opt$MCC]][[delta_opt$MRG]],
           opt_MCC_MRG = c(delta_opt$MCC, delta_opt$MRG),
           deltas = delta_df)

  return(opt_ERs)

}

# function to obtain delta between ERs + exons
# can be swapped for user-defined functions
.delta <- function(query, subject){

  hits <- GenomicRanges::findOverlaps(query = query, subject = subject)

  # obtain situation where 1 ER overlaps multiple exons...
  n_dis_exons_ab_1 <-
    hits %>%
    as.data.frame() %>%
    dplyr::group_by(queryHits) %>%
    dplyr::summarise(n_dis_exons = dplyr::n_distinct(subjectHits)) %>%
    dplyr::filter(n_dis_exons > 1)

  # ...and remove them
  hits <- hits[!(S4Vectors::queryHits(hits) %in% n_dis_exons_ab_1$queryHits)]

  delta_raw <-
    dplyr::bind_cols(query[S4Vectors::queryHits(hits)] %>%
                       as.data.frame(row.names = NULL)%>%
                       dplyr::select(seqnames, start, end),
                     subject[S4Vectors::subjectHits(hits)] %>%
                       as.data.frame(row.names = NULL) %>%
                       dplyr::select(seqnames, start, end)) %>%
    dplyr::mutate(start_diff = start - start1,
                  end_diff = end - end1,
                  delta = abs(start_diff) + abs(end_diff))

  delta_summarised <-
    dplyr::tibble(sum = sum(delta_raw$delta),
                  mean = mean(delta_raw$delta),
                  median = median(delta_raw$delta),
                  n_eq_0 = sum(delta_raw$delta == 0),
                  propor_eq_0 = mean(delta_raw$delta == 0))

  return(delta_summarised)

}
