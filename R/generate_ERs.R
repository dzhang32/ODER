#' Collapse bigwig files into a mean coverage
#'
#' \code{calc_mean_cov} takes as input a set of bws. Then uses the total AUC to
#' normalise coverage per sample and finally collapses the coverage for each
#' sample into a mean coverage.
#'
#' @param bw_paths paths to bigwig files.
#' @param aucs vector containing aucs matching the order of bigwig paths.
#' @param target_auc total auc to normalise all samples to. E.g. 40e6 * 100
#'   would be the estimated total auc for sample sequenced to 40 million reads
#'   of 100bp in length.
#' @param chr_to_filter chrs to obtain mean coverage for.
#'
#' @return list of Rle objects. Each containing the mean coverage for an entire
#'   chromosome.
#'
#' @export
calc_mean_cov <- function(bw_paths, aucs, target_auc, chr_to_filter){

  # chrs need to be in UCSC format for bigwig
  # obtain lengths from UCSC for loadCoverage
  chr_info <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC("hg38") %>%
    filter(UCSC_seqlevel %in% chr_to_filter)

  all_chrs_mean_cov <- list()

  print(stringr::str_c(Sys.time(), " - Obtaining mean coverage across ", length(bw_paths), " samples"))

  for(i in 1:nrow(chr_info)){

    print(stringr::str_c(Sys.time(), " - ", chr_info$UCSC_seqlevel[i]))

    chr_mean_cov <-
      derfinder::loadCoverage(files = bw_paths,
                              totalMapped = aucs, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                              targetSize = target_auc,
                              chr = chr_info$UCSC_seqlevel[i],
                              chrlen = chr_info$UCSC_seqlength[i],
                              inputType = "BigWig",
                              returnMean = T,
                              returnCoverage = F,
                              verbose = F,
                              cutoff = NULL) # setting cutoff as NULLhere and instead to be applied in findRegions()

    all_chrs_mean_cov[[chr_info$UCSC_seqlevel[i]]] <- chr_mean_cov$meanCoverage

  }

  return(all_chrs_mean_cov)

}

#' Define sets of ERs
#'
#' \code{gen_ERs} defines ERs across a inputted range of mean coverage cut-offs
#' (MCCs) and max region gaps (MRGs).
#'
#' @param chrs_mean_cov output of \code{\link{cal_mean_cov}}.
#' @param MCCs mean coverage cut-offs to apply.
#' @param MRGs max region gaps to apply.
#'
#' @return list containing sets of ERs, each generated using a particular
#'   combination of MCC and MRG.
#'
#' @export
gen_ERs <- function(chrs_mean_cov, MCCs, MRGs){

  ##### Create list to save ERs #####

  ERs_MCCs_MRGs <- list()

  for(j in 1:length(MCCs)){

    MCC_label <- stringr::str_c("MCC_", MCCs[j])

    for(k in 1:length(MRGs)){

      MRG_label <- stringr::str_c("MRG_", MRGs[k])

      ERs_MCCs_MRGs[[MCC_label]][[MRG_label]] <- list()

    }
  }

  ##### Generate ERs #####

  for(i in 1:length(chrs_mean_cov)){

    print(stringr::str_c(Sys.time(), " - Generating ERs for ", names(chrs_mean_cov)[i]))

    for(j in 1:length(MCCs)){

      MCC_label <- stringr::str_c("MCC_", MCCs[j])

      # generate ERs at particular MCC
      suppressMessages(
        ERs_MCC <- derfinder::findRegions(position = S4Vectors::Rle(TRUE, length(chrs_mean_cov[[i]])),
                                          fstats = chrs_mean_cov[[i]],
                                          chr = names(chrs_mean_cov)[i],
                                          cutoff = MCCs[j],
                                          maxRegionGap = 0L,
                                          maxClusterGap = length(chrs_mean_cov[[i]]), # setting this to chr length (ignore clusters) to reduce run time. No impact on ERs
                                          verbose = F)
      )

      for(k in 1:length(MRGs)){

        MRG_label <- stringr::str_c("MRG_", MRGs[k])

        # collapse ERs with less than a MRG apart
        ERs_MCCs_MRGs[[MCC_label]][[MRG_label]][[names(chrs_mean_cov)[i]]] <- ERs_MCC %>%
          GenomicRanges::reduce(min.gapwidth = MRGs[k])

      }
    }
  }

  ##### Merge ERs across chromosomes #####

  for(j in 1:length(MCCs)){

    MCC_label <- stringr::str_c("MCC_", MCCs[j])

    for(k in 1:length(MRGs)){

      MRG_label <- stringr::str_c("MRG_", MRGs[k])

      ERs_MCCs_MRGs[[MCC_label]][[MRG_label]] <-
        ERs_MCCs_MRGs[[MCC_label]][[MRG_label]] %>%
        GenomicRanges::GRangesList() %>%
        unlist() %>%
        sort()

    }
  }

  return(ERs_MCCs_MRGs)

}
