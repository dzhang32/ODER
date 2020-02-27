#' Collapse bigwig files into a mean coverage
#'
#' \code{cal_mean_cov} takes as input a set of bws. Then uses the total AUC to normalise coverage per sample and finally collapses the coverage for each sample into a mean coverage.
#'
#' @param bw_paths paths to bigwig files.
#' @param aucs vector containing aucs matching the order of bigwig paths.
#' @param target_auc total auc to normalise all samples to. E.g. 40e6 * 100 would be the estimated total AUC for sample sequenced to 40 million reads of 100bp in length.
#' @param chr_to_filter chrs to obtain mean coverage for. Default: chr1-22, chrX, chrY, chrM.
#'
#' @return list of Rle objects. Each containing the mean coverage for an entire chromosome.
#'
#' @export
cal_mean_cov <- function(bw_paths, aucs, target_auc, chr_to_filter = stringr::str_c("chr", c(1:22, "X", "Y", "M"))){

  # chrs need to be in UCSC format for bigwig, could be changed to be more flexible
  chr_info <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC("hg38") %>%
    filter(UCSC_seqlevel %in% chr_to_filter)

  all_chrs_mean_cov <- list()

  for(i in 1:nrow(chr_info)){

    print(stringr::str_c(Sys.time(), " - obtaining mean cov for ", chr_info$UCSC_seqlevel[i]))

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
                              cutoff = NULL) # setting cutoff as null here and instead will be applied downstream in findRegions()

    all_chrs_mean_cov[[chr_info$UCSC_seqlevel[i]]] <- chr_mean_cov$meanCoverage

  }

  return(all_chrs_mean_cov)

}

#' Generate sets of ERs
#'
#' \code{gen_ERs} defines ERs across a inputted range of mean coverage cut-offs (MCCs) and max region gaps (MRGs).
#'
#' @param chr_mean_cov output of \code{cal_mean_cov}.
#' @param MCCs cut-offs to apply.
#' @param MRGs max region gaps to apply.
#'
#' @return list containing sets of ERs, each generated using a particular combination of MCC and MRG.
#'
#' @export
def_ERs <- function(chr_mean_cov, MCCs, MRGs){

  ##### Create list to save ERs #####

  ERs_MCCs_MRGs <- list()

  for(j in seq_along(MCCs)){

    MCC_label <- stringr::str_c("MCC_", MCCs[j])

    for(k in seq_along(MRGs)){

      MRG_label <- stringr::str_c("MRG_", MRGs[k])

      ERs_MCCs_MRGs[[MCC_label]][[MRG_label]] <- list()

    }
  }

  ##### Generate ERs #####

  for(i in 1:length(chrs_mean_cov)){

    for(j in seq_along(MCCs)){

      MCC_label <- stringr::str_c("MCC_", MCCs[j])

      # generate ERs at particular MCC
      suppressMessages(
        ERs_MCC <-derfinder::findRegions(position = Rle(TRUE, length(chrs_mean_cov[[i]])),
                                         fstats = chrs_mean_cov[[i]],
                                         chr = names(chrs_mean_cov)[i],
                                         cutoff = MCCs[j],
                                         maxRegionGap = 0L,
                                         maxClusterGap = length(chrs_mean_cov[[i]]), # setting this to chr length (ignore clusters) to reduce run time. No impact on ERs
                                         verbose = T)
      )

      for(k in seq_along(MRGs)){

        MRG_label <- stringr::str_c("MRG_", MRGs[k])

        print(stringr::str_c(Sys.time(), " - Generating ERs for ", names(chrs_mean_cov)[i], " - MCC:", MCCs[j], ", MRG:", MRGs[k]))

        ERs_MCCs_MRGs[[MCC_label]][[MRG_label]][[names(chrs_mean_cov)[i]]] <- ERs_MCC %>%
          GenomicRanges::reduce(min.gapwidth = MRGs[k])

      }
    }
  }

  ##### Merge ERs across chromosomes #####

  for(j in seq_along(MCCs)){

    MCC_label <- stringr::str_c("MCC_", MCCs[j])

    for(k in seq_along(MRGs)){

      MRG_label <- stringr::str_c("MRG_", MRGs[k])

      ERs_MCCs_MRGs[[MCC_label]][[MRG_label]] <-
        ERs_MCC_MRGs[[MCC_label]][[MRG_label]] %>%
        GenomicRanges::GRangesList() %>%
        unlist() %>%
        sort()

    }
  }

  return(ERs_MCCs_MRGs)

}

optimise_ERs <- function(ERs_MCCs_MRGs, gr){



}

get_no_overlap_exons <- function(gtf_path){

  gtf_path <- "/data/re"


}
