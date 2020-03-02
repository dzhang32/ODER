library(tidyverse)
library(stringr)
# devtools::install_github("dzhang32/ODER")
# library(ODER)
devtools::load_all("~/packages/ODER")

# Load data -------------------------------------------------------------------------------------------

gtex_recount_meta_data_tidy <- read_delim("/data/recount/GTEx_SRP012682/gtex_recount_meta_data_tidy.txt", delim = "\t")

load("~/projects/OMIM_wd/results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_list.rda")

ens_gtf <- rtracklayer::import("/data/references/ensembl/gtf_gff3/v87/Homo_sapiens.GRCh38.87.gtf")
ens_TxDb <- AnnotationDbi::loadDb(file = "/data/references/ensembl/txdb_sqlite/v87/ensembl_grch38_v87_txdb.sqlite")

# Functions -------------------------------------------------------------------------------------------

mark_overlapping_genes_gr <- function(gr_1, gr_2, identical_gr = F, maxgap = -1L, minoverlap = 1L, ignore.strand = T, ...){

  gr_1_overlapping_counts <-
    gr_1 %>%
    countOverlaps(gr_2, maxgap = maxgap, minoverlap = minoverlap, ignore.strand = ignore.strand, ...)

  if(identical_gr == T){

    stopifnot(all(gr_1_overlapping_counts >= 1))

    gr_1$overlap_gr2 <-
      gr_1_overlapping_counts > 1

  }else{

    gr_1$overlap_gr2 <-
      gr_1_overlapping_counts > 0

  }

  num_overlapping_genes <- sum(gr_1$overlap_gr2)

  print(str_c(num_overlapping_genes, "/", length(gr_1), " (propor: ", round(num_overlapping_genes/length(gr_1), digits = 2), ")", " overlapping ranges.."))

  return(gr_1)

}

# Main ------------------------------------------------------------------------------------------------

##### Check non-overlapping exons match #####

exons_txdb <-
  ens_TxDb %>% exons(columns = c("EXONNAME", "GENEID", "TXNAME"))
exons_txdb <- mark_overlapping_genes_gr(exons_txdb, exons_txdb, identical_gr = T)
exons_txdb_no_overlap <- exons_txdb[exons_txdb$overlap_gr2 == F]

exons_gtf_no_overlap <- ODER::get_no_overlap_exons(gtf = ens_gtf, add_chr = T, ignore.strand = T)
exons_gtf_no_overlap <- exons_gtf_no_overlap %>% keepSeqlevels(str_c("chr", c(1:22, "X", "Y", "M")), pruning.mode = "coarse")

identical(ranges(sort(exons_txdb_no_overlap)), ranges(sort(exons_gtf_no_overlap)))

##### Check optimal ERs match #####

gtex_metadata_crbl <-
  gtex_recount_meta_data_tidy %>%
  filter(smtsd == "Brain - Cerebellum",
         smafrze == "USE ME")

opt_ERs_crbl <- ODER(bw_paths = str_c("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/",
                                      gtex_metadata_crbl$bigwig_file),
                     aucs = gtex_metadata_crbl$auc,
                     target_auc = (40e6 * 100),
                     chr_to_filter = c("chr21", "chr22"), # stringr::str_c("chr", c(1:22, "X", "Y", "M")),
                     MCCs = seq(1, 10, 0.2),
                     MRGs = seq(0, 100, by = 10),
                     gtf = ens_gtf,
                     ucsc_chr = T,
                     ignore.strand = T)


