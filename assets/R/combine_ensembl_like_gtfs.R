# libraries
library(tidyverse)
library(rtracklayer)

##
## read modified gtfs
##

gtf_mirbase <- rtracklayer::import(
  "../gtfs/database/mirbase/mirbase_v22_hsa_gene_level.gtf"
)

gtf_dash <- rtracklayer::import(
  "../gtfs/database/dashr/dashr_v2_hsa_trf_pirna.gtf"
  )

gtf_gencode <- rtracklayer::import(
  "../gtfs/database/gencode/gencode_v41_hsa.gtf.gz"
)

gtf_ncbi <- rtracklayer::import(
  "../gtfs/database/ncbi/ncbi_v110_rrna.gtf"
)

##
## remove overlapping rRNAs
## from gencode
##

# get unique
gtf_ncbi_rrna_unique <- IRanges::subsetByOverlaps(
  gtf_ncbi, gtf_gencode,
  invert = TRUE
  ) %>%
  as_tibble()

# get overlaps
gtf_gencode_rrna_overlaps <- IRanges::subsetByOverlaps(
  gtf_gencode, gtf_ncbi,
  invert = FALSE
  )

# remove overlaps from gencode
gtf_gencode_filtered <- IRanges::subsetByOverlaps(
    gtf_gencode, 
    gtf_gencode_rrna_overlaps,
    invert = TRUE
  )

##
## combine gtfs to a single file
##

# combine and sort
gtf_comprehensive <- c(
  gtf_mirbase,
  gtf_dash,
  gtf_gencode_filtered,
  gtf_ncbi
  ) %>% 
  sortSeqlevels() %>% 
  sort()

# export gtf
rtracklayer::export(
  gtf_comprehensive,
  "../gtfs/database/comprehensive_annotation.gtf"
)
# gzip
R.utils::gzip(
  "../gtfs/database/comprehensive_annotation.gtf",
  overwrite = TRUE
  )

##
## make biotype count round table for MGcount
##

small_biotypes <- c(
  "snRNA", "snoRNA", "scaRNA", "pre_miRNA", "miRNA",
  "Y_RNA", "siRNA", "vault_RNA", "ribozyme", "misc_RNA",
  "piRNA", "scRNA", "sRNA", "rRNA", "rRNA_pseudogene",
  "Mt_rRNA", "tRNA", "tRF", "tRF5", "tRF3", "Mt_tRNA"
)

btypes_crounds <- gtf_comprehensive %>% 
  as_tibble() %>% 
  group_by(gene_biotype) %>% 
  tally() %>% 
  dplyr::rename(
    biotype = "gene_biotype",
    assignation_round = "n"
  ) %>% 
  mutate(
    assignation_round = case_when(
      biotype %in% {{small_biotypes}} ~ "small",
      TRUE ~ "long"
    )
  )

# btypes_crounds.csv
write_csv(
  btypes_crounds,
  "../gtfs/database/btypes_crounds.csv"
)



