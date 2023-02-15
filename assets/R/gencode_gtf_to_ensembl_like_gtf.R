# libraries
library(tidyverse)
library(rtracklayer)

# download files
download.file( 
  url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz",
  destfile = "../gtfs/database/gencode/gencode.v41.primary_assembly.annotation.gtf.gz",
  method = "curl"
)

# read gff 
gtf_gencode <- rtracklayer::import(
  "../gtfs/database/gencode/gencode.v41.primary_assembly.annotation.gtf.gz"
  )

##
## re-format dashr gff annotations to ensembl like .gtf
##

# tabularise granges obj
gtf_gencode_tb <- gtf_gencode %>% 
  as_tibble()

# modify 
gtf_gencode_mod <- gtf_gencode_tb %>% 
  # remove readthrough and PAR transcripts
  filter(!tag %in% c("readthrough_transcript","PAR")) %>%
  # select relevant columns
  select(
    seqnames:gene_name, transcript_id:transcript_name, exon_id
    ) %>% 
  # rename type to biotype
  dplyr::rename(
    gene_biotype = "gene_type",
    transcript_biotype = "transcript_type"
  ) %>% 
  # remove unwanted biotypes
  filter(
    !gene_biotype %in% c("miRNA", "artifact", "TEC")
    ) %>%
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
  )

# export gtf
rtracklayer::export(
  gtf_gencode_mod, "../gtfs/database/gencode/gencode_v41_hsa.gtf"
)
# gzip
R.utils::gzip(
  "../gtfs/database/gencode/gencode_v41_hsa.gtf",
  overwrite=TRUE
)

