# libraries
library(tidyverse)
library(rtracklayer)

# download files
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz", 
  destfile = "../gtfs/database/ncbi/GRCh38_latest_genomic.gff.gz",
  method = "curl"
)

# read gff 
gff_ncbi <- rtracklayer::import(
  "../gtfs/database/ncbi/GRCh38_latest_genomic.gff.gz"
  )

##
## re-format ncbi gff annotations to ensembl like .gtf
##

# tabularise granges obj
gff_ncbi_tb <- gff_ncbi %>% 
  as_tibble()

##
## rRNAs

# get rRNA genes
gtf_ncbi_rrna <- gff_ncbi_tb %>% 
  filter(
    gene_biotype %in% "rRNA" & source %in% "BestRefSeq"
    ) %>%
  droplevels() %>% 
  select(
    seqnames:Dbxref,
    description:gene_biotype
    ) %>% 
  rowwise() %>% 
  mutate(
    Dbxref=gsub("GeneID:","",unlist(Dbxref)[1])
  )

# get chromosomes for rRNAs
gtf_ncbi_rrna_chr <- gff_ncbi_tb %>% 
  filter(!is.na(chromosome)) %>%
  select(
    seqnames, chromosome
  ) %>% 
  distinct() %>%
  # combine with rRNA information
  left_join(
    gtf_ncbi_rrna, .
  ) %>% 
  # remove unknown chromosomes
  filter(
    chromosome != "Unknown" & !is.na(chromosome)
  ) %>% 
  # rename seqnames
  mutate(
    seqnames = factor(
      paste0("chr", chromosome)
      )
  ) %>% 
  # add extra row for each rRNA to follow:
  # gene, (transcript, exon)x3
  mutate(dumy = 3) %>%
  uncount(dumy) %>% 
  # assign type
  mutate(
    ID = make.unique(ID, sep = "x"),
    type = case_when(
      grepl("x1$", ID) ~ "transcript",
      grepl("x2$", ID) ~ "exon",
      TRUE ~ "gene"
      )
    ) %>% 
  dplyr::rename(
    gene_id = "Dbxref",
    gene_name = "gene"
  ) %>%
  # add missing info
  mutate(
    transcript_id = ifelse(type != "gene", paste0(gene_id, ".t"), NA),
    transcript_biotype = ifelse(type != "gene", gene_biotype, NA),
    transcript_name = ifelse(type != "gene", gene_name, NA),
    exon_id = ifelse(type == "exon", paste0(transcript_id, ".e"), NA)
  ) %>% 
  # select columns
  select(seqnames:type, gene_id, gene_biotype, gene_name:exon_id) %>% 
  select(-chromosome) %>% 
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
  )

# export gtf
rtracklayer::export(
  gtf_ncbi_rrna_chr, "../gtfs/database/ncbi/ncbi_v110_rrna.gtf"
)

