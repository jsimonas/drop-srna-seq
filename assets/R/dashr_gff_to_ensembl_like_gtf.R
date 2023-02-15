# libraries
library(tidyverse)
library(rtracklayer)

# download files
download.file(
  url = "https://dashr2.lisanwanglab.org/downloads/dashr.v2.sncRNA.annotation.hg38.gff", 
  destfile = "../gtfs/database/dashr/dashr_v2_hsa.gff3"
)

# read gff 
gff_dash <- rtracklayer::import(
  "../gtfs/database/dashr/dashr_v2_hsa.gff3"
  )

##
## re-format dashr gff annotations to ensembl like .gtf
##

# tabularise granges obj
gff_dash_tb <- gff_dash %>% 
  as_tibble() %>% 
  filter(!grepl("_", seqnames)) %>%
  droplevels()

##
## tRFs

# make tRF gtf
gtf_dash_trf <- gff_dash_tb %>% 
  filter(type %in% c("tRF5","tRNA","tRF3")) %>%
  select(seqnames:ID) %>% 
  # add extra row for each tRF to follow:
  # gene, (transcript, exon)x3
  mutate(
    dumy = rep(c(3,2,2), nrow(.)/3)
  ) %>%
  uncount(dumy) %>% 
  # assign type
  mutate(
    ID = make.unique(ID, sep = "x"),
    type = case_when(
      grepl("tRF5$", ID) ~ "gene",
      grepl("tRF5x1$", ID) ~ "transcript",
      grepl("x1$", ID) ~ "exon",
      grepl("x2$", ID) ~ "exon",
      TRUE ~ "transcript"
    ),
    ID = gsub("x.$", "", ID),
    gene_id = gsub("-tRF5|-tRF3", "", ID),
    gene_name = gene_id
  ) %>% 
  # add gene coordinates
  group_by(gene_id) %>% 
  mutate(
    start = case_when(
      type == "gene" ~ min(start),
      TRUE ~ start
      ),
    end = case_when(
      type == "gene" ~ max(end),
      TRUE ~ end
    ), 
    width = case_when(
      type == "gene" ~ (max(end)-min(start))+1L,
      TRUE ~ width
    )
  ) %>% 
  ungroup() %>% 
  # add missing info
  mutate(
    source = "dashr_v2",
    gene_biotype = "tRNA",
    transcript_id = ifelse(type != "gene", ID, NA),
    transcript_biotype = ifelse(type != "gene", "tRNA", NA),
    transcript_name = ifelse(type != "gene", ID, NA),
    exon_id = ifelse(type == "exon", paste0(transcript_id, ".e"), NA)
  ) %>% 
  # select columns
  select(seqnames:type, gene_biotype, gene_id:exon_id) %>% 
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
  )

##
## piRNA

# make piRNA gtf
gtf_dash_pirna <- gff_dash_tb %>% 
  filter(type %in% "piRNA") %>%
  select(seqnames:ID, ncbiAccession) %>% 
  # modify duplicated piRNA ids and names
  mutate(
  ID = make.unique(ID, sep = "_"),
  ncbiAccession =  make.unique(ncbiAccession, sep = "_"),
  ) %>% 
  # add rows for each piRNA to follow:
  # gene, transcript, exon
  mutate(
    dumy = 3
  ) %>%
  uncount(dumy) %>% 
  # assign type
  mutate(
    type = rep(c("gene", "transcript", "exon"), nrow(.)/3)
  ) %>% 
  # add missing info
  mutate(
    source = "dashr_v2",
    gene_id = ncbiAccession,
    gene_name = gsub("piR", "pir", ID),
    gene_biotype = "piRNA",
    transcript_id = ifelse(type != "gene", paste0(ncbiAccession,".t"), NA),
    transcript_biotype = ifelse(type != "gene", "tRNA", NA),
    transcript_name = ifelse(type != "gene", ID, NA),
    exon_id = ifelse(type == "exon", paste0(ncbiAccession, ".e"), NA)
  ) %>% 
  # select columns
  select(-score:-ncbiAccession) %>% 
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
    )

# combine
gtf_dash <- c(gtf_dash_pirna, gtf_dash_trf)

# export gtf
rtracklayer::export(
  gtf_dash, "../gtfs/database/dashr/dashr_v2_hsa_trf_pirna.gtf"
)

