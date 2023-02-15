# libraries
library(tidyverse)
library(rtracklayer)

# download files
download.file(
  url = "https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3", 
  destfile = "../gtfs/database/mirbase/mirbase_v22_hsa.gff3"
  )

download.file(
  url = "https://mirgenedb.org/gff/hsa?sort=pos&all=1", 
  destfile = "../gtfs/database/mirgenedb/mirgenedb_v2_hsa.gff3"
)

# read mirbase gff 
gff_mirbase <- rtracklayer::import(
  "../gtfs/database/mirbase/mirbase_v22_hsa.gff3"
  )
# read mirgenedb gff 
gff_mirgenedb <- rtracklayer::import(
  "../gtfs/database/mirgenedb/mirgenedb_v2_hsa.gff3"
)

##
## re-format mirbase gff annotations to ensembl like .gtf
##

##
## gene level

# tabularise granges obj
gff_mirbase_tb <- gff_mirbase %>% 
  as_tibble() %>% rowid_to_column()

# make gtf
gtf_mirbase <- gff_mirbase_tb %>% 
  filter(type == "miRNA_primary_transcript") %>% 
  # modify duplicated pre-mirna ids and names
  mutate(
    ID = make.unique(ID, sep = "_"),
    Name =  make.unique(Name, sep = "-"),
  ) %>% 
  bind_rows(., anti_join(
    gff_mirbase_tb,., by = "rowid"
  )) %>%
  arrange(rowid) %>%
  # modify duplicated mature mirna names
  mutate(
    Name =  make.unique(Name, sep = "_")
  ) %>% 
  # add extra row for each mature
  mutate(
    dumy = if_else(type == "miRNA", 2,1)
  ) %>%
  uncount(dumy) %>% 
  # assign type
  mutate(
    Name = make.unique(Name, sep = "x"),
    type = case_when(
      type == "miRNA_primary_transcript" ~ "gene",
      grepl("x1", Name) ~ "exon",
      TRUE ~ "transcript"
    ),
    Name = gsub("x1", "", Name)
  ) %>% 
  # add missing info
  mutate(
    source = "mirbase_v22.1",
    gene_id = ifelse(grepl("MIM", ID), Derives_from, ID),
    gene_name = ifelse(!(grepl("-miR-", Name) | grepl("-let-.*[0-9]p", Name)), Name, NA),
    gene_biotype = "pre_miRNA",
    transcript_id = ifelse(grepl("MIM", ID), ID, NA),
    transcript_biotype = ifelse(!type %in% "gene", "miRNA", NA),
    transcript_name = ifelse(grepl("-miR-", Name) | grepl("-let-.*[0-9]p", Name), Name, NA),
    exon_id = ifelse(type == "exon", paste0(transcript_id, ".e"), NA)
  ) %>% 
  fill(gene_name) %>% 
  # remove unwanted columns
  select(-rowid, -ID:-Derives_from) %>%
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
  )

# export gtf
rtracklayer::export(
  gtf_mirbase, "../gtfs/database/mirbase/mirbase_v22_hsa_gene_level.gtf"
  )

##
## mature miRNA level
gtf_mirbase_mat <- gtf_mirbase %>% 
  as_tibble() %>% 
  # remove gene rows
  filter(!is.na(transcript_id)) %>% 
  # add row for each mature
  mutate(
    dumy= ifelse(is.na(exon_id),2,1)
  ) %>% 
  uncount(dumy) %>% 
  # modify gene ids and names
  group_by(transcript_id) %>%
  mutate(
    type = make.unique(type, sep = "x"),
    type = case_when(
      !grepl("x1", type) & !type=="exon" ~ "gene",
      TRUE ~ type
    ),
    type = gsub("x1", "", type),
    gene_id = transcript_id,
    gene_name = transcript_name,
    gene_biotype = "mature_miRNA"
  ) %>% 
  mutate(
    transcript_id = case_when(
      type == "gene" ~ NA_character_,
      TRUE ~ paste0(transcript_id, ".t")
    ),
    transcript_name = case_when(
      type == "gene" ~ NA_character_,
      TRUE ~ paste0(transcript_name, "t")
    )
  ) %>% 
  ungroup() %>% 
  # expand mature seq by 3 bases
  mutate(
    start = start-3,
    end = end+3
  ) %>% 
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
  )

# export gtf
rtracklayer::export(
  gtf_mirbase_mat, "../gtfs/database/mirbase/mirbase_v22_hsa_mature_level.gtf"
)


##
## re-format mirgenedb gff annotations to ensembl like .gtf
##

# tabularise granges obj
gff_mirgenedb_tb <- gff_mirgenedb %>% 
  as_tibble() %>% unnest(Alias)

# make gtf
gtf_mirgenedb <- gff_mirgenedb_tb %>% 
  # modify duplicated alias from mirbase
  mutate(
    Alias = make.unique(Alias, sep = ".")
  ) %>% 
  # add extra row for each mature
  mutate(
    dumy = if_else(type == "miRNA", 2,1)
  ) %>%
  uncount(dumy) %>% 
  # assign type
  mutate(
    ID = make.unique(ID, sep = "x"),
    type = case_when(
      type == "pre_miRNA" ~ "gene",
      grepl("x1", ID) ~ "exon",
      TRUE ~ "transcript"
    ),
    ID = gsub("x1", "", ID)
  ) %>% 
  # add missing info
  mutate(
    source = "mirgenedb_v2.1",
    gene_id = ifelse(grepl("MIM", Alias), NA, Alias),
    gene_name = ifelse(grepl("_pre", ID), ID, NA),
    gene_biotype = "pre_miRNA",
    transcript_id = ifelse(grepl("MIM", Alias), Alias, NA),
    transcript_biotype = "miRNA",
    transcript_name = ifelse(!grepl("_pre", ID), ID, NA),
    exon_id = ifelse(type == "exon", paste0(transcript_id, ".e"), NA)
  ) %>% 
  fill(c(gene_id, gene_name)) %>% 
  # remove unwanted columns
  select(-ID:-Alias) %>%
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
  )

# export gtf
rtracklayer::export(
  gtf_mirgenedb,
  "../gtfs/database/mirgenedb/mirgenedb_v2_hsa_gene_level.gtf"
)

##
## mature miRNA level
gtf_mirgenedb_mat <- gtf_mirgenedb %>% 
  as_tibble() %>% 
  # remove gene rows
  filter(!is.na(transcript_id)) %>% 
  # add row for each mature
  mutate(
    dumy= ifelse(is.na(exon_id),2,1)
  ) %>% 
  uncount(dumy) %>% 
  # modify gene ids and names
  group_by(transcript_id) %>%
  mutate(
    type = make.unique(type, sep = "x"),
    type = case_when(
      !grepl("x1", type) & !type=="exon" ~ "gene",
      TRUE ~ type
    ),
    type = gsub("x1", "", type),
    gene_id = transcript_id,
    gene_name = transcript_name,
    gene_biotype = "mature_miRNA"
  ) %>% 
  mutate(
    transcript_id = case_when(
      type == "gene" ~ NA_character_,
      TRUE ~ paste0(transcript_id, ".t")
    ),
    transcript_name = case_when(
      type == "gene" ~ NA_character_,
      TRUE ~ paste0(transcript_name, "t")
    )
  ) %>% 
  ungroup() %>% 
  # expand mature seq by 3 bases
  mutate(
    start = start-3,
    end = end+3
  ) %>% 
  # convert back to granges obj
  makeGRangesFromDataFrame(
    keep.extra.columns=TRUE
  )

# export gtf
rtracklayer::export(
  gtf_mirgenedb_mat,
  "../gtfs/database/mirgenedb/mirgenedb_v2_hsa_mature_level.gtf"
)


