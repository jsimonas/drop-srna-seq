# load libraries
library(MetBrewer)
library(tidyverse)
library(ggsankey)

# data
tb <- tibble(
  input = c(rep("bcl-files", 3), rep(NA_character_,3)),
  bcl2fastq = c("R1.fastq", "R2.fastq", "R3.fastq", rep(NA_character_,3)),
  seqkit = c(rep("bc.fastq", 2), rep("cdna.fastq",2), NA_character_, NA_character_),
  cutadapt = c(rep("bc.fastq", 2), rep("cdna.trimmed.fastq",2), "miR.gtf", NA_character_),
  star = c("miR-bam", rep("CBtagged.bam", 2),  rep("miR-bam",2), NA_character_),
  `umi-tools` = c("miR.mtx", rep("demuxUMI.bam", 2), "miR.features", "miR.barcodes", NA_character_),
  samtools = c(NA_character_, "CB.bam", "CB.fastq", rep(NA_character_,2), "annotations.gft"),
  mgcount = c("miR.gff", "CB.mtx", rep(NA_character_,3), "CB.mtx"),
  miraligner = c("CB.mtx", rep(NA_character_,3), "CB.mtx", NA_character_),
  ) %>% 
  make_long(input:miraligner) %>% 
  mutate(
    next_x = as.character(next_x),
    next_x = case_when(
      x == "samtools" & node == "CB.fastq" ~ "miraligner",
      TRUE ~ next_x
    ),
    next_node = case_when(
      x == "samtools" & node == "CB.fastq" ~ "CB.mtx",
      TRUE ~ next_node
    ),
    next_x = factor(next_x, levels = levels(x))
  ) %>% 
  filter(!is.na(node))

#figure
p1 <- ggplot(
  tb,
  aes(
    x = x, next_x = next_x, node = node,
    next_node = next_node, fill = node, label = node
    )
  ) +
  geom_sankey(
    flow.alpha = .45, space = 1.8,
    na.rm = TRUE, node.color = "gray30"
    ) +
  geom_sankey_text(
    size = 3.5, space = 1.8, hjust = 0, 
    position = position_nudge(x = 0.1), na.rm = TRUE
  ) +
#  scale_fill_viridis_d(option = 4, direction = -1)  +
  scale_fill_manual(values = met.brewer("Tam", n=19)) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major.x = element_line(color = "gray60", linetype = 2)
  )

ggsave(p1, filename = "../misc/pipe_scheme.png", width = 9, height = 2.5, dpi = 600)

