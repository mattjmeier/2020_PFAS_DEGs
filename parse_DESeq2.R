library(tidyverse)
data_dir <- "~/shared/data/2020_PFAS_files_for_anthony/" # set as required

tsv_files <- fs::dir_ls(data_dir, regexp = "\\.txt$", recurse=1)
test <- readr::read_tsv(tsv_files[1], col_names=T)

DESeq2Results <- tsv_files %>%
  map_dfr(read_tsv, col_names=T, .id="source") %>%
  filter(padj < 0.05) %>%
  mutate(contrast = stringr::str_replace(source, ".txt", "")) %>%
  mutate(contrast = stringr::str_replace(contrast, ".*\\/d", "")) %>%
  mutate(contrast = stringr::str_replace(contrast, "vsDMSO", "")) %>%
  dplyr::select(-source,-myfileName) %>%
  separate(contrast, into=c("timepoint","group"), sep="/") %>%
  extract(group, into=c("chemical","dose"), regex="([[:alpha:]]+)([[:digit:]\\.]+)") %>%
  mutate(group = paste0(timepoint,"_",chemical,"_",dose))

DEG_summary <- DESeq2Results %>% group_by(group) %>% dplyr::tally()

genes_per_group <- DESeq2Results %>% group_by(group) %>% dplyr::select(geneSymbol)

write.table(DEG_summary,
            "~/shared/data/2020_PFAS_files_for_anthony/DEG_summary.txt",
            sep="\t",
            quote=F,
            row.names=F)

write.table(genes_per_group,
            "~/shared/data/2020_PFAS_files_for_anthony/genes_per_group.txt",
            sep="\t",
            quote=F,
            row.names=F)
