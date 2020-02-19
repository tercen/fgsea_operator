library(tercen)
library(dplyr)
library(fgsea)

ctx = tercenCtx()

gene_factor_name_char <- ctx$rnames[[1]]
genes <- ctx$rselect(ctx$rnames)
genes_list <- as.character(unlist(genes))

pathways_df <- ctx$cselect(ctx$cnames)
pathways_df <- bind_cols(pathways_df, col_idx = (0:(nrow(pathways_df)-1)))
names(pathways_df) <- c("pathway", ".ci")

in_matrix <- as.tibble(ctx$as.matrix())
values <- (in_matrix != 0)
values <- unnest(as.tibble(values))


# values <- bind_cols(genes, in_table)
in_table <- ctx$select()

gene_df <- in_table %>% 
  select(.y, .ci, .ri) %>% 
  group_by(.ri) %>%
  summarise(rank = mean(.y)) %>%
  bind_cols(genes) 

rank_list <- gene_df$rank
names(rank_list) <- gene_df[[gene_factor_name_char]]

getGeneId <- function(alist) genes_list[alist]

pathway_list <- lapply(values, getGeneId)


names(pathway_list) <- unlist(pathways_df$pathway)

fgseaRes <- fgsea(pathway_list, rank_list, maxSize=500)


# fgseaRes <- fgsea(examplePathways, exampleRanks, maxSize=500)


result <- fgseaRes %>% select(pathway, NES) %>% merge(pathways_df, by = c("pathway"))


result %>%
  ctx$addNamespace() %>%
  ctx$save()

