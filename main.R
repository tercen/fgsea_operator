library(tercen)
library(dplyr)
library(fgsea)
library(tidyverse)


ctx = tercenCtx()

pathway_factor_name_char <- ctx$rnames[[1]]
pathway_df <- ctx$rselect(ctx$rnames) %>% mutate(.ri = 0:(nrow(.)-1))
gene_factor_name <- sym(ctx$labels[[1]])
gene_factor_char <- ctx$labels[[1]]

in_table <- ctx$select(list(".ci", ".ri", ".y", ctx$labels[[1]]))

gene_df <- in_table %>% 
  select(.y, !!gene_factor_name) %>%
  group_by(!!gene_factor_name) %>%
  summarise(value = mean(.y)) %>%
  column_to_rownames(ctx$labels[[1]])

rank_list <- gene_df$value
names(rank_list) <- rownames(gene_df)

getGeneId <- function(df,...){as.character(unlist((df[[gene_factor_char]])))}
pathway_list <- in_table %>% group_by(.ri) %>% group_split %>% map(getGeneId)
pathway_names <- pathway_df[[pathway_factor_name_char <- ctx$rnames[[1]]]]
names(pathway_list) <- pathway_names

fgseaRes <- fgsea(pathway_list, rank_list, maxSize=500, nperm=10000)

result <- fgseaRes %>% select(pathway, NES, padj) %>% left_join(pathway_df, by = c("pathway"= pathway_factor_name_char))

result %>%
  ctx$addNamespace() %>%
  ctx$save()
