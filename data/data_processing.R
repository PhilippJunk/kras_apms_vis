## 
## Data Processing for KRAS APMS Shiny App
## 

app_data_dir <- '/home/junkpp/work/other/2022-05-20_KRAS_APMS/kras_apms_vis/data/'

# data dir can be obtained from Zenodo: TODO
data_dir <- '/home/junkpp/work/12_CT_APMS/analysis/'

library(tidyverse)

###############################################################################
# Import data

# Diff analysis
df_gsea <- read.csv(str_glue('{data_dir}/outputs_02/df_gsea.csv')) %>%
  filter(p_adj <= 0.1)
df_anova <- read.csv(str_glue('{data_dir}/outputs_03/df_anova.csv')) %>% 
  filter(p_adj <= 0.1)

# LFQi
df_apms <- read.csv(str_glue('{data_dir}/outputs_01/df_apms.csv'))
df_sum <- read.csv(str_glue('{data_dir}/outputs_03/df_sum.csv')) %>%
  filter(id %in% unique(c(df_gsea$id, df_anova$id)))

# GO data
df_ontology <- read.csv(str_glue('{data_dir}/outputs_03/df_ontology.csv')) %>%
  filter(id %in% unique(c(df_gsea$id, df_anova$id)))
df_annotation <- read.csv(str_glue('{data_dir}/outputs_03/df_annotation.csv')) %>%
  filter(id %in% unique(c(df_gsea$id, df_anova$id))) %>%
  mutate(definition = replace_na(definition, ''))

# GO semantic analysis
load(str_glue('{data_dir}/outputs_04/dist_mat.Rdata'))
load(str_glue('{data_dir}/outputs_04/clusters.Rdata'))
df_wordcloud <- read.csv(str_glue('{data_dir}/outputs_04/df_wordcloud.csv')) %>%
  filter(padj < 0.05)

###############################################################################
# save output
save(df_apms, df_sum, df_ontology, df_annotation, df_gsea, df_anova, 
     dist_mat, clusters, df_wordcloud,
     file = str_glue('{app_data_dir}/data.Rdata'))
