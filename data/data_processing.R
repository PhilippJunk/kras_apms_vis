## 
## Data Processing for KRAS APMS Shiny App
## Written by Philipp Junk, 2022
## Data by Camille Ternet
## 

library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(MsCoreUtils)
library(limma)

library(tidyverse)
library(furrr)

future::plan(multisession, workers = 7)

# TODO download raw datasets from zenodo to data/raw/ and read from them
if (.Platform$OS.type == 'unix') {
  setwd('/home/junkpp/work/12_CT_APMS')
} else if (.Platform$OS.type == 'windows') {
  setwd("D:/Daten/Work/Dublin/CloudStation/12_CT_APMS")
}


###############################################################
# files and databases data

# TODO put in RAW
# get hgnc database file
if (.Platform$OS.type == 'unix') {
  hgnc_database <- '/home/junkpp/work/12_CT_APMS/id_mapping/hgnc_complete_set.txt'
} else {
  hgnc_database <- 'D:/Daten/Work/Dublin/CloudStation/12_CT_APMS/id_mapping/hgnc_complete_set.txt'
}

###############################################################
# APMS data
df_apms <- read_csv('data/2021-12-14/raw_biolrepl_long.csv') %>% 
  left_join(read.csv('data/annotation.csv') %>%
              mutate(label = str_c(mut_status, condition, concentration, biol_repl, sep='_'),
                     group = str_c(mut_status, condition, concentration, sep='_')) %>%
              select(label, group, mut_status, condition, concentration, biol_repl) %>%
              distinct,
            by = 'label') %>% 
  left_join(read_tsv(hgnc_database, show_col_types = F) %>% 
              select(symbol, entrez_id), 
            by=c('hgnc' = 'symbol')) %>%
  mutate(entrez_id = as.character(entrez_id)) %>%
  mutate(mut_status = case_when(mut_status == '12c' ~ 'g12c',
                                mut_status == '12d' ~ 'g12d',
                                mut_status == '12v' ~ 'g12v',
                                T ~ mut_status)) %>%
  mutate(group = str_c(mut_status, condition, concentration, sep='_'),
         label = str_c(mut_status, condition, concentration, biol_repl, sep='_')) %>%
  select(-biol_repl)



###############################################################
# Annotate with ontology data

# Ontologies used:
# - SysGO
# - GO Biological Processes
# 
# For each ontology, all gene sets are extracted and checked whether
# their components are present in the APMS. Only gene sets with 
# present components are retained. 
#
# SysGO is extracted from a pre-processed file
# GO is extracted from the org.Hs.eg.db and the GO.db packages

# get sysgo preprocessed file
if (.Platform$OS.type == 'unix') {
  sysgo_database <- '/home/junkpp/work/other/2022-03-25_UpdateSysGO/sysgo_annotation.csv'
} else {
  sysgo_database <- 'D:/Daten/Work/Dublin/CloudStation/other/2022-03-25_UpdateSysGO/sysgo_annotation.csv'
}

# obtain GO BP annotation
annotation <- bind_rows(
  AnnotationDbi::keys(org.Hs.eg.db::org.Hs.egGO2ALLEGS) %>%
    map(function(go_id) {
      go_obj <- GO.db::GOTERM[[go_id]]
      term <- go_obj@Term
      ontology <- go_obj@Ontology
      definition <- go_obj@Definition
      if (ontology != 'BP') {
        return(NULL)
      } else {
        tibble(id = go_id,
               process = term,
               ontology = ontology,
               definition = definition)
      }
    }) %>%
    compact %>%
    bind_rows,
  read.csv(sysgo_database,header=T, sep=';')
) %>%
  mutate(definition = replace_na(definition, ''))
# str(annotation)



## build ontology data frame
ontology <- bind_rows(
  # SYSGO
  read.csv('D:/Daten/Work/Dublin/CloudStation/other/2022-03-25_UpdateSysGO/sysgo_filtered.csv',
           header=T, sep=';') %>% 
    select(-process) %>%
    group_by(id) %>%
    nest %>%
    mutate(data = map(data, function(df) {
      df %>%
        mutate(n_all = length(unique(hgnc))) %>%
        filter(hgnc %in% unique(df_apms$hgnc)) %>%
        mutate(n_found = length(unique(hgnc)))
    })) %>%
    unnest(cols = c(data)) %>%
    ungroup %>%
    select(id, hgnc, n_all, n_found) %>%
    mutate(ontology = 'SysGO'),
  # GO
  annotation %>% 
    filter(ontology == 'BP') %>% 
    pull(id) %>% 
    map(function(go_id) {
      # extract all entrez ids for a process
      proteins <- org.Hs.eg.db::org.Hs.egGO2ALLEGS[[go_id]] %>%
      unname %>%
      unique
      n_all <- length(proteins)
      # filter entrez ids by presence in APMS data
      proteins <- proteins[proteins %in% unique(df_apms$entrez_id)]
      # generate output
      if (length(proteins) == 0) {
        return(NULL)
      } else {
        return(
          tibble(
            id = go_id,
            entrez_id = proteins,
            n_all = n_all,
            n_found = length(proteins)
          )
        )
      }
    }) %>%
    compact %>%
    bind_rows %>%
    left_join(df_apms %>% select(hgnc, entrez_id) %>% distinct, 
              by = 'entrez_id') %>%
    select(id, hgnc, n_all, n_found) %>%
    mutate(ontology = 'BP')
  )

# str(ontology)


#####################
#####################
# Ancestor graphs
# 
# sysgo_graph <- bind_rows(
#   # process3 to process2
#   xlsx_cells(sysgo_database) %>%
#     filter(sheet == 'function & location') %>%
#     filter(row >= 2) %>%
#     select(character, row, col) %>%
#     group_by(row) %>%
#     summarize(symbol = character[col == 1],
#               process_1 = character[col == 2],
#               process_2 = character[col == 3],
#               process_3 = character[col == 4]) %>%
#     mutate(process_1 = case_when(process_1 == '1_26_29_adapters/scaffolds' ~ '1_26_29_adaptors/scaffolds',
#                                  T ~ process_1)) %>%
#     ungroup %>%
#     select(process_3, process_2) %>%
#     rename(from = process_3,
#            to = process_2)%>%
#     group_by(from, to) %>%
#     summarise(count = n()) %>%
#     ungroup,
#   xlsx_cells(sysgo_database) %>%
#     filter(sheet == 'function & location') %>%
#     filter(row >= 2) %>%
#     select(character, row, col) %>%
#     group_by(row) %>%
#     summarize(symbol = character[col == 1],
#               process_1 = character[col == 2],
#               process_2 = character[col == 3],
#               process_3 = character[col == 4]) %>%
#     mutate(process_1 = case_when(process_1 == '1_26_29_adapters/scaffolds' ~ '1_26_29_adaptors/scaffolds',
#                                  T ~ process_1)) %>%
#     ungroup %>%
#     select(-row) %>%
#     select(process_2, process_1) %>%
#     rename(from = process_2,
#            to = process_1) %>%
#     group_by(from, to) %>%
#     summarise(count = n()) %>%
#     ungroup
# )
# 
# 
# edges <- sysgo_graph
# nodes <- tibble(
#   id = unique(c(edges$to, edges$from))
# ) %>%
#   mutate(
#     # x = seq(from=1, to=n()),
#     # y = seq(from=1, to=n()),
#     label = id)
# 
# visNetwork(nodes, edges, width = '100%') %>%
#   # visIgraphLayout(
#   #   layout = 'layout_with_graphopt'
#   # ) %>%
#   visHierarchicalLayout() %>%
#   visPhysics(stabilization = FALSE) %>%
#   visEdges(
#     arrows = 'to',
#     smooth = F) %>%
#   visNodes(
#     shape = "square",
#     color = list(
#       background = "#0085AF",
#       border = "#013848",
#       highlight = "#FF8000"
#     ),
#     shadow = FALSE,
#     physics = FALSE
#   )



#######
####### continue ANOVA calculation
#######

###############################################################
# collapse processes on summed LFQ intensities per go process per sample

# TODO discuss filtering with CK
df_sum <- ontology %>%
  filter(n_all <= 1000) %>%
  filter(n_found >= 10) %>% 
  filter(n_found/n_all > 0.05) %>%
  inner_join(df_apms, by = 'hgnc') %>% 
  group_by(id, label, mut_status, condition, concentration) %>%
  summarise(sum_LFQ = sum(LFQ), 
            ontology = unique(ontology),
            .groups = 'keep') %>%
  mutate(sum_LFQ = log2(sum_LFQ)) %>%
  ungroup 

# str(df_sum)

###############################################################
# Perform ANOVA on processes

df_anova <- df_sum %>%
  pull(id) %>%
  unique %>%
  future_map(function(proc) {
    # create temporary dataframe
    df_proc <- df_sum %>%
      filter(id == proc) 
    # check levels in all factors:
    levels_mut_status <- df_proc$mut_status %>% unique %>% length
    levels_condition <- df_proc$condition %>% unique %>% length
    levels_concentration <- df_proc$concentration %>% unique %>% length    
    # construct formula based on number of levels in each factor
    formula_components <- character()
    if (levels_mut_status > 1) {
      formula_components <- c(formula_components, 'mut_status')
    }
    if (levels_condition > 1) {
      formula_components <- c(formula_components, 'condition') 
    }
    if (levels_concentration > 1) {
      formula_components <- c(formula_components, 'concentration')
    }
    # if no level was present, don't perform ANOVA
    if (length(formula_components) == 0) {
      return(NULL)
    }
    # construct formula
    formula_anova <-  str_c(
      c('sum_LFQ', str_c(formula_components, collapse='*')), 
      collapse='~') %>%
      as.formula()
    
    # perform ANOVA
    anova <- df_proc %>%
      aov(data=., formula_anova)
    # return results from ANOVA
    anova %>% summary %>% pluck(1) %>%
      mutate(Factor = anova %>% summary %>% pluck(1) %>% rownames,
             id = proc) %>%
      rename(p_val = 'Pr(>F)')
  }) %>%
  compact %>%
  bind_rows %>%
  mutate(p_adj = p.adjust(p_val, 'holm'),
         Factor = str_trim(Factor))

str(df_anova)


df_tukey <- df_anova %>%
  filter(p_adj < 0.05) %>%
  pull(id) %>%
  unique %>% 
  future_map(function(proc) {
    # create temporary dataframe
    df_proc <- df_sum %>%
      filter(id == proc) 
    # check levels in all factors:
    levels_mut_status <- df_proc$mut_status %>% unique %>% length
    levels_condition <- df_proc$condition %>% unique %>% length
    levels_concentration <- df_proc$concentration %>% unique %>% length    
    # construct formula based on number of levels in each factor
    formula_components <- character()
    if (levels_mut_status > 1) {
      formula_components <- c(formula_components, 'mut_status')
    }
    if (levels_condition > 1) {
      formula_components <- c(formula_components, 'condition') 
    }
    if (levels_concentration > 1) {
      formula_components <- c(formula_components, 'concentration')
    }
    # if no level was present, don't perform ANOVA
    if (length(formula_components) == 0) {
      return(NULL)
    }
    # construct formula
    formula_anova <-  str_c(
      c('sum_LFQ', str_c(formula_components, collapse='*')), 
      collapse='~') %>%
      as.formula()
    # extract signficiant terms
    which_terms <- df_anova %>%
      filter(id == proc) %>%
      filter(p_adj < 0.05) %>%
      pull(Factor)
    
    # perform tukey post_hoc test
    tukey <- df_proc %>%
      aov(data=., formula_anova) %>%
      TukeyHSD(which = which_terms)
    
    # collect tukey post hoc output
    out <- names(tukey) %>%
      map(function(factr) {
        tukey[[factr]] %>%
          as.data.frame() %>%
          mutate(term = factr, 
                 Comparison = rownames(tukey[[factr]]))
      }) %>%
      bind_rows %>% 
      rename(p_adj = 'p adj') %>% 
      mutate(id = proc) %>%
      filter(!is.na(diff))
      
  }) %>%
  compact %>%
  bind_rows %>%
  mutate(p_adj = p.adjust(p_adj, method = 'fdr'),
         group1 = str_extract(Comparison, '(?<=-)[:graph:]+'),
         group2 = str_extract(Comparison, '[:graph:]+(?=-)')) %>%
  select(id, term, group1, group2, p_adj) %>%
  filter(p_adj <= 0.1)


str(df_tukey)

# # old approach
# df_anova <- df_sum %>% 
#   pull(id) %>%
#   unique %>%
#   future_map(function(proc) {
#     # create temporary dataframe
#     df_proc <- df_sum %>%
#       filter(id == proc)
#     # check levels in all factors:
#     levels_mut_status <- df_proc$mut_status %>% unique %>% length
#     levels_condition <- df_proc$condition %>% unique %>% length
#     levels_concentration <- df_proc$concentration %>% unique %>% length    
#     # construct formula based on number of levels in each factor
#     formula_components <- character()
#     if (levels_mut_status > 1) {
#       formula_components <- c(formula_components, 'mut_status')
#     }
#     if (levels_condition > 1) {
#       formula_components <- c(formula_components, 'condition') 
#     }
#     if (levels_concentration > 1) {
#       formula_components <- c(formula_components, 'concentration')
#     }
#     # if no level was present, don't perform ANOVA
#     if (length(formula_components) == 0) {
#       return(NULL)
#     }
#     # construct formula
#     formula_anova <-  str_c(
#       c('sum_LFQ', str_c(formula_components, collapse='*')), 
#       collapse='~') %>%
#       as.formula()
#     
#     # perform ANOVA
#     anova <- df_proc %>%
#       aov(data=., formula_anova)
#     # get significant factors from ANOVA
#     Factors <- anova %>% summary %>% pluck(1) %>%
#       mutate(Factor = anova %>% summary %>% pluck(1) %>% rownames) %>% 
#       rename(p_val = 'Pr(>F)') %>% 
#       filter(p_val < 0.05) %>% 
#       pull(Factor) %>%
#       str_trim
#     if (length(Factors) == 0) {
#       return(NULL)
#     }
#     # extract residual standard deviation from anova object
#     # links
#     # https://www.graphpad.com/support/faqid/1564/
#     # http://www.balkinresearchmethods.com/Balkin_Research_Methods/Research_Methods_and_Statistics_files/Effect%20Size%20slides.pdf
#     # https://stats.stackexchange.com/questions/422151/calculate-cohens-d-for-pairwise-tukey-tests-in-a-1-way-anova
#     # https://www.reddit.com/r/statistics/comments/hfj65x/q_how_to_calculate_cohen_d_effect_size_on_tukey/
#     # see also BruceR's EMMEANS code
#     # https://github.com/psychbruce/bruceR/blob/82be838bcc9aa19ba21916a3f59e0f34a180fd2d/R/bruceR_stats_03_manova.R#L963
#     pooled_sd <- sigma(anova)
#     # perform Tukey post hoc test on significant factors
#     tukey <- TukeyHSD(anova, which = Factors)
#     # collect post hoc output
#     out <- names(tukey) %>%
#       map(function(factr) {
#         tukey[[factr]] %>%
#           as.data.frame() %>%
#           mutate(Factor = factr, 
#                  Comparison = rownames(tukey[[factr]]))
#       }) %>%
#       bind_rows %>% 
#       rename(p_adj = 'p adj') %>% 
#       mutate(id = proc) %>%
#       filter(!is.na(diff)) %>%
#       mutate(cohens_d = round(abs(diff)/pooled_sd, 3)) %>%
#       select(id, Factor, Comparison, p_adj, cohens_d)
#     return(out)
#   }) %>%
#   compact %>%
#   bind_rows() %>%
#   mutate(p_adj = p.adjust(p_adj, method = 'bonferroni')) 
# 
# # str(df_anova)

#######
####### GSEA
#######

# TODO move calculation for differential expression and GSEA here as well!
df_gsea <- bind_rows(
  read.csv('D:/Daten/Work/Dublin/CloudStation/12_CT_APMS/analysis/gse_gobp.csv',
           header=T, sep=';') %>%
    mutate(ontology = 'BP'),
  read.csv('D:/Daten/Work/Dublin/CloudStation/12_CT_APMS/analysis/gse_sysgo.csv',
           header=T, sep=';') %>%
    mutate(ontology = 'SysGO')
) %>% 
  select(ID, NES, pvalue, contrast, ontology) %>%
  rename(id = ID) %>%
  # adjust pvalue
  mutate(p_adj = p.adjust(pvalue, method = 'fdr')) %>%
  filter(p_adj <= 0.1) %>% 
  # retrieve which condition is enriched compared to which other condition
  mutate(enriched_a = str_extract(contrast, '[:graph:]+(?=_vs)'),
         enriched_b = str_extract(contrast, '(?<=vs_)[:graph:]+')) %>%
  mutate(enriched_a = str_replace_all(enriched_a, 'X', 'g'),
         enriched_b = str_replace_all(enriched_b, 'X', 'g')) %>%
  mutate(enriched_in = case_when(NES > 0 ~ enriched_a,
                                 NES < 0 ~ enriched_b),
         enriched_against = case_when(NES > 0 ~ enriched_b,
                                      NES < 0 ~ enriched_a),
         NES = abs(NES)) %>% 
  select(-enriched_a, -enriched_b, -contrast)
  

# TODO maybe look at core enriched genes





#######
####### Outputs

app_data_dir <- 'D:/Daten/Work/Dublin/CloudStation/12_CT_APMS/kras_apms/data/'

annotation %>%
  write_csv(str_c(app_data_dir, 'df_annotation.csv'))

ontology %>%
  write_csv(str_c(app_data_dir, 'df_ontology.csv'))

df_apms %>%
  write_csv(str_c(app_data_dir, 'df_apms.csv'))

df_sum %>%
  write_csv(str_c(app_data_dir, 'df_sum.csv'))

df_tukey %>%
  # pre-filtering to reduce file size 
  filter(p_adj <= 0.1) %>%
  write_csv(str_c(app_data_dir, 'df_anova.csv'))
  

df_gsea %>%
  write_csv(str_c(app_data_dir, 'df_gsea.csv'))


## Supplimentary table for Christina
df_tukey %>% 
  filter(p_adj < 0.05) %>% 
  count(id, term) %>% 
  inner_join(df_annotation %>% select(id, process, ontology), by='id') %>% 
  pivot_wider(names_from=term, values_from=n, values_fill=0) %>% 
  select(ontology, id, process, condition, mut_status, concentration, everything()) %>% 
  arrange(ontology, -(condition+mut_status+concentration)) %>% 
  write_delim('D:/Daten/Work/Dublin/CloudStation/12_CT_APMS/analysis/supp_table_anova.csv', delim=';')

# #####################################################
# #### TEST PCA of GO terms
# # Do I retain information?
# 
# df_wide <- df_sum %>%
#   filter(condition %in% c('unstim', 'dmog')) %>%
#   filter(mut_status == '12c') %>% 
#   filter(concentration != '20') %>% 
#   filter(ontology == 'BP') %>%
#   select(-ontology) %>%
#   pivot_wider(names_from = id, values_from = sum_LFQ, values_fill = 0)
# 
# dim(df_wide)
# 
# pca <- df_wide %>%
#   select(where(is.numeric)) %>%
#   select(where(~ sd(.x) > 0)) %>%
#   prcomp(., scale=TRUE, center=TRUE)
# 
# # variance explained for first 2 components
# pc1_varexpl <- round(summary(pca)$importance[2,1] * 100, 2)
# pc2_varexpl <- round(summary(pca)$importance[2,2] * 100, 2)
# 
# # full pca
# p_pca_go <- pca %>% broom::augment(df_wide) %>%
#   ggplot(aes(x = .fittedPC1, y = .fittedPC2,
#              color=condition, shape = concentration)) +
#   #geom_density_2d() +
#   geom_point(size=4) +
#   labs(x = str_c('PC1 (', pc1_varexpl, ' %)'),
#        y = str_c('PC2 (', pc2_varexpl, ' %)')
#   ) +
#   facet_wrap(vars(mut_status)) +
#   theme_bw()
# 
# 
# p_pca_go
# 
# # loading
# loadings <- pca$rotation %>% 
#   as.data.frame %>% 
#   select(PC1, PC2) %>%
#   mutate(id = rownames(pca$rotation)) 
# 
# str(loadings)
# 
# loadings %>%
#   arrange(desc(abs(PC1))) %>% head
# 
# loadings %>%
#   arrange(desc(abs(PC2))) %>% head
# 
# ggplot(loadings, aes(x = PC1, y = PC2)) + geom_point()
# loadings seem to be pretty shit
# correlations of only around 0.06


# This looks really cool, good conservation of variance
# NICE!