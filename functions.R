library(ComplexHeatmap)
library(tidyverse)

get_go_term = function(go_id) {
  term = suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = go_id, columns = "TERM")$TERM)
  term[is.na(term)] = "NA"
  term
}

# click_action = function(df, output) {
#   output[["go_info"]] = renderUI({
#     if(!is.null(df)) {
#       go_id1 = rownames(full_dist_mat)[df$row_index]
#       go_id2 = colnames(full_dist_mat)[df$column_index]
#       
#       HTML(str_glue(
#         "<pre>
# ## Row GO ID
# <a href='http://amigo.geneontology.org/amigo/term/{go_id1}' target='_blank'>{go_id1}</a>: {get_go_term(go_id1)}
# 
# ## Column GO ID:
# <a href='http://amigo.geneontology.org/amigo/term/{go_id2}' target='_blank'>{go_id2}</a>: {get_go_term(go_id2)}
# </pre>"
#       ))
#     }
#   })
# }
# 
# brush_action = function(df, output) {
#   output[["go_info"]] = renderUI({
#     if(!is.null(df)) {
#       row_index = unique(unlist(df$row_index))
#       column_index = unique(unlist(df$column_index))
#       go_id1 = rownames(full_dist_mat)[row_index]
#       go_id2 = colnames(full_dist_mat)[column_index]
#       
#       go_id = union(go_id1, go_id2)
#       
#       go_text = str_glue("<a href='http://amigo.geneontology.org/amigo/term/{go_id}' target='_blank'>{go_id}</a>: {get_go_term(go_id)} <button id='{go_id}' class='go_sel_button'>Select</button>") %>% 
#         str_c(collapse='\n')
#       HTML(str_glue("<pre>{go_text}</pre>"))
#     }
#   })
# }

###############################################################################
# make custom heatmap

extract_gsea_condition <- function(df_gsea, reference = 'unstim_none', all_terms) {
  df <- bind_rows(
    df_gsea %>% 
      filter(str_detect(enriched_against, str_glue('{reference}$'))) %>%
      filter(str_extract(enriched_in, '[:alnum:]+(?=_)') == str_extract(enriched_against, '[:alnum:]+(?=_)')) %>%
      mutate(group = enriched_in,
             score = NES * -log10(p_adj)),
    df_gsea %>% filter(str_detect(enriched_in, str_glue('{reference}$'))) %>%
      filter(str_extract(enriched_in, '[:alnum:]+(?=_)') == str_extract(enriched_against, '[:alnum:]+(?=_)')) %>%
      mutate(group = enriched_against,
             score = -(sign(NES) * -log10(p_adj)))
  ) %>% 
    filter(ontology == 'BP') %>%
    select(id, group, score) %>% 
    full_join(tibble(id = all_terms), by='id') %>% 
    pivot_wider(names_from = group, values_from = score, values_fill = 0) %>% 
    select(-'NA') %>% 
    filter(id %in% all_terms)
  mat <- df %>% select(where(is.numeric)) %>% as.matrix
  rownames(mat) <- df$id
  mat <- mat[all_terms,,drop=F]
  mat
}

extract_gsea_mutstatus <- function(df_gsea, reference = 'wt', all_terms) {
  df <- bind_rows(
    df_gsea %>% 
      filter(str_detect(enriched_against, str_glue('^{reference}(?=_)'))) %>%
      filter(str_extract(enriched_in, '(?<=_)[:alnum:]+_[:alnum:]+') == str_extract(enriched_against, '(?<=_)[:alnum:]+_[:alnum:]+')) %>%
      mutate(group = enriched_in,
             score = NES * -log10(p_adj)),
    df_gsea %>% filter(str_detect(enriched_in, str_glue('^{reference}(?=_)'))) %>%
      filter(str_extract(enriched_in, '(?<=_)[:alnum:]+_[:alnum:]+') == str_extract(enriched_against, '(?<=_)[:alnum:]+_[:alnum:]+')) %>%
      mutate(group = enriched_against,
             score = -(sign(NES) * -log10(p_adj)))
  ) %>% 
    filter(ontology == 'BP') %>%
    select(id, group, score) %>% 
    full_join(tibble(id = all_terms), by='id') %>% 
    pivot_wider(names_from = group, values_from = score, values_fill = 0) %>% 
    select(-'NA') %>% 
    filter(id %in% all_terms)
  mat <- df %>% select(where(is.numeric)) %>% as.matrix
  rownames(mat) <- df$id
  mat <- mat[all_terms,,drop=F]
  mat
}

extract_anova_condition <- function(df_anova, reference = 'unstim', all_terms) {
  ## TODO calculate score based on difference in mean and pvalue in the future to indicate direction!
  df <- df_anova %>%
    filter(str_detect(id, '^GO')) %>% # TODO not necessary if only GO 
    filter(term == 'condition') %>%  
    filter(group1 == reference | group2 == reference) %>% 
    mutate(group = case_when(group1 == reference ~ group2,
                             group2 == reference ~ group1)) %>% 
    mutate(score = -log10(p_adj)) %>%
    select(id, group, score) %>%
    full_join(tibble(id = all_terms), by='id') %>% 
    pivot_wider(names_from = group, values_from = score, values_fill = 0) %>% 
    select(-'NA') %>% 
    filter(id %in% all_terms)
  mat <- df %>% select(where(is.numeric)) %>% as.matrix
  rownames(mat) <- df$id
  mat <- mat[all_terms,,drop=F]
  # t(mat)
}

extract_anova_mutstatus <- function(df_anova, reference = 'wt', all_terms) {
  ## TODO calculate score based on difference in mean and pvalue in the future to indicate direction!
  df <- df_anova %>%
    filter(str_detect(id, '^GO')) %>% # TODO not necessary if only GO 
    filter(term == 'mut_status') %>%  
    filter(group1 == reference | group2 == reference) %>% 
    mutate(group = case_when(group1 == reference ~ group2,
                             group2 == reference ~ group1)) %>% 
    mutate(score = -log10(p_adj)) %>%
    select(id, group, score) %>%
    full_join(tibble(id = all_terms), by='id') %>% 
    pivot_wider(names_from = group, values_from = score, values_fill = 0) %>% 
    select(-'NA') %>% 
    filter(id %in% all_terms)
  mat <- df %>% select(where(is.numeric)) %>% as.matrix
  rownames(mat) <- df$id
  mat <- mat[all_terms,,drop=F]
  # t(mat)
}

# TODO transform data matrices into discrete matrix
# TODO annotation 
# TODO integrate into shiny
# TODO 
# TODO rearrange columns for GSEA

custom_ht_clusters = function(
    dist_mat, 
    cl,
    df_gsea,
    df_anova,
    ref_gsea_condition = NULL,
    ref_gsea_mut_status = NULL,
    ref_anova_condition = NULL,
    ref_anova_mut_status = NULL,
    reduce_matrix = F)
{
  
  # settings
  col = c('white', 'black')
  min_term = 10
  
  # init
  mat_list <- list()
  
  # get data if references are specified
  if (!is.null(ref_gsea_condition)) {
    mat_list$gsea_condition <- extract_gsea_condition(
      df_gsea %>% filter(p_adj<0.05),
      ref_gsea_condition,
      all_terms = rownames(dist_mat))
  }
  if (!is.null(ref_gsea_mut_status)) {
    mat_list$gsea_mut_status <- extract_gsea_mutstatus(
      df_gsea %>% filter(p_adj<0.05), 
      ref_gsea_mut_status,
      all_terms = rownames(dist_mat))
  }
  if (!is.null(ref_anova_condition)) {
    mat_list$anova_condition <- extract_anova_condition(
      df_anova %>% filter(p_adj<0.05), 
      ref_anova_condition,
      all_terms = rownames(dist_mat))
  }
  if (!is.null(ref_anova_mut_status)) {
    mat_list$anova_mut_status <- extract_anova_mutstatus(
      df_anova %>% filter(p_adj<0.05), 
      ref_anova_mut_status,
      all_terms = rownames(dist_mat))
  }
  
  # reduce matrix by showing only observations found in all sub-heatmaps from observations
  if (reduce_matrix) {
    if (length(mat_list) > 0) {
      # retrieve terms to keep from matrices
      show_terms <- mat_list %>% 
        map(function(m) {m %>% rowSums %>% `==`(0)}) %>% 
        transpose %>% 
        discard(function(l) {flatten_lgl(l) %>% all}) %>% 
        names
    } else {
      show_terms <- rownames(dist_mat)
    }
    # filter distance matrix and clustering
    cl <- cl[rownames(dist_mat) %in% show_terms]
    dist_mat <- dist_mat[rownames(dist_mat) %in% show_terms,
                         colnames(dist_mat) %in% show_terms]
    # filter score matrices
    for (i in 1:length(mat_list)) {
      mat <- mat_list[[i]]
      mat_list[[i]] <- mat[rownames(mat) %in% show_terms,]
    }
  }
  
  # color function for similarity matrix
  col_fun = circlize::colorRamp2(
    seq(0, quantile(dist_mat[dist_mat > 0], 0.975), length = length(col)), col)
  
  # generate order of GO terms
  cl = as.vector(cl)
  cl_tb = table(cl)
  cl[as.character(cl) %in% names(cl_tb[cl_tb < min_term])] = 0
  cl = factor(cl, levels = c(setdiff(sort(cl), 0), 0))
  cl = factor(cl, levels = c(setdiff(names(sort(table(cl), decreasing = TRUE)), 0), 0))
  od2 = unlist(lapply(levels(cl), function(le) {
    l = cl == le
    if(sum(l) <= 1) {
      return(which(l))
    } else {
      mm = dist_mat[l, l, drop = FALSE]
      which(l)[hclust(stats::dist(mm))$order]
    }
  }))
  
  # create and annotate similarity heatmap
  ht = Heatmap(
    dist_mat, 
    col = col_fun,
    name = "Similarity", 
    column_title = NULL,
    show_row_names = FALSE, 
    show_column_names = FALSE,
    show_row_dend = FALSE, 
    show_column_dend = FALSE,
    row_order = od2, 
    column_order = od2,
    border = "#404040", 
    row_title = NULL,
    use_raster = T, 
    width = 10, height = 10) + 
    NULL
  
  ht@ht_list[[1]]@heatmap_param$post_fun = function(ht) {
    decorate_heatmap_body("Similarity", {
      grid.rect(gp = gpar(fill = NA, col = "#404040"))
      cl = factor(cl, levels = unique(cl[od2]))
      tbcl = table(cl)
      ncl = length(cl)
      x = cumsum(c(0, tbcl))/ncl
      grid.segments(x, 0, x, 1, default.units = "npc", gp = gpar(col = "#404040"))
      grid.segments(0, 1 - x, 1, 1 - x, default.units = "npc", gp = gpar(col = "#404040"))
    })
  }
  
  # create data heatmaps
  make_data_heatmap <- function(mat, name, width) {
    Heatmap(
      mat, 
      col = circlize::colorRamp2(c(-1, 0, 1), c(scales::muted('blue'), 'white', scales::muted('red'))),
      name = name, 
      column_title = NULL,
      show_row_names = FALSE, 
      show_column_names = TRUE,
      show_row_dend = FALSE, 
      show_column_dend = FALSE,
      row_order = od2, 
      cluster_rows=FALSE, 
      cluster_columns = FALSE,
      border = "#404040", 
      row_title = NULL,
      use_raster = T, 
      width = width) + 
      NULL
  }
  
  if (!is.null(ref_gsea_condition)) {
    ht_gsea_condition <- make_data_heatmap(mat_list$gsea_condition, 'GSEA_condition', 5)
  } else {
    ht_gsea_condition <- NULL
  }
  
  if (!is.null(ref_gsea_mut_status)) {
    ht_gsea_mut_status <- make_data_heatmap(mat_list$gsea_mut_status, 'GSEA_mut_status', 5)
  } else {
    ht_gsea_mut_status <- NULL
  }
  
  if (!is.null(ref_anova_condition)) {
    ht_anova_condition <- make_data_heatmap(mat_list$anova_condition, 'ANOVA_condition', 1)
  } else {
    ht_anova_condition <- NULL 
  }
  
  if (!is.null(ref_anova_mut_status)) {
    ht_anova_mut_status <- make_data_heatmap(mat_list$anova_mut_status, 'ANOVA_mut_status', 1)
  } else {
    ht_anova_mut_status <- NULL
  }
  
  ht_list <- (ht_anova_mut_status + (ht_anova_condition + (ht_gsea_mut_status + (ht_gsea_condition + ht))))
  
  return(invisible(ht_list))
}
