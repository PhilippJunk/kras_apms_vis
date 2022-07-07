library(ComplexHeatmap)
library(tidyverse)

get_go_term <- function(go_id, df_annotation) {
  left_join(tibble(id = go_id), df_annotation, by = 'id') %>% 
    mutate(process = replace_na(process, 'NA')) %>%
    pull(process)
}

###############################################################################
# make custom heatmap

extract_gsea_condition <- function(df_gsea, reference = 'unstim_none', all_terms) {
  df <- bind_rows(
    df_gsea %>% 
      filter(str_detect(enriched_against, str_glue('{reference}$'))) %>%
      filter(str_extract(enriched_in, '[:alnum:]+(?=_)') == str_extract(enriched_against, '[:alnum:]+(?=_)')) %>%
      mutate(group = enriched_in,
             score = 1),
    df_gsea %>% filter(str_detect(enriched_in, str_glue('{reference}$'))) %>%
      filter(str_extract(enriched_in, '[:alnum:]+(?=_)') == str_extract(enriched_against, '[:alnum:]+(?=_)')) %>%
      mutate(group = enriched_against,
             score = -1)
  ) %>% 
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
             score = 1),
    df_gsea %>% filter(str_detect(enriched_in, str_glue('^{reference}(?=_)'))) %>%
      filter(str_extract(enriched_in, '(?<=_)[:alnum:]+_[:alnum:]+') == str_extract(enriched_against, '(?<=_)[:alnum:]+_[:alnum:]+')) %>%
      mutate(group = enriched_against,
             score = -1)
  ) %>% 
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
  df <- df_anova %>% 
    filter(term == 'condition') %>%
    filter(group_higher == reference | group_lower == reference) %>%
    mutate(group = case_when(group_higher == reference ~ group_lower,
                             group_lower == reference ~ group_higher),
           score = case_when(group_higher == reference ~ -1,
                             group_lower == reference ~ 1)) %>% 
    select(id, group, score) %>%
    full_join(tibble(id = all_terms), by='id') %>% 
    pivot_wider(names_from = group, values_from = score, values_fill = 0) %>% 
    select(-'NA') %>% 
    filter(id %in% all_terms)
  mat <- df %>% select(where(is.numeric)) %>% as.matrix
  rownames(mat) <- df$id
  mat <- mat[all_terms,,drop=F]
}

extract_anova_mutstatus <- function(df_anova, reference = 'wt', all_terms) {
  df <- df_anova %>% 
    filter(term == 'mut_status') %>%
    filter(group_higher == reference | group_lower == reference) %>%
    mutate(group = case_when(group_higher == reference ~ group_lower,
                             group_lower == reference ~ group_higher),
           score = case_when(group_higher == reference ~ -1,
                             group_lower == reference ~ 1)) %>% 
    select(id, group, score) %>%
    full_join(tibble(id = all_terms), by='id') %>% 
    pivot_wider(names_from = group, values_from = score, values_fill = 0) %>% 
    select(-'NA') %>% 
    filter(id %in% all_terms)
  mat <- df %>% select(where(is.numeric)) %>% as.matrix
  rownames(mat) <- df$id
  mat <- mat[all_terms,,drop=F]
}

# TODO annotation 
# TODO rearrange columns for GSEA
# TODO adapt widths
# TODO figure out a way to show cluster wordclouds. Somehow. 
# (thinking about this one, it could be a pre_arranged grid object on the right side of the heatmap,
# with the individual heatmaps height-stacked according to the number of terms / cluster)

custom_ht_clusters = function(
    dist_mat, 
    cl,
    df_gsea,
    df_anova,
    ref_gsea_condition = NULL,
    ref_gsea_mut_status = NULL,
    ref_anova_condition = NULL,
    ref_anova_mut_status = NULL,
    reduce_matrix = F,
    draw_word_cloud = F) {
  
  # settings
  col = c('white', 'black')
  min_term = 20
  
  # stat = 'count'
  stat = 'pvalue'
  min_stat = ifelse(stat == "count", 5, 0.05)
  exclude_words = NULL
  max_words = 10
  word_cloud_grob_param = list()
  fontsize_range = c(4, 16)
  bg_gp = gpar(fill = "#FFFFFF", col = "#AAAAAA")
  
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
  if (reduce_matrix & length(mat_list) > 0) {
    # retrieve terms to keep from matrices
    show_terms <- mat_list %>% 
      map(function(m) {m %>% abs %>% rowSums %>% `==`(0)}) %>% 
      transpose %>% 
      discard(function(l) {flatten_lgl(l) %>% all}) %>% 
      names
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
    heatmap_legend_param = list(title = 'GO Semantic\nSimilarity'),
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
  
  # word cloud
  if(draw_word_cloud) {
    go_id = rownames(dist_mat)
    
    align_to = split(seq_along(cl), cl)
    go_id = split(go_id, cl)
    
    align_to = align_to[names(align_to) != "0"]
    go_id = go_id[names(go_id) != "0"]
    
    if(length(align_to)) {
      ht = ht + rowAnnotation(keywords = anno_word_cloud_from_GO(align_to, go_id = go_id, stat = stat, min_stat = min_stat,
                                                                 exclude_words = exclude_words, max_words = max_words, word_cloud_grob_param = word_cloud_grob_param, 
                                                                 fontsize_range = fontsize_range, bg_gp = bg_gp))
    } 
  }
  
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
  make_data_heatmap <- function(mat, name, 
                                width, 
                                plot_main_legend=FALSE, 
                                top_anno = NULL, 
                                bottom_anno = NULL) 
    {
    Heatmap(
      mat, 
      col = c('-1' = scales::muted('blue'), '0' = 'white', '1' = scales::muted('red')),
      show_heatmap_legend = plot_main_legend,
      heatmap_legend_param = list(
        at = c(-1, 0, 1), title = 'Comparison to\nreference', 
        legend_gp = gpar(fill = c(c(scales::muted('blue'), 'white', scales::muted('red')))),
        labels = c('sig. lower', 'no sig. diff', 'sig. higher')),
      top_annotation = top_anno,
      bottom_annotation = bottom_anno,
      name = name, 
      column_title = NULL,
      show_row_names = FALSE, 
      show_column_names = FALSE,
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
  
  # create legends
  if (any(!is.na(c(ref_gsea_mut_status, ref_gsea_condition, 
                   ref_anova_mut_status, ref_anova_condition)))) {
    df_col_mut_status <- tibble(mut_status = c('wt', 'g12c', 'g12d', 'g12v'),
                                labels = c('WT', 'G12C', 'G12D', 'G12V'),
                                colors = c('#A98743', '#437C90', '#255957', '#8F2844'))
    df_col_condition <- tibble(condition = c('unstim', 'dmog', 'egf', 'il6', 'pge2', 'tnfa'),
                               labels = c('Unstim', 'DMOG', 'EGF', 'IL6', 'PGE2', 'TNFA'),
                               colors = c("#001524","#12616D","#75964A","#A1869E","#FF7D00","#78290F"))
    df_col_concentration <- tibble(concentration = c('none', '20', '200'),
                                   labels = c('None', '20 ng/ml', '200 ng/ml'),
                                   colors = scales::brewer_pal(palette=11)(3))
    
    annotation_lgd_list <- list(
      'Mutation Status' = list(
        title = 'Mutation Status',
        at = df_col_mut_status$mut_status, 
        labels = df_col_mut_status$labels,
        legend_gp = gpar(fill = df_col_mut_status$colors)),
      'Condition' = list(
        title = 'Condition',
        at = df_col_condition$condition,
        labels = df_col_condition$labels,
        legend_gp = gpar(fill = df_col_condition$colors)
      ),
      'Concentration' = list(
        title = 'Concentration',
        at = df_col_concentration$concentration,
        labels = df_col_concentration$labels,
        legend_gp = gpar(fill = df_col_concentration$colors)
      ))
    plot_anno_legend_mut <- T
    plot_anno_legend_cond <- T
    plot_anno_legend_conc <- T
  }
  
  # create data heatmaps
  plot_main_legend <- TRUE
  if (!is.null(ref_gsea_condition)) {
    # create annotations
    mut_status_top <- str_extract(colnames(mat_list$gsea_condition), '^[:alnum:]+(?=_)')
    condition_top <- str_extract(colnames(mat_list$gsea_condition), '(?<=_)[:alnum:]+(?=_)')
    concentration_top <- str_extract(colnames(mat_list$gsea_condition), '(?<=_)[:alnum:]+$')
    top_anno <- HeatmapAnnotation(
      'Mutation Status' = mut_status_top, 
      'Condition' = condition_top, 
      'Concentration' = concentration_top,
      col = list('Mutation Status' = structure(df_col_mut_status$colors, names=df_col_mut_status$mut_status),
                 'Condition' = structure(df_col_condition$colors, names=df_col_condition$condition),
                 'Concentration' = structure(df_col_concentration$colors, names=df_col_concentration$concentration)),
      annotation_legend_param = annotation_lgd_list[c(plot_anno_legend_mut, plot_anno_legend_cond, 
                                                      plot_anno_legend_conc)],
      show_legend = any(plot_anno_legend_mut, plot_anno_legend_cond, plot_anno_legend_conc),
      border = T, show_annotation_name = F)
    
    
    mut_status_bottom <- mut_status_top
    condition_bottom <- rep(str_extract(ref_gsea_condition,'^[:alnum:]+(?=_)'), 
                            length(colnames(mat_list$gsea_condition)))
    concentration_bottom <- rep(str_extract(ref_gsea_condition, '(?<=_)[:alnum:]+$'), 
                                length(colnames(mat_list$gsea_condition)))
    bottom_anno <- HeatmapAnnotation(
      'Mutation Status' = mut_status_bottom, 
      'Condition' = condition_bottom, 
      'Concentration' = concentration_bottom,
      col = list('Mutation Status' = structure(df_col_mut_status$colors, names=df_col_mut_status$mut_status),
                 'Condition' = structure(df_col_condition$colors, names=df_col_condition$condition),
                 'Concentration' = structure(df_col_concentration$colors, names=df_col_concentration$concentration)),
      show_legend = F, 
      border = T, show_annotation_name = F)
    
    # create heatmap
    ht_gsea_condition <- make_data_heatmap(mat_list$gsea_condition, 
                                           'GSEA_condition', 5, 
                                           plot_main_legend,
                                           top_anno, bottom_anno)
    # update legend settings
    plot_main_legend <- FALSE
    plot_anno_legend_mut <- FALSE
    plot_anno_legend_cond <- FALSE
    plot_anno_legend_conc <- FALSE
  } else {
    ht_gsea_condition <- NULL
  }
  
  if (!is.null(ref_gsea_mut_status)) {
    # create annotations
    mut_status_top <- str_extract(colnames(mat_list$gsea_mut_status), '^[:alnum:]+(?=_)')
    condition_top <- str_extract(colnames(mat_list$gsea_mut_status), '(?<=_)[:alnum:]+(?=_)')
    concentration_top <- str_extract(colnames(mat_list$gsea_mut_status), '(?<=_)[:alnum:]+$')
    top_anno <- HeatmapAnnotation(
      'Mutation Status' = mut_status_top, 
      'Condition' = condition_top, 
      'Concentration' = concentration_top,
      col = list('Mutation Status' = structure(df_col_mut_status$colors, names=df_col_mut_status$mut_status),
                 'Condition' = structure(df_col_condition$colors, names=df_col_condition$condition),
                 'Concentration' = structure(df_col_concentration$colors, names=df_col_concentration$concentration)),
      annotation_legend_param = annotation_lgd_list[c(plot_anno_legend_mut, plot_anno_legend_cond, 
                                                      plot_anno_legend_conc)],
      show_legend = any(plot_anno_legend_mut, plot_anno_legend_cond, plot_anno_legend_conc),
      border = T, show_annotation_name = F)
    
    
    mut_status_bottom <-rep(ref_gsea_mut_status, length(colnames(mat_list$gsea_mut_status)))
    condition_bottom <- condition_top
    concentration_bottom <- concentration_top
    bottom_anno <- HeatmapAnnotation(
      'Mutation Status' = mut_status_bottom, 
      'Condition' = condition_bottom, 
      'Concentration' = concentration_bottom,
      col = list('Mutation Status' = structure(df_col_mut_status$colors, names=df_col_mut_status$mut_status),
                 'Condition' = structure(df_col_condition$colors, names=df_col_condition$condition),
                 'Concentration' = structure(df_col_concentration$colors, names=df_col_concentration$concentration)),
      show_legend = FALSE,
      border = T, show_annotation_name = F)
    
    ht_gsea_mut_status <- make_data_heatmap(mat_list$gsea_mut_status, 
                                            'GSEA_mut_status', 5, 
                                            plot_main_legend,
                                            top_anno, bottom_anno)
    plot_main_legend <- FALSE
    plot_anno_legend_mut <- FALSE
    plot_anno_legend_cond <- FALSE
    plot_anno_legend_conc <- FALSE
  } else {
    ht_gsea_mut_status <- NULL
  }
  
  if (!is.null(ref_anova_condition)) {
    # create annotations
    condition_top <- colnames(mat_list$anova_condition)
    top_anno <- HeatmapAnnotation(
      'Condition' = condition_top, 
      col = list('Condition' = structure(df_col_condition$colors, names=df_col_condition$condition)),
      annotation_legend_param = annotation_lgd_list[c(F, plot_anno_legend_cond, F)],
      show_legend = plot_anno_legend_cond,
      border = T, show_annotation_name = F)
    
    condition_bottom <- rep(ref_anova_condition, length(condition_top))
    bottom_anno <- HeatmapAnnotation(
      'Condition' = condition_bottom, 
      col = list('Condition' = structure(df_col_condition$colors, names=df_col_condition$condition)),
      annotation_legend_param = annotation_lgd_list[c(F, plot_anno_legend_cond, F)],
      show_legend = FALSE,
      border = T, show_annotation_name = F)
    
    
    ht_anova_condition <- make_data_heatmap(mat_list$anova_condition, 
                                            'ANOVA_condition', 1,
                                            plot_main_legend,
                                            top_anno, bottom_anno)
    plot_main_legend <- FALSE
    plot_anno_legend_cond <- FALSE
  } else {
    ht_anova_condition <- NULL 
  }
  
  if (!is.null(ref_anova_mut_status)) {
    # create annotations
    mut_status_top <- colnames(mat_list$anova_mut_status)
    top_anno <- HeatmapAnnotation(
      'Mutation Status' = mut_status_top, 
      col = list('Mutation Status' = structure(df_col_mut_status$colors, names=df_col_mut_status$mut_status)),
      annotation_legend_param = annotation_lgd_list[c(plot_anno_legend_mut, F, F)],
      show_legend = plot_anno_legend_mut,
      border = T, show_annotation_name = F)
    
    mut_status_bottom <- rep(ref_anova_mut_status, length(mut_status_top))
    bottom_anno <- HeatmapAnnotation(
      'Mutation Status' = mut_status_bottom, 
      col = list('Mutation Status' = structure(df_col_mut_status$colors, names=df_col_mut_status$mut_status)),
      annotation_legend_param = annotation_lgd_list[c(plot_anno_legend_mut, F, F)],
      show_legend = FALSE,
      border = T, show_annotation_name = F)
    
    ht_anova_mut_status <- make_data_heatmap(mat_list$anova_mut_status, 
                                             'ANOVA_mut_status', 1, 
                                             plot_main_legend,
                                             top_anno, bottom_anno)
    plot_main_legend <- FALSE
    plot_anno_legend_mut <- FALSE
  } else {
    ht_anova_mut_status <- NULL
  }
  
  ht_list <- (ht_anova_mut_status + (ht_anova_condition + (ht_gsea_mut_status + (ht_gsea_condition + ht))))
  
  return(invisible(ht_list))
}

make_ht <- function(dist_mat, clusters, df_gsea, df_anova, ht_settings) {
  custom_ht_clusters(
    dist_mat, clusters, df_gsea, df_anova,
    ref_gsea_condition = ht_settings$ref_gsea_cond,
    ref_gsea_mut_status = ht_settings$ref_gsea_mut,
    ref_anova_condition = ht_settings$ref_anova_cond,
    ref_anova_mut_status = ht_settings$ref_anova_mut,
    reduce_matrix = ht_settings$reduce,
    draw_word_cloud = FALSE)
}
