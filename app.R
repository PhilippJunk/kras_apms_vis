##
## Interactive Visualization KRAS APMS
##

# TODO next steps
# - try caching heatmap to potentially improve performance
# - write help (either popup, or maybe a couple readthedocs pages?)

library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(tidyverse)
library(ggwordcloud)
library(shiny)
library(shinydashboard)
library(shinyjs)

# necessary for upload to shinyapps.io
# library(BiocManager)
# options(repos = BiocManager::repositories())

# helper function for heatmap generation
source('functions.R')

# load data
load('data/data.Rdata')

# action function for InteractiveComplexHeatmap: 
# needs to be defined here
click_action = function(df, output) {
  output[["go_info_heatmap"]] = renderUI({
    if(!is.null(df)) {
      df <- as_tibble(df)
      
      # distinguish between click on similarity heatmap and click on data heatmaps
      if (df$heatmap == 'Similarity') {
        go_id1 = df$row_label
        go_id2 = df$column_label
        
        HTML(str_glue(
          "<pre>
## Row GO ID
<a href='http://amigo.geneontology.org/amigo/term/{go_id1}' target='_blank'>{go_id1}</a>: {get_go_term(go_id1, df_annotation)} <button id='{go_id1}' class='go_sel_button'>Select</button>

## Column GO ID:
<a href='http://amigo.geneontology.org/amigo/term/{go_id2}' target='_blank'>{go_id2}</a>: {get_go_term(go_id2, df_annotation)} <button id='{go_id2}' class='go_sel_button'>Select</button>
</pre>"
        ))
      } else {
        go_id1 = df$row_label

        HTML(str_glue(
          "<pre>
## Row GO ID
<a href='http://amigo.geneontology.org/amigo/term/{go_id1}' target='_blank'>{go_id1}</a>: {get_go_term(go_id1, df_annotation)} <button id='{go_id1}' class='go_sel_button'>Select</button>
</pre>"
        ))
      }
      

    }
  })
}

brush_action = function(df, output) {
  output[["go_info_heatmap"]] = renderUI({
    if(!is.null(df)) {
      df <- as_tibble(df)
      # extract GO IDs from brushed area
      go_id <- c(df %>% filter(heatmap == 'Similarity') %>%
                   select(row_label, column_label) %>% 
                   as.list %>% 
                   flatten %>% 
                   flatten_chr,
                 df %>% filter(heatmap != 'Similarity') %>%
                   select(row_label) %>% 
                   as.list %>%
                   flatten %>% 
                   flatten_chr)

      go_text = str_glue("<a href='http://amigo.geneontology.org/amigo/term/{go_id}' target='_blank'>{go_id}</a>: {get_go_term(go_id, df_annotation)} <button id='{go_id}' class='go_sel_button'>Select</button>") %>%
        str_c(collapse='\n')
      HTML(str_glue("<pre>{go_text}</pre>"))
    }
  })
}

# color schemes 
conditions <- c('unstim', 'dmog', 'egf', 'il6', 'pge2', 'tnfa')
conditions_hq <- c('unstim.', 'DMOG', 'EGF', 'IL-6', 'PGE2', 'TNF\u03B1')
conditions_colors <- c("#001524","#12616D","#75964A","#A1869E","#FF7D00","#78290F")

mut_stats <- c('wt', 'g12c', 'g12d', 'g12v')
mut_stats_hq <- c('WT', 'G12C', 'G12D', 'G12V')
mut_stats_colors <- c('#A98743', '#437C90', '#255957', '#8F2844')

# get choices for inputs
choices_clusters <- tibble(cluster = clusters) %>%
  count(cluster, name = 'count') %>% 
  arrange(-count) %>%
  filter(count >= 20) %>%
  pull(cluster)
  
choices_mut_status <- df_apms %>%
  pull(mut_status) %>%
  unique %>%
  set_names(nm = str_to_upper(.))

choices_condition <- df_apms %>% 
  pull(condition) %>%
  unique %>%
  set_names(nm = case_when(. != 'unstim' ~ str_to_upper(.), 
                           T ~ str_to_title(.)))

choices_cond_conc <- df_apms %>%
  select(concentration, condition) %>%
  distinct %>%
  mutate(choice = str_c(condition, concentration, sep='_')) %>%
  pull(choice) %>%
  set_names(nm = case_when(. != 'unstim_none' ~ str_to_upper(.), 
                           T ~ 'Unstim'))
  
choices_anova_factors <- df_anova %>%
  group_by(term) %>%
  summarise(n_components = str_count(unique(term), ':')) %>%
  ungroup %>%
  arrange(n_components, term) %>% 
  pull(term)

# set default settings for initializing
default_seltype <- 'process'
default_id <- 'GO:1903578'
default_mut_status <- choices_mut_status
default_cond_conc <- choices_cond_conc
default_anova_factors <- c('condition', 'mut_status', 'concentration')
default_anova_pval <- 0.05
default_gsea_pval <- 0.05


# TODO think about potential icons??
ui <- dashboardPage(
  dashboardHeader(
    title = 'KRAS APMS Visualization',
    dropdownMenu(type = 'notifications', headerText = 'Further links', 
                 icon = icon('info'), badgeStatus = NULL, # TODO add link to article
                 notificationItem(text = 'Source code', icon = icon('github'), href = 'https://github.com/PhilippJunk/kras_apms_vis'), #TODO link
                 notificationItem(text = 'Contact Scientific', icon = icon('envelope'), href = NULL),
                 notificationItem(text = 'Contact Technical', icon = icon('envelope'), href = 'mailto:philipp.junk@ucdconnect.ie?subject=Shiny App KRAS APMS'))), # TODO link
  dashboardSidebar(
    sidebarMenu(
      menuItem(text = 'Overview',
               tabName = 'overview',
               icon = NULL,
               menuSubItem(text = 'Heatmap',
                           tabName = 'overview-heatmap',
                           icon = NULL),
               menuSubItem(text = 'Semantic Clusters',
                           tabName = 'overview-clusters',
                           icon = NULL)
               ),
      menuItem(text = 'Invidivial GO Term',
               tabName = 'specific',
               icon = NULL,
               menuSubItem(text = 'Summary',
                           tabName = 'specific-summary',
                           icon = NULL),
               menuSubItem(text = 'Summed Intensities',
                           tabName = 'specific-sum',
                           icon = NULL),
               menuSubItem(text = 'Individual Intensities',
                           tabName = 'specfic-indiv',
                           icon = NULL)
               ),
      menuItem(text = 'Settings',
               tabName = 'settings',
               icon = icon('gear', verify_fa = FALSE)),
      menuItem(text = 'Help',
               tabName = 'help',
               icon = icon('question')),
      hr(),
      h4('Ontology Control Panel', style = 'text-align:center'),
      # Radio buttons on how to select GO terms
      radioButtons(inputId = 'seltype', label = 'Selection Type', 
                   choiceValues = c('id', 'process'), 
                   choiceNames = c('ID', 'Process'),
                   selected = default_seltype, inline = TRUE),
      # Input for selection of GO terms
      selectInput(inputId = 'id', label = 'ID/Process', 
                  choices = default_id, selected = default_id,
                  multiple = FALSE),
      # History of selected GO terms
      actionButton('id_history', 'View History', icon = icon('history')),
      # Random selection of GO term
      actionButton('id_random', 'Random GO term', icon = icon('dice')),
      hr(),
      h4('Data Control Panel', style = 'text-align:center'),
      # Input for selection of mutation status
      selectInput(inputId = 'mut_status', label = 'Selection Mutation Status',
                  choices = choices_mut_status, 
                  selected = default_mut_status, multiple = TRUE),
      # Input for selection of condition/concentration 
      selectInput(inputId = 'cond_conc', label = 'Select Condition/Concentration',
                  choices = choices_cond_conc, selected = default_cond_conc, 
                  multiple = TRUE)
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'custom_css.css'),
              tags$script(src = "custom_button.js")),
    tabItems(
      tabItem(tabName = 'overview-heatmap',
              h2('Overview Semantic Distance Heatmap'),
              fluidRow(box(width = 12, title = 'Heatmap Controls',
                           collapsible = TRUE, collapsed = TRUE,
                           selectizeInput(inputId = 'ht_ref_anova_mut',
                                          label = 'Reference group: ANOVA mutation status',
                                          choices = choices_mut_status, selected = NULL,
                                          options = list(maxItems = 1, 
                                                         onInitialize = I('function() { this.setValue(""); }'))),
                           selectizeInput(inputId = 'ht_ref_anova_cond',
                                          label = 'Reference group: ANOVA condition',
                                          choices = choices_condition, selected = NULL,
                                          options = list(maxItems = 1, 
                                                         onInitialize = I('function() { this.setValue(""); }'))),
                           selectizeInput(inputId = 'ht_ref_gsea_mut',
                                          label = 'Reference group: GSEA mutation status',
                                          choices = choices_mut_status, selected = NULL,
                                          options = list(maxItems = 1, 
                                                         onInitialize = I('function() { this.setValue(""); }'))),
                           selectizeInput(inputId = 'ht_ref_gsea_cond',
                                          label = 'Reference group: GSEA condition/concentration',
                                          choices = choices_cond_conc, 
                                          selected = 'unstim_none', 
                                          options = list(maxItems = 1)),
                           hr(),
                           radioButtons(inputId = 'ht_reduce', 
                                        label = 'Only show GO terms that are signficantly different in the data shown in the heatmap?',
                                        choiceValues = c(TRUE, FALSE), choiceNames = c('Yes', 'No'), 
                                        selected = TRUE, inline = TRUE),
                           radioButtons(inputId = 'ht_clusters', label = 'Show all clusters?',
                                        choiceValues = c(TRUE, FALSE), choiceNames = c('Yes', 'No'), 
                                        selected = TRUE, inline = TRUE),
                           conditionalPanel(condition = 'input.ht_clusters == "FALSE"',
                                            selectInput(inputId = 'ht_cluster_sel',
                                                        label = 'Select clusters to show',
                                                        choices = choices_clusters,
                                                        selected = 1, multiple = T)),
                           fluidRow(column(2 ,actionButton(inputId = 'ht_apply', label = 'Render heatmap')),
                                    column(10, hidden(span(style = 'color:red', id = 'error_ht',
                                                           'Rendering the heatmap with the current settings not possible. Please choose other settings.')))))),
              fluidRow(box(originalHeatmapOutput("ht", title = NULL, width = '1000px'),
                           id = 'orig_ht', width = 12, title = "Full heatmap")),
              fluidRow(box(subHeatmapOutput("ht", title = NULL),
                           id = 'sub_ht', width = 12, title = "Sub-heatmap")),
              shinyjs::hidden(fluidRow(title = 'OutputPanel', HeatmapInfoOutput('ht', title=NULL))),
              fluidRow(box(htmlOutput("go_info_heatmap"),
                           title = 'GO terms in selected area', width = 12, 
                           collapsible = T))),
      tabItem(tabName = 'overview-clusters',
              h2('Overview GO Semantic Clusters'),
              uiOutput('cluster_wc_tabs'),
              fluidRow(box(htmlOutput('go_info_clusters'),
                           title = 'GO terms in selected cluster', width = 12, 
                           collapsible = T))),
      tabItem(tabName = 'specific-summary',
              h2('Summary Selected GO term'),
              fluidRow(box(htmlOutput('ontology_info'),
                           title = 'Process Information', width = 12)),
              fluidRow(box(plotOutput('plot_proteins_overview'),
                           downloadButton('dl_png_proteins_overview', label = 'Download PNG'),
                           downloadButton('dl_csv_proteins_overview', label = 'Download CSV'),
                           title = 'Samples per protein', width = 12)),
              fluidRow(box(plotOutput('plot_goprocess_info'),
                           downloadButton('dl_png_goprocess_info', label = 'Download PNG'),
                           downloadButton('dl_csv_goprocess_info', label = 'Download CSV'),
                           title = 'Number of identified proteins', width = 12))),
      tabItem(tabName = 'specific-sum',
              h2('Summed LFQ intensities'),
              fluidRow(box(plotOutput('plot_lfqsum'),
                           downloadButton('dl_png_lfqsum', label = 'Download PNG'),
                           downloadButton('dl_csv_lfqsum', label = 'Download CSV'),
                           title = 'Summed LFQ intensities', width = 12)),
              fluidRow(box(
                title = 'Differential Analysis: ANOVA', width = 12, collapsible = T,
                sidebarLayout(mainPanel(dataTableOutput('table_anova')),
                              sidebarPanel(sliderInput(inputId = 'anova_pval', 
                                                       label = 'Set cutoff for adjusted p-value',
                                                       min = 0.01, max = 0.1,
                                                       value = default_anova_pval),
                                           selectInput(inputId = 'anova_factors',
                                                       label = 'Select terms to display in table',
                                                       choices = choices_anova_factors,
                                                       selected = default_anova_factors,
                                                       multiple = TRUE),
                                           downloadButton('dl_csv_anova', label = 'Download CSV'))))),
              fluidRow(box(
                title = 'Differential Analysis: GSEA', width = 12, collapsible = T,
                sidebarLayout(mainPanel(dataTableOutput('table_gsea')),
                              sidebarPanel(sliderInput(inputId = 'gsea_pval',
                                                       label = 'Set cutoff for adjusted p-value',
                                                       min = 0.01, max = 0.1, 
                                                       value = default_gsea_pval),
                                           downloadButton('dl_csv_gsea', label = 'Download CSV')))))),
      tabItem(tabName = 'specfic-indiv',
              h2('Individual Proteins LFQ intensities'),
              fluidRow(box(
                sidebarLayout(mainPanel(plotOutput('plot_proteins')),
                              sidebarPanel(selectizeInput(inputId = 'indiv_proteins',
                                           label = 'Select proteins to plot',
                                           choices = character(),
                                           options = list(maxItems = 10)))),
                downloadButton('dl_png_proteins', label = 'Download PNG'),
                downloadButton('dl_csv_proteins', label = 'Download CSV'),
                width = 12))),
      tabItem(tabName = 'settings',
              h2('Settings'),
              fluidRow(column(width = 6,
                              box(numericInput('settings_plot_height', 'Set height of plot downloads (in mm)', 
                                               value = 150, max = 300),
                                  numericInput('settings_plot_width', 'Set width of plot downloads (in mm)', 
                                               value = 500, max = 1000),
                                  numericInput('settings_plot_dpi', 'Set resolution of plot downloads (DPI)', 
                                               value = 100, max = 300),
                                  width = NULL, title = 'Download settings')),
                       column(width = 6,
                              box(numericInput('settings_n_history', 'Set number of elements shown',
                                               value = 20, min=1),
                                  width = NULL, title = 'History Settings'),
                              box(numericInput('settings_n_proteins', 'Set number of proteins shown',
                                               value = 50, min=1),
                                  width = NULL, title = 'Settings "Samples per Proteins"')))),
      tabItem(tabName = 'help',
              h2('Help'))
    )
  )
)

# Server logic
server <- function(input, output, session) {
  ##################################################################
  ## REACTIVE VALUES: general data frames
  
  # reactive subset of df_apms 
  dfr_apms <- reactive({
    # return data frame derived from df_apms filtered by
    # - selected GO process
    # - selected mutations
    # - selected combinations of conditions/concentrations
    validate(need(input$id, 'Please select an ontology term.'),
             need(input$mut_status, 'Please select at least one mutation status'),
             need(input$cond_conc, 'Please select at least one condition'))
    df_ontology %>%
      filter(id == input$id) %>%
      select(hgnc) %>%
      left_join(df_apms, by = 'hgnc') %>%
      filter(mut_status %in% input$mut_status) %>%
      filter(str_glue('{condition}_{concentration}') %in% input$cond_conc)
  })
  
  # reactive subset of df_sum
  dfr_sum <- reactive({
    # return data frame derived from df_sum filtered by
    # - selected GO process
    # - selected mutations
    # - selected combinations of conditions/concentrations
    validate(need(input$id, 'Please select an ontology term.'),
             need(input$mut_status, 'Please select at least one mutation status'),
             need(input$cond_conc, 'Please select at least one condition'))
    df_sum %>%
      filter(id == input$id) %>%
      filter(mut_status %in% input$mut_status) %>%
      filter(str_glue('{condition}_{concentration}') %in% input$cond_conc)
  })
  
  # reactive subset of df_anova
  dfr_anova <- reactive({
    # returns data frame derived from df_anova filtered by
    # - selected GO process
    # - selected p_value
    # - selected interactions
    # TODO filtering by data input panel
    validate(need(input$id, 'Please select an ontology term.'))
    df_anova %>%
      filter(id == input$id) %>%
      filter(p_adj <= input$anova_pval) %>%
      filter(term %in% input$anova_factors) %>%
      arrange(factor(term, levels = input$anova_factors)) %>%
      group_by(term) %>%
      arrange(p_adj, .by_group = T) %>%
      ungroup %>%
      select(term, group_higher, group_lower, estimate, p_adj)
  })
  
  # reactive subset of df_gsea
  dfr_gsea <- reactive({
    # return data frame derived from df_gsea filtered by
    # - selected GO process
    # - selected p_value
    # TODO filtering by data input panel
    validate(need(input$id, 'Please select an ontology term.'))
    df_gsea %>%
      filter(id == input$id) %>%
      filter(p_adj <= input$gsea_pval) %>%
      group_by(enriched_in) %>%
      arrange(desc(NES), .by_group = T) %>%
      ungroup %>%
      select(id, enriched_in, enriched_against, NES, p_adj)
  })
  
  ##################################################################
  ## REACTIVE VALUES: data frames for plots
  
  dfr_plot_proteins_overview <- reactive({
    dfr_apms() %>%
      group_by(hgnc) %>%
      summarise(count = n()) %>%
      ungroup %>%
      arrange(desc(count)) %>%
      mutate(hgnc = factor(hgnc, levels = hgnc))
  })
  
  dfr_plot_goprocess_info <- reactive({
    bind_rows(
      tibble(
        group = 'Whole Set',
        condition = 'whole_set',
        count = df_annotation %>% filter(id == input$id) %>% 
          pull(n_all) %>% head(1)
      ),
      dfr_apms() %>%
        group_by(group, mut_status, condition, concentration) %>%
        summarise(count = n()) %>%
        ungroup %>%
        arrange(-count) %>%
        select(group, condition, count)
    ) %>%
      mutate(group = factor(group, levels = unique(group)))
  })
  
  dfr_plot_proteins <- reactive({
    validate(need(input$indiv_proteins, 'Please select at least one protein to plot.'))
    dfr_apms() %>%
      filter(hgnc %in% input$indiv_proteins) %>%
      mutate(hgnc = factor(hgnc, levels = input$indiv_proteins)) %>%
      mutate(mut_status = str_to_upper(mut_status))
  })
  
  
  
  ##################################################################
  ## REACTIVE VALUES: plots
  
  # heatmap
  ht <- reactive({
    make_ht(dist_mat, clusters, df_gsea, df_anova, ht_settings())
  })
  
  # construct plot for overview over individual proteins
  reac_plot_proteins_overview <- reactive({
    validate(need(input$settings_n_proteins, 'Please provide the number of proteins to show in Settings.'))
    df <- dfr_plot_proteins_overview()
    n_total <- df$hgnc %>% unique %>% length
    df <- df %>% 
      slice_max(count, n=settings_n_proteins())
    n_here <- min(settings_n_proteins(), n_total)
    
    ggplot(df, aes(x = hgnc, y = count)) +
      geom_bar(stat = 'identity', color = 'black', fill = 'gray') +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      labs(x = 'HGNC', y = 'Number of samples',
           caption = str_glue('Showing {n_here} of {n_total} proteins.')) + 
      theme_minimal() +
      theme(plot.background = element_rect(fill = 'white', color = 'white')) +
      NULL
  })
  
  # construct plot for information on GO process
  reac_plot_goprocess_info <- reactive({
    # TODO include custom color set
    dfr_plot_goprocess_info() %>%
      ggplot(aes(x = group, y = count, fill = condition)) +
      geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha=0.7) +
      scale_fill_manual(values = c(conditions_colors, 'white'), 
                        breaks = c(conditions, 'whole_set'), 
                        labels = c(conditions_hq, 'Whole Set')) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      labs(x = 'Group', y = 'Number of proteins', full = 'Condition') +
      theme_minimal() +
      theme(plot.background = element_rect(fill = 'white', color = 'white')) +
      NULL
  })
  
  # construct plot of sum of LFQ intensities 
  reac_plot_lfqsum <- reactive({
    dfr_sum() %>%
      # mutate(mut_status = str_to_upper(mut_status)) %>%
      mutate(condition = case_when(condition == 'unstim' ~ condition,
                                   T ~ str_to_upper(condition))) %>%
      ggplot(aes(x = mut_status, y = sum_LFQ, fill=mut_status)) +
      geom_boxplot(color = 'black', alpha=0.8) +
      geom_point(size = 2, position = position_jitter(height=0, width=0.2)) +
      facet_grid(cols = vars(condition, concentration)) +
      scale_x_discrete(guide = guide_axis(angle = 45),
                       breaks = mut_stats,
                       labels = mut_stats_hq) +
      scale_fill_manual(values = mut_stats_colors,
                        breaks = mut_stats, 
                        labels = mut_stats_hq) + 
      labs(x = 'Mutation Status', y = 'log2(Sum(LFQ))', fill = 'Mutation Status') +
      theme_bw(base_size = 15) +
      NULL
  })
  
  # construct plot for individual proteins
  reac_plot_proteins <- reactive({
    plot <- dfr_plot_proteins() %>%
      ggplot(aes(x = str_glue('{condition}_{concentration}'), 
                 y = log2(LFQ), fill = condition)) +
      geom_boxplot(color = 'black', alpha=0.8) +
      geom_point(position = position_jitter(height=0, width=0.2)) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      scale_fill_manual(values = conditions_colors, 
                        breaks = conditions, 
                        labels = conditions_hq) +
      labs(x = 'Condition/Concentration', y = 'log2(LFQ)', fill = 'Condition') +
      theme_bw() +
      NULL
    
    # adjust faceting depending on number of proteins to show
    if (length(input$indiv_proteins) > 1 & length(input$mut_status) > 1) {
      plot <- plot + facet_grid(hgnc ~ mut_status, scales = 'free')
    } else if (length(input$mut_status) > 1) {
      plot <- plot + facet_grid(. ~ mut_status, scales = 'free')
    } else if(length(input$indiv_proteins) > 1) {
      plot <- plot + facet_grid(hgnc ~ .)
    }
    plot
  })
  
  ##################################################################
  ## REACTIVE VALUES: others
  
  # history of selected ids
  id_history <- reactiveVal(value = NULL)
  
  # choices for id selection
  choices_id <- reactive({
    validate(need(input$mut_status, 'Please select at least one mutation status'),
             need(input$cond_conc, 'Please select at least one condition'))
    df_sum %>%
      filter(mut_status %in% input$mut_status) %>%
      filter(str_glue('{condition}_{concentration}') %in% input$cond_conc) %>%
      pull(id) %>% 
      unique
  })
  
  # manually validate settings because numericInput does not
  settings_plot_height <- reactive({
    max_value <- 300
    min(max_value, input$settings_plot_height)
  })
  
  settings_plot_width <- reactive({
    max_value <- 1000
    min(max_value, input$settings_plot_width)
  })
  
  settings_plot_dpi <- reactive({
    max_value <- 300
    min(max_value, input$settings_plot_dpi)
  })
  
  settings_n_history <- reactive({
    min_value <- 1
    max(min_value, input$settings_n_history)
  })
  
  settings_n_proteins <- reactive({
    min_value <- 1
    max(min_value, input$settings_n_proteins)
  })
  
  ##################################################################
  ## RENDERED ELEMENTS
  
  # render cluster word clouds in tab overview
  output$cluster_wc_tabs <- renderUI({
    tabs <- tibble(cluster = clusters) %>% 
      count(cluster, name = 'count') %>% 
      filter(count >= 20) %>%
      arrange(-count) %>%
      pull(cluster) %>% 
      map(function(cl) {
        n_terms <- sum(clusters == cl)
        p_wordcloud <- df_wordcloud %>% 
          filter(cluster == cl) %>% 
          slice_min(padj, n=20) %>%
          ggplot(aes(label = keyword, size = -log10(padj), color=keyword)) +
          geom_text_wordcloud() +
          scale_size_area(max_size = 30) +
          labs(caption = str_glue('{n_terms} GO terms in cluster {cl}')) +
          theme_minimal() + NULL
        tab_title = str_glue('Cluster {cl}')
        tabPanel(tab_title, value = cl, 
                 renderPlot({p_wordcloud}))
      })
    do.call(tabBox, c(tabs, list(width = 12, title = 'Cluster Wordclouds', id='cluster_wc')))
  })
  
  # render all go terms in currently selected cluster
  output$go_info_clusters <- renderUI({
    # get currently active tab and extract GO terms
    go_id <- rownames(dist_mat)[clusters == as.numeric(input$cluster_wc)]
    if(length(go_id > 0)) {
      go_text = str_glue("<a href='http://amigo.geneontology.org/amigo/term/{go_id}' target='_blank'>{go_id}</a>: {get_go_term(go_id, df_annotation)} <button id='{go_id}' class='go_sel_button'>Select</button>") %>%
        str_c(collapse='\n')
      HTML(str_glue("<pre>{go_text}</pre>"))
    }
  })
  
  # render information about ontology term currently displayed
  output$ontology_info <- renderUI({
    # get related GO terms
    relatives <- tibble(dist = dist_mat[rownames(dist_mat) == input$id,], 
                        id = colnames(dist_mat)) %>% 
      arrange(-dist) %>% 
      filter(id != input$id) %>% head(5) %>% pull(id)
    relatives_text = str_glue("<a href='http://amigo.geneontology.org/amigo/term/{relatives}' target='_blank'>{relatives}</a>: {get_go_term(relatives, df_annotation)} <button id='{relatives}' class='go_sel_button'>Select</button>") %>%
        str_c(collapse='\n')
    col2 <- HTML(
      h5('Closest relatives in data set') %>% as.character, 
      str_glue("<pre>{relatives_text}</pre>"))
    
    # get information
    annotation <- df_annotation %>%
      filter(id == input$id) %>% head(1)
    col1 <- HTML(
      h4(annotation$id) %>% as.character,
      h4(annotation$process) %>% as.character,
      p(annotation$definition) %>% as.character,
      hr() %>% as.character,
      p('Size whole set: ', annotation$n_all) %>% as.character,
      p('Number found in APMS: ', annotation$n_found) %>% as.character)

    # assemble UI
    fluidRow(column(6, col1),
             column(6, col2))
  })
  
  # render plot for information on GO process
  output$plot_goprocess_info <- renderPlot({
    reac_plot_goprocess_info()
  })
  
  # render plot of sum of LFQ intensities 
  output$plot_lfqsum <- renderPlot({
    reac_plot_lfqsum()
  })

  # render plot for overview over individual proteins
  output$plot_proteins_overview <- renderPlot({
    reac_plot_proteins_overview()
  })
  
  # render plot for individual proteins
  output$plot_proteins <- renderPlot({
    reac_plot_proteins()
  })
  
  # render ANOVA table
  output$table_anova <- renderDataTable({
    # TODO potentially rename estimate into something more meaningful
    dfr_anova() %>%
      mutate(estimate = round(estimate, 2)) %>%
      mutate(p_adj = scales::scientific(p_adj)) %>%
      rename('Adj. P-value' = p_adj)
    },
    options = list(
      paging = FALSE
    )
  )
  
  # render GSEA table
  output$table_gsea <- renderDataTable({
    dfr_gsea() %>%
      mutate(p_adj = scales::scientific(p_adj)) %>%
      mutate(NES = round(NES, 2)) %>%
      rename('Adj. P-value' = p_adj)
  },
  options = list(
    paging = FALSE
  ))
  
  ##################################################################
  ## DOWNLOAD HANDLERS
  
  # for goprocess_info
  output$dl_png_goprocess_info <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_proteinsPerCondition.png'),
    content = function(con) {ggsave(con, reac_plot_goprocess_info(), device = 'png',
                                    height = settings_plot_height(),
                                    width = settings_plot_width(),
                                    dpi = settings_plot_dpi(),
                                    units = 'mm', limitsize = F)},
    contentType = 'image/png'
  )
  
  output$dl_csv_goprocess_info <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_proteinsPerCondition.csv'),
    content = function(con) {write.csv(dfr_plot_goprocess_info(), con)},
    contentType = 'text/csv'
  )
  
  # for lfqsum
  output$dl_png_lfqsum <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_sumLFQ.png'),
    content = function(con) {ggsave(con, reac_plot_lfqsum(), device = 'png',
                                    height = settings_plot_height(),
                                    width = settings_plot_width(),
                                    dpi = settings_plot_dpi(),
                                    units = 'mm', limitsize = F)},
    contentType = 'image/png'
  )
  
  output$dl_csv_lfqsum <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_sumLFQ.csv'),
    content = function(con) {write.csv(dfr_sum(), con)},
    contentType = 'text/csv'
  )
  
  # for proteins_overview
  output$dl_png_proteins_overview <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_conditionsPerProtein.png'),
    content = function(con) {ggsave(con, reac_plot_proteins_overview(), device = 'png',
                                    height = settings_plot_height(),
                                    width = settings_plot_width(),
                                    dpi = settings_plot_dpi(),
                                    units = 'mm', limitsize = F)},
    contentType = 'image/png'
  )
  
  output$dl_csv_proteins_overview <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_conditionsPerProtein.csv'),
    content = function(con) {write.csv(dfr_plot_proteins_overview(), con)},
    contentType = 'text/csv'
  )
  
  # for proteins
  output$dl_png_proteins <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_indivProteins.png'),
    content = function(con) {ggsave(con, reac_plot_proteins(), device = 'png',
                                    height = settings_plot_height(),
                                    width = settings_plot_width(),
                                    dpi = settings_plot_dpi(),
                                    units = 'mm', limitsize = F)},
    contentType = 'image/png'
  )
  
  output$dl_csv_proteins <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_indivProteins.csv'),
    content = function(con) {write.csv(dfr_plot_proteins(), con)},
    contentType = 'text/csv'
  )
  
  # for anova
  output$dl_csv_anova <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_anova.csv'),
    content = function(con) {write.csv(dfr_anova(), con)},
    contentType = 'text/csv'
  )
  
  # for gsea
  output$dl_csv_gsea <- downloadHandler(
    filename = str_glue('{str_replace(input$id, ":", "")}_gsea.csv'),
    content = function(con) {write.csv(dfr_gsea(), con)},
    contentType = 'text/csv'
  )
  
  ##################################################################
  ## OBSERVERS
  
  # add trigger to extract settings for plotting the heatmap
  # tied to ht_apply action button
  ht_settings <- bindEvent({
    reactive({
      list(ref_anova_mut = if(input$ht_ref_anova_mut != '') {input$ht_ref_anova_mut} else {NULL},
           ref_anova_cond = if(input$ht_ref_anova_cond != '') {input$ht_ref_anova_cond} else {NULL},
           ref_gsea_mut = if(input$ht_ref_gsea_mut != '') {input$ht_ref_gsea_mut} else {NULL},
           ref_gsea_cond = if(input$ht_ref_gsea_cond != '') {input$ht_ref_gsea_cond} else {NULL},
           reduce = as.logical(input$ht_reduce),
           all_clusters = as.logical(input$ht_clusters),
           selected_clusters = as.numeric(input$ht_cluster_sel))
    })
  },
  input$ht_apply,
  ignoreNULL = FALSE)
  
  # catch settings where heatmap returns empty and then does not update
  observe({
    if (!is.null(ht())) {
      shinyjs::show('orig_ht')
      shinyjs::show('sub_ht')
      shinyjs::hide('error_ht')
      
      makeInteractiveComplexHeatmap(
        input, output, session, ht(), "ht",
        click_action = click_action, brush_action = brush_action)
    }
    else {
      shinyjs::hide('orig_ht')
      shinyjs::hide('sub_ht')
      shinyjs::show('error_ht')
    }
  })
  
  # observer for updating history of selected ids
  bindEvent({
    observe({
      id_history(c(input$id, id_history()))
    })
  },
  input$id)
  
  # add observers to selection of processes
  bindEvent({
    observe({
      choices <- choices_id()
      current <- input$id
      
      # add names depending on which seltype has been selected
      if (input$seltype == 'process') {
        choices_names <- left_join(
          tibble(id = choices),
          df_annotation, 
          by = 'id') %>% 
          mutate(name = str_glue('{id} {process}')) %>%
          pull(name)
        names(choices) <- choices_names
      }

      # check if current selection still in choices
      if (current %in% choices) {
        selected <- current
      } else {
        selected <- choices[1]
      }

      updateSelectInput(inputId = 'id', choices = choices, selected = selected)
    })},
    input$seltype,
    choices_id(),
    ignoreNULL = FALSE
  )
  
  # add trigger to id history
  bindEvent({
    observe({
      showModal(modalDialog(
        p(str_glue('Showing the last {settings_n_history()} selected GO IDs.')),
        {
          go_id <- id_history() %>% head(settings_n_history())
          if(length(go_id > 0)) {
            go_text = str_glue("<a href='http://amigo.geneontology.org/amigo/term/{go_id}' target='_blank'>{go_id}</a>: {get_go_term(go_id, df_annotation)} <button id='{go_id}' class='go_sel_button'>Select</button>") %>%
              str_c(collapse='\n')
            HTML(str_glue("<pre>{go_text}</pre>"))
          }},
        title = 'GO ID selection history',
        easyClose = TRUE,
        fade = TRUE,
        size = 'l'))
    })
  },
  input$id_history)
  

  # add trigger to random id selection 
  bindEvent({
    observe({
      updateSelectInput(inputId = 'id', selected = sample(choices_id(), 1))
    })
  },
  input$id_random)
  
  # add observation to selection of individual proteins to visualize
  # Updated based on which GO process is selected
  observe({
    df_proteins <- dfr_apms() %>%
      group_by(hgnc) %>%
      summarise(count = n()) %>%
      ungroup %>%
      arrange(desc(count))
    updateSelectizeInput(
       inputId = 'indiv_proteins',
       choices = df_proteins$hgnc,
       selected = df_proteins$hgnc[1],
       server = TRUE
     )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

