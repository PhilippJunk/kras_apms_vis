##
## Interactive Visualization KRAS APMS
## Based on data by Camille Ternet
## Developed by Philipp Junk, 2022
##

library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(tidyverse)
library(shiny)
library(shinyjs)
library(shinyalert)

# helper function for heatmap generation
source('functions.R')

# action function for InteractiveComplexHeatmap: 
# need to be defined here
click_action = function(df, output) {
  output[["go_info"]] = renderUI({
    if(!is.null(df)) {
      go_id1 = rownames(full_dist_mat)[df$row_index]
      go_id2 = colnames(full_dist_mat)[df$column_index]

      HTML(str_glue(
        "<pre>
## Row GO ID
<a href='http://amigo.geneontology.org/amigo/term/{go_id1}' target='_blank'>{go_id1}</a>: {get_go_term(go_id1)} <button id='{go_id1}' class='go_sel_button'>Select</button>

## Column GO ID:
<a href='http://amigo.geneontology.org/amigo/term/{go_id2}' target='_blank'>{go_id2}</a>: {get_go_term(go_id2)} <button id='{go_id2}' class='go_sel_button'>Select</button>
</pre>"
      ))
    }
  })
}

brush_action = function(df, output) {
  output[["go_info"]] = renderUI({
    if(!is.null(df)) {
      row_index = unique(unlist(df$row_index))
      column_index = unique(unlist(df$column_index))
      go_id1 = rownames(full_dist_mat)[row_index]
      go_id2 = colnames(full_dist_mat)[column_index]

      go_id = union(go_id1, go_id2)

      go_text = str_glue("<a href='http://amigo.geneontology.org/amigo/term/{go_id}' target='_blank'>{go_id}</a>: {get_go_term(go_id)} <button id='{go_id}' class='go_sel_button'>Select</button>") %>%
        str_c(collapse='\n')
      HTML(str_glue("<pre>{go_text}</pre>"))
    }
  })
}


# TODO look into HDF5 files for saving big data sets; or Rdata files

# TODO import and pre-process data

# what data do I need:
# - APMS data
# - APMS data summed up on GO processes
# - GO processes
# - GO process annotation
# - Statistics from ANOVA/Tukey

# setwd("/home/junkpp/work/other/2022-05-20_KRAS_APMS/kras_apms_vis/")

df_apms <- read.csv('data/df_apms.csv', header = T)
df_sum <- read.csv('data/df_sum.csv', header = T)
df_ontology <- read.csv('data/df_ontology.csv', header = T)
df_annotation <- read.csv('data/df_annotation.csv', header = T)
df_anova <- read.csv('data/df_anova.csv', header = T)
df_gsea <- read.csv('data/df_gsea.csv', header = T)

# TODO PUT OTHER DATA OBJECTS IN HERE AS WELL, CURRENTLY ONLY 
# full_dist_mat and cluster_all from SimplifyEnrichment
load('data/data.Rdata')

# create heatmap for tests
ht <- custom_ht_clusters(
  full_dist_mat, 
  cluster_all, 
  df_gsea = df_gsea,
  df_anova = df_anova,
  ref_gsea_condition = 'unstim_none', 
  # ref_anova_condition = 'unstim',
  # ref_gsea_mut_status = 'wt', 
  # ref_anova_mut_status = 'wt',
  reduce_matrix = T)

# color schemes from Camille's thesis:
# TODO filter out what I don't need, maybe include more
# TODO I need color codes for mutation status
conditions <- c('unstim', 'dmog', 'egf', 'il6', 'pge2', 'tnfa')
conditions_hq <- c('unstim.', 'DMOG', 'EGF', 'IL-6', 'PGE2', 'TNF\u03B1')
condition_colors <- c("#001524","#12616D","#75964A",
                      "#A1869E","#FF7D00","#78290F")

mut_status <- c('wt', '12c', '12d', '12v')
mut_status_hq <- c('WT', 'G12C', 'G12D', 'G12V')
mut_status_colors <- c() # TODO

# get choices for inputs
choices_conditions <- df_apms %>%
  select(concentration, condition) %>%
  distinct %>%
  mutate(choice = case_when(
    condition == 'unstim' ~ condition,
    T ~ str_c(condition, concentration, sep='_'))) %>%
  pull(choice)
  
choices_anova_factors <- df_anova %>%
  group_by(term) %>%
  summarise(n_components = str_count(unique(term), ':')) %>%
  ungroup %>%
  arrange(n_components, term) %>% 
  pull(term)

# set default settings for initializing
default_ontology <- 'BP'
default_seltype <- 'process'
default_id <- 'GO:0061621'
default_mut_status <- c('WT', 'G12C', 'G12D', 'G12V')
default_panel <- 'sum_lfq'
default_anova_factors <- c('condition', 'mut_status', 'concentration')
default_anova_pval <- 0.05
default_gsea_pval <- 0.05

# Define UI for application that draws a histogram
ui <- fluidPage(
  useShinyjs(),
  tags$script(src = "custom_button.js"),
  
  # Application title
  titlePanel("KRAS APMS Data"),
  
  ###################################
  # TEST HEATMAP
  # fluidRow(
  #   title = "Original heatmap", width = 4, solidHeader = TRUE, status = "primary",
  #   originalHeatmapOutput("ht", title = NULL)
  # ),
  # 
  # fluidRow(
  #   title = "Sub-heatmap", width = 4, solidHeader = TRUE, status = "primary",
  #   subHeatmapOutput("ht", title = NULL)
  # ),
  # 
  # htmlOutput("go_info"),
  # 
  # 
  # hr(),
  ###################################
  
  # Panel with selection and information for processes
  fluidRow(
    column(
      6,
      h4('Process Information'),
      htmlOutput('ontology_info')
    ),
    column(
      3,
      h4('Ontology Control Panel'),
      # Radio buttons for which type of ontology to use
      radioButtons(
        inputId = 'ontology',
        label = 'Ontology', 
        selected = default_ontology,
        choiceValues = c('BP',
                         'SysGO'),
        choiceNames = c('GO - Biological Process',
                        'SysGO')
      ),
      # Radio buttons on how to select GO terms
      radioButtons(
        inputId = 'seltype',
        label = 'Selection Type',
        choiceValues = c('id', 'process'),
        choiceNames = c('ID', 'Process'),
        selected = default_seltype,
        inline = TRUE
      ),
      # Input for selection of GO terms
      selectInput(
        inputId = 'id',
        label = 'ID/Process',
        choices = c(default_id),
        selected = default_id,
        multiple = FALSE
      )
    ),
    column(
      3,
      h4('Data Control Panel'),
      # Input for selection of mutation status
      selectInput(
        inputId = 'mut_status',
        label = 'Selection Mutation Status',
        choices = c('WT', 'G12C', 'G12D', 'G12V'),
        selected = default_mut_status,
        multiple = TRUE
      )
	  # TODO maybe rebuild?
    )
  ),
  
  hr(),
  
  # main information output
  tabsetPanel(
    id = 'main_panel',
    selected = default_panel,
    # Overview GO processes (for now)
    tabPanel(
      'GO SUMMARY HEATMAP',
      value = 'go_summary_heatmap',
      fluidRow(title = "Original heatmap", originalHeatmapOutput("ht", title = NULL)),
      fluidRow(title = "Sub-heatmap", subHeatmapOutput("ht", title = NULL)),
      shinyjs::hidden(fluidRow('OutputPanel', HeatmapInfoOutput('ht', title=NULL))),
      htmlOutput("go_info")
    ),
    
    # Visualization GO Process
    tabPanel(
      'GO Process',
      value = 'go_process',
      plotOutput('plot_goprocess_info')
    ),
    
    # Visualization sum LFQ
    tabPanel(
      'Sum LFQ Intensity', 
      value = 'sum_lfq',
      plotOutput('plot_lfqsum')
    ),
    
    # Visualization overview proteins
    tabPanel(
      'Individual Proteins 1',
      values = 'proteins_overview',
      plotOutput('plot_proteins_overview')
    ),
    
    # Visulization individual proteins 
    tabPanel(
      'Individual Proteins',
      value = 'proteins',
      sidebarLayout(
        mainPanel(
          plotOutput('plot_proteins')
        ),
        sidebarPanel(
          selectizeInput(
            inputId = 'indiv_proteins',
            label = 'Select proteins to plot',
            choices = character(),
            options = list(maxItems = 10))
          )
        )
    ),
    
    # Table of ANOVA/Tukey
    tabPanel(
      'Statistics',
      value = 'stats',
      sidebarLayout(
        mainPanel(
          dataTableOutput('table_anova')
        ),
        sidebarPanel(
          sliderInput(
            inputId = 'anova_pval',
            label = 'Set cutoff for adjusted p-value',
            min = 0.01,
            max = 0.1,
            value = default_anova_pval
          ),
          selectInput(
            inputId = 'anova_factors',
            label = 'Select terms to display in table',
            choices = choices_anova_factors,
            selected = default_anova_factors,
            multiple = TRUE
          )
        )
      )
    ),
    
    # Table of GSEA
    tabPanel(
      'GSEA',
      value = 'gsea',
      sidebarLayout(
        mainPanel(
          dataTableOutput('table_gsea')
        ),
        sidebarPanel(
          sliderInput(
            inputId = 'gsea_pval',
            label = 'Set cutoff for adjusted p-value',
            min = 0.01,
            max = 0.1,
            value = default_gsea_pval
          )
        )
      )
    )
  ),
    
  hr(),
   
  # additional features
  fluidRow(
    column(
       1,
       actionButton('button_saveplot', 'Save Plot'),
       offset = 1
    ),
    column(
       1, 
       actionButton('button_reset', 'Reset')
    ),
    column(
       1,
       actionLink('button_help', 'Help'),
       offset = 7
    ),
    column(
       1,
       actionLink('button_impressum', 'Impressum')
    )
  )
   
)

# Server logic
server <- function(input, output, session) {
  ##################################################################
  # TEST HEATMAP
  
  makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
                                click_action = click_action, brush_action = brush_action)

  ##################################################################
  ## REACTIVE VALUES: data frames
  
  # reactive subset of df_apms 
  dfr_apms <- reactive({
    # return data frame derived from df_apms filtered by
    # - selected GO process
    # - selected mutations
    validate(need(input$id, 'Please select an ontology term.'),
             need(input$mut_status, 'Please select at least one mutation status'))
    df_ontology %>%
      filter(id == input$id) %>%
      select(hgnc) %>%
      left_join(df_apms, by = 'hgnc') %>%
      filter(mut_status %in% str_to_lower(input$mut_status)) 
  })
  
  # reactive subset of df_sum
  dfr_sum <- reactive({
    # return data frame derived from df_sum filtered by
    # - selected GO process
    # - selected mutations
    validate(need(input$id, 'Please select an ontology term.'),
             need(input$mut_status, 'Please select at least one mutation status'))
    df_sum %>%
      filter(id == input$id) %>%
      filter(mut_status %in% str_to_lower(input$mut_status)) 
  })
  
  # reactive subset of df_anova
  dfr_anova <- reactive({
    # returns data frame derived from df_anova filtered by
    # - selected GO process
    # - selected p_value
    # - selected interactions
    validate(need(input$id, 'Please select an ontology term.'))
    df_anova %>%
      filter(id == input$id) %>%
      filter(p_adj <= input$anova_pval) %>%
      filter(term %in% input$anova_factors) %>%
      arrange(factor(term, levels = input$anova_factors)) %>%
      group_by(term) %>%
      arrange(p_adj, .by_group = T) %>%
      ungroup
  })
  
  # reactive subset of df_gsea
  dfr_gsea <- reactive({
    # return data frame derived from df_gsea filtered by
    # - selected GO process
    # - selected p_value
    validate(need(input$id, 'Please select an ontology term.'))
    df_gsea %>%
      filter(id == input$id) %>%
      filter(p_adj <= input$gsea_pval) %>%
      group_by(enriched_in) %>%
      arrange(desc(NES), .by_group = T) %>%
      ungroup
  })
  
  ##################################################################
  ## REACTIVE VALUES: plots
  
  # construct plot for information on GO process
  reac_plot_goprocess_info <- reactive({
    # TODO include custom color set
    df <- bind_rows(
      tibble(
        label = 'Whole Set',
        condition = 'whole_set',
        count = df_ontology %>% filter(id == input$id) %>% 
          pull(n_all) %>% head(1)
      ),
      dfr_apms() %>%
        group_by(label,group, mut_status, condition, concentration) %>%
        summarise(count = n()) %>%
        ungroup %>%
        arrange(-count) %>%
        select(label, condition, count)
    ) %>%
      mutate(label = factor(label, levels = unique(label)))
    df %>%
      ggplot(aes(x = label, y = count, fill = condition)) +
      geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha=0.7) +
      scale_fill_manual(values = c(condition_colors, 'white'), 
                        breaks = c(conditions, 'whole_set'), 
                        labels = c(conditions_hq, 'Whole Set')) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      theme_minimal() +
      NULL
  })
  
  # construct plot of sum of LFQ intensities 
  reac_plot_lfqsum <- reactive({
    dfr_sum() %>%
      ggplot(aes(x = mut_status, y = sum_LFQ, fill=mut_status)) +
      geom_boxplot(color = 'black', alpha=0.8) +
      geom_point(size = 2, position = position_jitter(height=0, width=0.2)) +
      facet_grid(cols = vars(condition, concentration)) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      theme_bw(base_size = 15) +
      NULL
  })
  
  # construct plot for overview over individual proteins
  reac_plot_proteins_overview <- reactive({
    dfr_apms() %>%
      group_by(hgnc) %>%
      summarise(count = n()) %>%
      ungroup %>%
      arrange(desc(count)) %>%
      mutate(hgnc = factor(hgnc, levels = hgnc)) %>% 
      ggplot(aes(x = hgnc, y = count)) +
      geom_bar(stat = 'identity', color = 'black', fill = 'gray') +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      labs(x = 'HGNC', y = 'Number of samples',
           caption = 'Showing X of X proteins.') + # TODO
      theme_minimal() +
      NULL
  })
  
  # construct plot for individual proteins
  reac_plot_proteins <- reactive({
    validate(need(input$indiv_proteins, 'Please select at least one protein to plot.'))
    plot <- dfr_apms() %>%
      filter(hgnc %in% input$indiv_proteins) %>%
      mutate(hgnc = factor(hgnc, levels = input$indiv_proteins)) %>%
      ggplot(aes(x = str_glue('{condition}_{concentration}'), 
                 y = log2(LFQ), fill = condition)) +
      geom_boxplot(color = 'black', alpha=0.8) +
      geom_point(position = position_jitter(height=0, width=0.2)) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      scale_fill_manual(values = condition_colors, 
                        breaks = conditions, 
                        labels = conditions_hq) +
      theme_bw() +
      NULL
    
    # TODO also adjust based on number of mutations
    # adjust faceting depending on number of proteins to show
    if (length(input$indiv_proteins) > 1) {
      plot <- plot + 
        facet_grid(rows=vars(hgnc), cols = vars(mut_status),
                   scales = 'free')
    } else {
      plot <- plot + 
        facet_grid(cols = vars(mut_status),
                   scales = 'free')
    }
    plot
  })
  
  ##################################################################
  ## REACTIVE VALUES: other
  
  # TODO
  # choices for input$id
  # reac_input_id_choices <- reactive({
  #    
  # })
  
  ##################################################################
  ## RENDERED ELEMENTS
  
  # render information about ontology term currently displayed
  output$ontology_info <- renderUI({
    annotation <- df_annotation %>%
      filter(id == input$id) %>% head(1)
    n_all <- df_ontology %>% 
      filter(id == input$id) %>%
      pull(n_all) %>% head(1)
    n_found <- df_ontology %>% 
      filter(id == input$id) %>%
      pull(n_found) %>% head(1)
    HTML(str_c(str_c('<strong>', annotation$id, '</strong>'),
               str_c('<strong>', annotation$process, '</strong>'),
               annotation$definition,
               '<hr>',
               str_c('Size whole set: ', n_all),
               str_c('Number found in APMS: ', n_found),
               sep = '<br/>'))
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
    dfr_anova() %>%
      select(term, group1, group2, p_adj) %>%
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
      select(id, enriched_in, enriched_against, NES, p_adj) %>%
      rename('Adj. P-value' = p_adj)
  },
  options = list(
    paging = FALSE
  ))
  
  ##################################################################
  ## OBSERVERS
  
  # add observers to selection of processes
  # update based on process selection
  observe({
    df_temp <- df_annotation %>%
      filter(ontology == input$ontology) %>%
      select(id, process) %>%
      filter(id %in% unique(df_sum$id)) %>%
      distinct
    
    choice_values <- df_temp %>% pull(id)
    if (input$seltype == 'process') {
      choice_names <- str_c(choice_values, 
                            df_temp %>% pull(process),
                            sep=' ')
    } else {
      choice_names <- df_temp %>% pull(id)
    }
    names(choice_values) <- choice_names
     
    # keep selection if possible
    # get currently selected
    selected_curr <- input$id
    if (selected_curr %in% df_temp$id) {
      # retain selection
      selected <- df_temp %>%
        filter(id == selected_curr) %>%
        pull(id)
    } else if (selected_curr %in% df_temp$process) {  # TODO is this necessary?
      # retain selection
      selected <- df_temp %>% 
        filter(process == selected_curr) %>%
        pull(id)
    } else {
      # random selection
      selected <- choice_values[1]
    }

    updateSelectInput(
      inputId = 'id',
      choices = choice_values,
      selected = selected,
    )
  })
  
  # add observation to selection of individual proteins to plot
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
  
  # add triggers to bottom row of action buttons/links
  # TODO reset and savePlot
  observeEvent(input$button_reset, {
    # reset everything back to defaults
    updateRadioButtons(inputId = 'ontology', selected = default_ontology)
    updateRadioButtons(inputId = 'seltype', selected = default_seltype)
    updateSelectizeInput(inputId = 'id', selected = default_id)
    updateSelectInput(inputId = 'mut_status', selected = default_mut_status)
  }) 
  
  # for modals, consider html = TRUE and maybe custom icons?
  observeEvent(input$button_help, {
    shinyalert(
      title = 'Help Page',
      text = 'TODO',
      type = 'info'
    )
  })
  
  observeEvent(input$button_impressum, {
    shinyalert(
      title = 'Impressum',
      text = 'TODO',
      type = 'info'
    )
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

