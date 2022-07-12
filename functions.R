library(tidyverse)

get_go_term <- function(go_id, df_annotation) {
  left_join(tibble(id = go_id), df_annotation, by = 'id') %>% 
    mutate(process = replace_na(process, 'NA')) %>%
    pull(process)
}