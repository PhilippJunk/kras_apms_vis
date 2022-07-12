# Shiny visualization for KRAS APMS data

## Access app on [shinyapps.io](https://shinyapps.io/)

## Run app

Due to the unfortunate circumstance that the interactive heatmap requires some computational resources to render, deployment on `shinyapps.io` sadly is not an option. However, running the app locally is straightforward. 

First, you need to install `R version >= 4.2`. R can be installed from this [web site]().


Then, in R, you need to install the following dependencies:

```R
# tidyverse
if (!require(tidyverse)) {
  install.packages('tidyverse')
  library(tidyverse)
}

# CRAN
c('ggwordcloud', 'shiny', 'shinydashboard', 'shinyjs', 'shinyalert', 'BiocManager') %>%
  map(function(pkg) {
    if (!require(pkg, character.only = T)) {
      install.packages(pkg)
    }
    library(pkg, character.only = T)
  }) %>% invisible

# Bioconductor
c('ComplexHeatmap', 'InteractiveComplexHeatmap') %>%
  map(function(pkg) {
    if(!require(pkg, character.only = T)) {
      BiocManager::install(pkg)
    }
    library(pkg, character.only = T)
  }) %>% invisible

```

Afterwards, the app can easily be run from your local R session:

```R
shiny::runGithub

```
