![GitHub](https://img.shields.io/github/license/philippjunk/kras_apms_vis)
[![](https://img.shields.io/badge/Shiny-shinyapps.io-blue?style=flat&labelColor=white&logo=RStudio&logoColor=blue)](https://pjunk.shinyapps.io/kras_apms_vis/)

# Interactive visualization for KRAS APMS data


# Accessing the app

## Run app remotely

There is a hosted version of the app available [here](https://pjunk.shinyapps.io/kras_apms_vis/) or via the shinyapp-io badge on the top of the page. 

However, **some functionality had to be removed** for the app to run on the available resources. Specifically, making the overview heatmap fully interactive had to be removed. If you want to run the app with the full functionality, please run it locally on your computer with minimal effort:

## Run app locally

Running the app with the full functionally on a local computer is straightforward. First, you need to install `R version >= 4.2`. R can be installed from this [web site](https://www.r-project.org/).

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
shiny::runGitHub('philippjunk/kras_apms_vis')
```
