[![](https://img.shields.io/badge/Shiny-shinyapps.io-blue?style=flat&labelColor=white&logo=RStudio&logoColor=blue)](https://pjunk.shinyapps.io/kras_apms_vis/)
![GitHub](https://img.shields.io/github/license/philippjunk/kras_apms_vis)

# Interactive visualization for KRAS APMS data

This repository contains the source code for a R-shiny app interactively visualizing the data and functional analysis of the KRAS interactome in different genetic and culture contexts. Briefly, Caco-2 cells were transfected with different flag-KRAS plasmids (WT, G12C, G12D, G12V) and cultured with different stimuli (unstimulated, DMOG, EGF, IL6, PGE2, TNFa) at different concentrations. Affinity purification mass spectroscopy was performed to determine the interactome. 

Functional analysis with the aim of linking changes in the interactome to functional outcomes was performed in two different ways: A) differential interaction analysis using `limma` between meaningful pairwise contrasts, followed by a gene set enrichment analysis; B) summing up intensities on functional ontology terms (GO Biological Process) and then analysing differences by ANOVAs and Tukey post hoc tests. The outputs of both analysis were combined by a semantic similarity analysis using `SimplifyEnrichment`. This shiny app is an interactive access to the results of these different analysis as well as the underlying data.

For more details, please refer to our manuscript: COMING SOON

The raw data can be accessed here: COMING SOON

Our data analysis pipeline is available on Zenodo: COMING SOON

# Accessing the app

## Run app remotely

There is a hosted version of the app available [here](https://pjunk.shinyapps.io/kras_apms_vis/) or via the shinyapps.io badge on the top of the README. 

However, **some functionality had to be removed** for the app to run on the available resources. Specifically, making the overview heatmap fully interactive had to be removed. If you want to run the app with the full functionality, please run it locally on your computer with minimal effort:

## Run app locally

Running the app with the full functionally on a local computer is straightforward. First, you need to install `R version >= 4.2`. `R` can be installed from this [web site](https://www.r-project.org/).

Then, in `R`, you need to install the following dependencies:

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

Afterwards, the app can easily be run from your local `R` session:

```R
shiny::runGitHub('philippjunk/kras_apms_vis')
```
