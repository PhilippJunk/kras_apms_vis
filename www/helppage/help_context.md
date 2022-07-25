This is an R Shiny app for the interactive visualization of data and functional analysis of the KRAS interactome in different genetic and culture contexts. Briefly, Caco-2 cells were transfected with different flag-KRAS plasmids (WT, G12C, G12D, G12V) and cultured with different stimuli (unstimulated, DMOG, EGF, IL6, PGE2, TNFa) at different concentrations. Affinity purification mass spectroscopy was performed to determine the interactome. 

Functional analysis with the aim of linking changes in the interactome to functional outcomes was performed in two different ways: A) differential interaction analysis using `limma` between meaningful pairwise contrasts, followed by a gene set enrichment analysis (GSEA); B) summing up intensities on functional ontology terms (GO Biological Process) and then analyzing differences by ANOVAs and Tukey post hoc tests. The outputs of both analysis were combined by a semantic similarity analysis using `SimplifyEnrichment`. This shiny app is an interactive access to the results of these different analysis as well as the underlying data.

For more details, please refer to our manuscript: COMING SOON

The raw data is available on the PRoteomics IDEntification (PRIDE) database with the identifier [PXD035399](https://www.ebi.ac.uk/pride/archive/projects/PXD035399).

Our data analysis pipeline is available on [Zenodo](https://doi.org/10.5281/zenodo.6896565).
