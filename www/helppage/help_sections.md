The visualization for the analysis is split in two parts. The first part shows the overview over the functional analysis by GO semantic similarity. The second part focuses on visualization and statistics for individual GO terms.

#### Overview

In order to understand how the different GO terms showing up by the different analysis are connected, an overview is provided using GO semantic similarity. The results from the semantic similarity analysis are shown in two ways, as a heatmap with data from the ANOVA and GSEA analysis projected on it, or by focusing on the clusters of similar GO terms. 

The central part of the heatmap overview is the distance matrix between GO terms.
Additionally, data from both the ANOVA analysis and the GSEA is projected onto the heatmap. Since effects are mostly seen on the level of mutation status or condition, these are the main effects that are focused on. 
This data projection works the following way: for each visualized effect, a reference group is chosen. Then, all meaningful comparisons against this reference group are visualized. For the ANOVA analysis, this means all main effects against the reference group. For the GSEA analysis, this means that all contrasts are shown against the reference group where only the category of interest is varied. For example if the reference group is WT in the category mutation status, then all contrasts shown vary in mutation status, with one of them always being WT, whereas in all contrasts the condition and concentration are kept same. 
To reduce the complexity, the data heatmaps are kept simple and only show discrete output values, namely whether there is a significant effect with an adjusted p-value of lower than 0.05 and in which direction (meaning, whether it is significantly higher or significantly lower than the reference group) the effect is going. 
For the interactive version of the heatmap, the user is able to freely set the reference groups. For the static version of the heatmap, the reference groups are always WT for mutation status and unstimulated for condition and concentration. 
Which groups are compared can be seen from the color annotation on top of and below the data heatmaps which indicates the mutation status, condition and concentration for the reference group (bottom) and compared sample (top).

Secondly, the clusters obtained from performing binary cut clustering on the distance matrix can be explored in detail. This includes word clouds of the clusters, as well as a list of all containing GO terms. 

#### Individual GO Term

There are several ways to select a GO term of interest. It could come up from the **overview** analysis, or by a keyword search in the **ontology control panel**.
After a GO term is selected, there are information and visualization about this specific GO term to be explored. Firstly, name and description of the GO term of interest are displayed, together with how many proteins are part of this term and how many of them are found in the APMS data set. Additionally, semantically close GO term are displayed as well. 

Next, there are some visualizations that give an impression how the data coverage for this GO term is. **Samples per protein** are visualized, as well as the **number of identified proteins per sample**. 

On the next page, the **summed up LFQ intensities** for the groups of interest are shown, together with the **results of the statistical analysis**. 

Finally, there is the option of inspect the **LFQ intensities for individual proteins** in the ontology.
