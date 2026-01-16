### Introduction

This repository contains a pseudobulk differential expression analysis on a single-cell PBMC atlas from Gong et al., using data accessible via CellXGene and GEO (GSE275067). The code is organized in the code/ folder with two primary files: a Jupyter notebook for Python preprocessing for celltypist labeling and per-donor aggregation into pseudobulks, and a R Markdown  that performs limma-voom DE, visualization, and gene list export. The dataset comprises 235 donors, profiled with 10x Genomics Flex scRNA-seq on peripheral blood mononuclear cells. The broader study profiled peripheral immunity across adulthood using scRNA-seq, proteomics, and flow cytometry, and reported age-associated transcriptional remodeling in T cell subsets and shifts in immune composition with age. Here, I focused on a paired pseudobulk comparison between IL1B+ CD14 monocytes and GZMK- CD56dim NK cells. 

### Methods:

Data were downloaded from GEO (GSE275067). Single-cell count matrices were read with anndata v0.12.6 and processed in Python v3.13.9 using scanpy v1.11.1, celltypist v1.7.1, pandas v2.3.3, and numpy v2.3.5. For cell type labeling, counts were lightly filtered at the gene level, where genes expressed in at least 20 cells were retained, then normalized only for annotation using normalize_total to 1e4 and log1p, and annotated using the AIFI_L3 CellTypist model. After annotation, raw counts were restored for pseudobulk aggregation. Pseudobulks were constructed per donor and cell type by summing raw counts across cells of the same annotated type; I had at least 20 cells per pseudobulk to include a donor×cell type in downstream analyses.

Model resources and labeling steps were guided by the authors’ public documentation: https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/models/

Aggregated pseudobulk counts and metadata were exported to CSV and analyzed in R v4.4.3. Differential expression and visualization used edgeR v4.4.2, limma v3.62.2, ggplot2 v4.0.1, and dplyr v1.1.4. For differential expression, I focused on IL1B+ CD14 monocytes versus GZMK- CD56dim NK cells, and these cell types were selected because they had some of the highest cell numbers across donors. I then applied a paired limma-voom design (~ sample_id + cell_type) with TMM normalization, filterByExpr to remove lowly expressed genes, voom transformation, and empirical Bayes moderation. Donors lacking one of the two target cell types were excluded to ensure a paired comparison. This cell type–level pseudobulk approach was chosen because sequencing runs were pooled, and donor level aggregation provides a robust framework for between cell type contrasts while controlling for inter-donor variability.

### Results

Across all donors, the annotated cell counts included 552,613 GZMK- CD56dim NK cells and 443,051 IL1B+ CD14 monocytes. The volcano plot (Figure 1) shows a broad separation of genes with significant differences between the two cell types, with most genes passing an FDR threshold (padj < 0.05) and exhibiting nonzero log2 fold changes, and comparatively few nonsignificant genes. A PCA using voom log-CPM further supports this separation: the top 30 most significant genes cleanly cluster pseudobulks by cell type (Figure 2).

### Figure 1, Volcano plot comparing IL1B+ CD14 monocytes versus GZMK- CD56dim NK cells

!(Volcano plot)[images/VolcanoPlot.png]

### Figure 2, PCA of pseudobulks using the top 30 most significant genes

!(PCA)[images/PCA.png]

The strong differential signal is consistent with the distinct biological roles of these cell types and the large numbers of cells available for aggregation. IL1B+ CD14 monocytes are a pro-inflammatory subset of CD14 monocytes implicated in arthritic pathology (see Allen AIFI cell type description: https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/cell-type-descriptions/monocytes/). GZMK- CD56dim NK cells represent CD56dim NK populations that did not meet marker thresholds for other NK subtypes (see NK cell description: https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/cell-type-descriptions/nk-cells-and-ilcs/). 

After identifying differentially expressed genes between IL1B+ CD14 monocytes and GZMK- CD56dim NK cells, I performed pathway enrichment using Enrichr against Reactome 2024, separately on up- and down-regulated gene lists. The top upregulated pathways included Neutrophil Degranulation, Immune System, and Innate Immune System, reflecting the inflammatory and innate immune functions characteristic of monocytes. The top downregulated pathways included Gene Expression (Transcription), RNA Polymerase II Transcription, and the Generic Transcription Pathway, alongside cell cycle and DNA repair processes, consistent with reduced proliferative/transcriptional programs relative to NK cell states in this contrast. These enrichment results align with the biological roles of the two cell types and reinforce the strong separation observed in the volcano plot and PCA.

### Table 1, Upregulated pathways

| Index | Name                                                                             | P-value   | Adjusted p-value | Odds Ratio | Combined score |
|-------|----------------------------------------------------------------------------------|-----------|------------------|------------|----------------|
| 1     | Neutrophil Degranulation                                                         | 7.816e-59 | 1.566e-55        | 4.59       | 614.52         |
| 2     | Immune System                                                                    | 7.297e-49 | 7.308e-46        | 2.01       | 222.95         |
| 3     | Innate Immune System                                                             | 1.961e-43 | 1.310e-40        | 2.37       | 232.74         |
| 4     | Cytokine Signaling in Immune System                                              | 1.663e-24 | 8.328e-22        | 2.15       | 117.97         |
| 5     | Membrane Trafficking                                                             | 2.702e-23 | 1.083e-20        | 2.27       | 118.11         |
| 6     | Signal Transduction                                                              | 1.701e-20 | 5.680e-18        | 1.51       | 68.78          |
| 7     | Signaling by Interleukins                                                        | 3.808e-19 | 1.090e-16        | 2.38       | 100.87         |
| 8     | Metabolism                                                                       | 4.313e-18 | 1.080e-15        | 1.52       | 60.59          |
| 9     | Diseases of Signal Transduction by Growth Factor Receptors and Second Messengers | 8.502e-18 | 1.892e-15        | 2.30       | 90.52          |
| 10    | Signaling by Rho GTPases, Miro GTPases and RHOBTB3                               | 2.998e-16 | 6.005e-14        | 1.92       | 68.68          |


### Table 2, Downregulated pathways

| Index | Name                                                | P-value   | Adjusted p-value | Odds Ratio | Combined score |
|-------|-----------------------------------------------------|-----------|------------------|------------|----------------|
| 1     | Gene Expression (Transcription)                     | 1.816e-41 | 3.602e-38        | 2.05       | 192.58         |
| 2     | RNA Polymerase II Transcription                     | 3.718e-37 | 3.688e-34        | 2.08       | 174.61         |
| 3     | Generic Transcription Pathway                       | 1.044e-36 | 6.901e-34        | 2.14       | 177.48         |
| 4     | Cell Cycle                                          | 3.505e-31 | 1.738e-28        | 2.56       | 179.36         |
| 5     | Cell Cycle, Mitotic                                 | 1.056e-21 | 4.189e-19        | 2.38       | 114.71         |
| 6     | DNA Repair                                          | 1.898e-19 | 6.277e-17        | 2.82       | 121.75         |
| 7     | Activation of the Pre-Replicative Complex           | 4.447e-13 | 1.260e-10        | 23.27      | 661.91         |
| 8     | Epigenetic Regulation of Gene Expression            | 6.638e-12 | 1.646e-9         | 2.27       | 58.37          |
| 9     | Activation of ATR in Response to Replication Stress | 2.030e-11 | 4.475e-9         | 12.02      | 296.03         |
| 10    | S Phase                                             | 3.771e-11 | 7.482e-9         | 2.98       | 71.41          |

The top most upregulated pathways were neutrophil degranulation, immune system, and innate immune system. The top most downregulated pathways were gene expression (transcription),  RNA Polymerase II transcription, and generic transcription pathway. These pathways are consistent with the given cell types. 

### Conclusion:
This pseudobulk analysis, built on the Allen Institute for Immunology’s single-cell atlas and focused on the follow-up cohort, demonstrates transcriptional differences between IL1B+ CD14 monocytes and GZMK- CD56dim NK cells. By aggregating counts at the donor×cell-type level and using a paired design, I obtained clear DE signals and pathway-level insights that are coherent with known immunobiology. Future directions include expanding pseudobulk DE to additional cell types to map similarity/dissimilarity across the immune landscape, and incorporating age or clinical covariates into the design to probe interactions.

Works cited:
Gong, Q., Sharma, M., Glass, M.C. et al. Multi-omic profiling reveals age-related immune dynamics in healthy adults. Nature 648, 696–706 (2025). https://doi.org/10.1038/s41586-025-09686-5
