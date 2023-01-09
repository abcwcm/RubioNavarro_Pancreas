- **GEO for single-cell RNA-seq**: [GSE203151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203151)
  - Relevant RDS files referenced in scripts can be found here:
    - [sce_integratedData_HFD_WT_2019-09.rds](https://wcm.box.com/s/woo4m5g39gtox94rh8hhx8ugv6wz39v0)
    - [sce_integratedData_betaCells_HFD_WT_2019-09.rds](https://wcm.box.com/s/ncr8szvwoeu1ip3uxtizwk3vsxzzczc6)

The main goal was to identify difference in the transcriptomes of pancreatic beta cells from pre-diabetic mice fed with regular diet (=RD) or high-fat diet (=HFD). To this end, the whole pancreata were digested and subjected to single-cell sequencing; beta cells were identified in silico based on marker gene expression patterns and clustering.

The raw reads were aligned and processed with the CellRanger pipeline (v. 2.1.0) using the **mouse transcriptome and genome version mm10**. Subsequent analyses were performed in R following the recommendations of Amezquita et al. (https://osca.bioconductor.org/) using numerous functions provided in the R packages scater and scran.
Based on the calculation of outliers, we removed cells with fewer than 1,000 UMI. In addition, cells with more than 7.5% mitochondrial reads were removed as well as genes that were expressed in fewer than 5 cells of the same sample type.

We then processed and integrated the different samples using Seurat version 3.1 following the recommendations from the Satija Lab’s vignette (https://satijalab.org/seurat/v3.1/integration.html).
More specifically, read counts were first normalized using `SCTransform` for each sample individually.
The different samples were then integrated using the top 3,000 most informative genes before performing various dimensionality reduction steps including PCA and UMAP.
A shared nearest neighbor graph was constructed using Seurat’s `FindNeighbors` function with default settings (e.g. k=20) using the first 20 principal components. Subsequent clustering was performed with Seurat’s `FindClusters` function.
In order to find major subclusters of  cells that are transcriptionally distinct and not rare (<1% of beta cells), we set the resolution parameter to 0.2.
For visualizations and assessments of normalized expression values, the SCTransform-normalized (log-transformed) expression values were used unless noted otherwise.

To remove putative doublets that may have resulted from the simultaneous capture of cells representing two different cell types within the same GEM, we additionally employed normalized expression thresholds for each of the marker genes and cells that expressed more than one marker gene above the given threshold were removed from the downstream analyses focusing on beta cells.

To determine genes that show significant expression differences that can be attributed to the condition (RD or HFD), we first extracted all the beta cells from both conditions and re-calculated the most variable genes, PCA and UMAP coordinates as well as clusters (resolution parameter = 0.2).

For each cluster, we then used the `findMarkers` function of the scran package, comparing cells of the RD condition to those of the HFD condition and performing gene-wise one-tailed t-tests testing for up-regulation.
Differentially expressed genes (DEG) were defined by an FDR threshold of 1%.
Gene enrichment analyses to investigate canonical pathways, biological processes or transcription factors were performed using the Ingenuity Pathway Analysis tool (Qiagen). 
