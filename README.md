# RubioNavarro: Mouse Pancreas

Scripts related to single-cell and bulk RNA-seq data of mouse pancreata as presented by Rubio-Navarro et al.

- **GEO for bulk RNA-seq**: [GSE205023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205023)
- **GEO for single-cell RNA-seq**: [GSE203151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203151)

## Single-cell RNA-seq data generation

Pancreatic islets were isolated and cultured overnight in RPMI 1640 supplemented with 10% FBS and penicillin/streptomycin. Then, islets were dissociated with 0.05% Trypsin (Corning) for 10 minutes at 37C. Dissociated cells were passed through a 35 um nylon mesh and resuspended in 1x PBS + 0.05% BSA. Cellular viability was higher than 95% as determined by using Trypan Blue assay and doublets represented less than 5% of total cells. Single cell suspensions were sequenced by the Chromium Single Cell 3’ Reagent Kit v3 (10x Genomics) and 10X Genomics’ Chromium Controller at the Weill Cornell Medicine Genomics Core Facility.

## Single-cell RNA-seq data analysis

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

## Bulk RNA-seq library generation and analysis

After RNA isolation, total RNA integrity was analyzed using an Agilent 2100 Bioanalyzer and concentrations were measured using the NanoDrop Spectrophotometer (Thermo Fisher Scientific). Preparation of the RNA sample library was performed by the Genomics Core Laboratory at Weill Cornell Medicine.
cDNA was generated using the SMARTer v4 Ultra Low Input RNA Kit (Clontech, # 63488) and Nextera XT kit (Nextera XT DNA Library Preparation Kit (Illumina, San Diego, CA) was used to make cDNA libraries suitable for Illumina sequencing.
The normalized cDNA libraries were pooled and sequenced on an Illumina HiSeq4000 sequencer at 50 pair-end cycles.
Transcriptome Data Analysis Sequencing reads were mapped with STAR v2.6.0c with default parameters to the **mouse reference genome 76 (GRCm38.p6)**.
Fragments per gene were counted with `featureCounts` v1.6.2 with respect to Ensembl annotations.
Differentially expressed genes between pairwise comparisons were identified by Wald tests using `DESeq2` v1.26.0 77 and only Benjamini–Hochberg corrected P-values < 0.05 were considered statistically significant.
Biological analyses, including canonical pathways, biological processes or transcription factors were performed using the Ingenuity Pathway Analysis tool (Qiagen). Base-2 log-transformed counts per million (CPM) values were used for heatmap plots of bulk RNA-seq data, which were centered and scaled by row. 

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)

Data analyses were performed by Paul Zumbo & Friederike Dündar at the [Applied Bioinformatics Core](https://abc.med.cornell.edu/) of Weill Cornell Medicine. 
Don't hesitate to get in touch, either via issues raised here or via email to abc at med.cornell.edu
