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

## Software versions

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Catalina 10.15.4 ##
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib ##
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8 ##
## attached base packages:
## [1] stats4 parallel stats
## [8] methods   base
##
## other attached packages:
##  [1] ReactomePA_1.30.0           clusterProfiler_3.14.3     
##  [3] destiny_3.0.1               dittoSeq_1.1.2             
##  [5] slingshot_1.4.0             princurve_2.1.4            
##  [7] patchwork_1.0.0             scater_1.14.6              
##  [9] ggplot2_3.3.0               SingleCellExperiment_1.8.0 
## [11] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
## [13] BiocParallel_1.20.1         matrixStats_0.56.0         
## [15] Biobase_2.46.0              GenomicRanges_1.38.0       
## [17] GenomeInfoDb_1.22.1         IRanges_2.20.2             
## [19] S4Vectors_0.24.4            BiocGenerics_0.32.0        
## [21] data.table_1.12.8           magrittr_1.5               
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.1           ggthemes_4.2.0           tidyr_1.0.2             
##   [4] bit64_0.9-7              knitr_1.28               irlba_2.3.3             
##   [7] RCurl_1.98-1.1           cowplot_1.0.0            RSQLite_2.2.0           
##  [10] RANN_2.6.1               europepmc_0.3            proxy_0.4-24            
##  [13] future_1.16.0            bit_1.1-15.2             enrichplot_1.6.1        
##  [16] webshot_0.5.2            xml2_1.3.1               assertthat_0.2.1        
##  [19] viridis_0.5.1            xfun_0.13                hms_0.5.3               
##  [22] evaluate_0.14            DEoptimR_1.0-8           fansi_0.4.1             
##  [25] progress_1.2.2           caTools_1.18.0           readxl_1.3.1            
##  [28] htmlwidgets_1.5.1        igraph_1.2.5             DBI_1.1.0               
##  [31] purrr_0.3.3              ellipsis_0.3.0           RSpectra_0.16-0         
##  [34] dplyr_0.8.5              backports_1.1.6          vctrs_0.2.4             
##  [37] TTR_0.23-6               ROCR_1.0-7               abind_1.4-5             
##  [40] RcppEigen_0.3.3.7.0      withr_2.1.2              ggforce_0.3.1           
##  [43] triebeard_0.3.0          robustbase_0.93-6        checkmate_2.0.0         
##  [46] vcd_1.4-7                sctransform_0.2.1        scran_1.14.6            
##  [49] xts_0.12-0               prettyunits_1.1.1        cluster_2.1.0           
##  [52] DOSE_3.12.0              lazyeval_0.2.2           ape_5.3                 
##  [55] laeken_0.5.1             crayon_1.3.4             labeling_0.3            
##  [58] edgeR_3.28.1             pkgconfig_2.0.3          tweenr_1.0.1            
##  [61] nlme_3.1-147             vipor_0.4.5              drat_0.1.5              
##  [64] nnet_7.3-13              rlang_0.4.5              globals_0.12.5          
##  [67] lifecycle_0.2.0          rsvd_1.0.3               cellranger_1.1.0        
##  [70] polyclip_1.10-0          RcppHNSW_0.2.0           lmtest_0.9-37           
##  [73] graph_1.64.0             Matrix_1.2-18            urltools_1.7.3          
##  [76] carData_3.0-4            boot_1.3-24              zoo_1.8-7               
##  [79] beeswarm_0.2.3           ggridges_0.5.2           pheatmap_1.0.12         
##  [82] png_0.1-7                viridisLite_0.3.0        bitops_1.0-6            
##  [85] KernSmooth_2.23-16       blob_1.2.1               DelayedMatrixStats_1.8.0
##  [88] stringr_1.4.0            qvalue_2.18.0            readr_1.3.1             
##  [91] gridGraphics_0.5-0       reactome.db_1.70.0       scales_1.1.0            
##  [94] memoise_1.1.0            graphite_1.32.0          plyr_1.8.6              
##  [97] hexbin_1.28.1            ica_1.0-2                gplots_3.0.3            
## [100] gdata_2.18.0             zlibbioc_1.32.0          compiler_3.6.2          
## [103] lsei_1.2-0               ABCutilities_0.3.3       dqrng_0.2.1             
## [106] kableExtra_1.1.0         RColorBrewer_1.1-2       pcaMethods_1.78.0       
## [109] fitdistrplus_1.0-14      cli_2.0.2                XVector_0.26.0          
## [112] listenv_0.8.0            pbapply_1.4-2            ggplot.multistats_1.0.0 
## [115] MASS_7.3-51.5            tidyselect_1.0.0         stringi_1.4.6           
## [118] forcats_0.5.0            highr_0.8                yaml_2.2.1              
## [121] GOSemSim_2.12.1          BiocSingular_1.2.2       locfit_1.5-9.4          
## [124] ggrepel_0.8.2            grid_3.6.2               fastmatch_1.1-0         
## [127] tools_3.6.2              future.apply_1.5.0       rio_0.5.16              
## [130] rstudioapi_0.11          foreign_0.8-76           gridExtra_2.3           
## [133] smoother_1.1             scatterplot3d_0.3-41     farver_2.0.3            
## [136] Rtsne_0.15               ggraph_2.0.2             digest_0.6.25           
## [139] rvcheck_0.1.8            BiocManager_1.30.10      Rcpp_1.0.4              
## [142] car_3.0-8                RcppAnnoy_0.0.16         httr_1.4.1              
## [145] AnnotationDbi_1.48.0     npsurv_0.4-0             scABC2_0.3.1            
## [148] colorspace_1.4-1         reticulate_1.15          rvest_0.3.5             
## [151] ranger_0.12.1            splines_3.6.2            statmod_1.4.34          
## [154] uwot_0.1.8               graphlayouts_0.6.0       sp_1.4-1                
## [157] ggplotify_0.0.5          plotly_4.9.2.1           jsonlite_1.6.1          
## [160] tidygraph_1.1.2          R6_2.4.1                 pillar_1.4.3            
## [163] htmltools_0.4.0          glue_1.4.0               VIM_5.1.1               
## [166] BiocNeighbors_1.4.2      class_7.3-16             codetools_0.2-16        
## [169] fgsea_1.12.0             tsne_0.1-3               lattice_0.20-41         
## [172] tibble_3.0.0             curl_4.3                 ggbeeswarm_0.6.0        
## [175] leiden_0.3.3             gtools_3.8.2             zip_2.0.4               
## [178] GO.db_3.10.0             openxlsx_4.1.4           survival_3.1-12         
## [181] limma_3.42.2             rmarkdown_2.1            munsell_0.5.0           
## [184] e1071_1.7-3              DO.db_2.9                GenomeInfoDbData_1.2.2  
## [187] haven_2.2.0              reshape2_1.4.4           gtable_0.3.0            
## [190] Seurat_3.1.5
```
