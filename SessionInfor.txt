R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.0

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SeuratWrappers_0.3.2        monocle3_1.3.7              DropletUtils_1.22.0         SoupX_1.6.2                
 [5] scater_1.30.1               SingleR_2.4.1               scRNAseq_2.16.0             scran_1.30.2               
 [9] scuttle_1.12.0              SingleCellExperiment_1.24.0 celldex_1.12.0              SummarizedExperiment_1.32.0
[13] Biobase_2.62.0              GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0             
[17] S4Vectors_0.40.2            BiocGenerics_0.48.1         MatrixGenerics_1.14.0       matrixStats_1.3.0          
[21] cellhashR_1.1.0             sctransform_0.4.1           dittoSeq_1.14.2             Seurat_5.0.3               
[25] SeuratObject_5.0.2          sp_2.1-4                    lubridate_1.9.3             forcats_1.0.0              
[29] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.5                
[33] tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1               tidyverse_2.0.0            

loaded via a namespace (and not attached):
  [1] ProtGenerics_1.34.0           spatstat.sparse_3.0-3         bitops_1.0-7                 
  [4] httr_1.4.7                    RColorBrewer_1.1-3            tools_4.3.1                  
  [7] utf8_1.2.4                    R6_2.5.1                      HDF5Array_1.30.1             
 [10] lazyeval_0.2.2                uwot_0.1.16                   rhdf5filters_1.14.1          
 [13] withr_3.0.0                   prettyunits_1.2.0             gridExtra_2.3                
 [16] progressr_0.14.0              cli_3.6.2                     spatstat.explore_3.2-7       
 [19] fastDummies_1.7.3             spatstat.data_3.0-4           ggridges_0.5.6               
 [22] pbapply_1.7-2                 Rsamtools_2.18.0              R.utils_2.12.3               
 [25] parallelly_1.37.1             limma_3.58.1                  rstudioapi_0.16.0            
 [28] RSQLite_2.3.6                 BiocIO_1.12.0                 generics_0.1.3               
 [31] ica_1.0-3                     spatstat.random_3.2-3         Matrix_1.6-5                 
 [34] ggbeeswarm_0.7.2              fansi_1.0.6                   abind_1.4-5                  
 [37] R.methodsS3_1.8.2             lifecycle_1.0.4               yaml_2.3.8                   
 [40] edgeR_4.0.16                  rhdf5_2.46.1                  SparseArray_1.2.4            
 [43] BiocFileCache_2.10.2          Rtsne_0.17                    grid_4.3.1                   
 [46] blob_1.2.4                    promises_1.3.0                dqrng_0.4.0                  
 [49] ExperimentHub_2.10.0          crayon_1.5.2                  miniUI_0.1.1.1               
 [52] lattice_0.22-6                beachmat_2.18.1               cowplot_1.1.3                
 [55] GenomicFeatures_1.54.4        KEGGREST_1.42.0               pillar_1.9.0                 
 [58] knitr_1.45                    metapod_1.10.1                boot_1.3-30                  
 [61] rjson_0.2.21                  future.apply_1.11.2           codetools_0.2-20             
 [64] leiden_0.4.3.1                glue_1.7.0                    remotes_2.5.0                
 [67] data.table_1.15.4             vctrs_0.6.5                   png_0.1-8                    
 [70] rmdformats_1.0.4              spam_2.10-0                   gtable_0.3.5                 
 [73] cachem_1.1.0                  xfun_0.44                     S4Arrays_1.2.1               
 [76] mime_0.12                     survival_3.6-4                pheatmap_1.0.12              
 [79] statmod_1.5.0                 bluster_1.12.0                interactiveDisplayBase_1.40.0
 [82] fitdistrplus_1.1-11           ROCR_1.0-11                   nlme_3.1-164                 
 [85] bit64_4.0.5                   progress_1.2.3                filelock_1.0.3               
 [88] RcppAnnoy_0.0.22              irlba_2.3.5.1                 vipor_0.4.7                  
 [91] KernSmooth_2.23-24            colorspace_2.1-0              DBI_1.2.2                    
 [94] tidyselect_1.2.1              bit_4.0.5                     compiler_4.3.1               
 [97] curl_5.2.1                    BiocNeighbors_1.20.2          xml2_1.3.6                   
[100] DelayedArray_0.28.0           plotly_4.10.4                 rtracklayer_1.62.0           
[103] bookdown_0.39                 scales_1.3.0                  lmtest_0.9-40                
[106] rappdirs_0.3.3                digest_0.6.35                 goftest_1.2-3                
[109] minqa_1.2.6                   spatstat.utils_3.0-4          rmarkdown_2.27               
[112] XVector_0.42.0                htmltools_0.5.8.1             pkgconfig_2.0.3              
[115] lme4_1.1-35.3                 sparseMatrixStats_1.14.0      ensembldb_2.26.0             
[118] dbplyr_2.5.0                  fastmap_1.2.0                 rlang_1.1.3                  
[121] htmlwidgets_1.6.4             shiny_1.8.1.1                 DelayedMatrixStats_1.24.0    
[124] zoo_1.8-12                    jsonlite_1.8.8                BiocParallel_1.36.0          
[127] R.oo_1.26.0                   BiocSingular_1.18.0           RCurl_1.98-1.14              
[130] magrittr_2.0.3                GenomeInfoDbData_1.2.11       dotCall64_1.1-1              
[133] patchwork_1.2.0               Rhdf5lib_1.24.2               munsell_0.5.1                
[136] Rcpp_1.0.12                   viridis_0.6.5                 reticulate_1.36.1            
[139] stringi_1.8.4                 zlibbioc_1.48.2               MASS_7.3-60.0.1              
[142] AnnotationHub_3.10.1          plyr_1.8.9                    parallel_4.3.1               
[145] listenv_0.9.1                 ggrepel_0.9.5                 deldir_2.0-4                 
[148] Biostrings_2.70.3             splines_4.3.1                 tensor_1.5                   
[151] hms_1.1.3                     locfit_1.5-9.9                igraph_2.0.3                 
[154] spatstat.geom_3.2-9           RcppHNSW_0.6.0                biomaRt_2.58.2               
[157] reshape2_1.4.4                ScaledMatrix_1.10.0           XML_3.99-0.16.1              
[160] BiocVersion_3.18.1            evaluate_0.23                 BiocManager_1.30.23          
[163] nloptr_2.0.3                  tzdb_0.4.0                    httpuv_1.6.15                
[166] RANN_2.6.1                    polyclip_1.10-6               future_1.33.2                
[169] scattermore_1.2               rsvd_1.0.5                    xtable_1.8-4                 
[172] AnnotationFilter_1.26.0       restfulr_0.0.15               RSpectra_0.16-1              
[175] later_1.3.2                   viridisLite_0.4.2             GenomicAlignments_1.38.2     
[178] memoise_2.0.1                 beeswarm_0.4.0                AnnotationDbi_1.64.1         
[181] cluster_2.1.6                 timechange_0.3.0              globals_0.16.3  
