# HeLa-SVFmix-scseq
This script is for the test of HeLA-SVF mix cells
#This script is copied from SVF_myscseq/src/1-preprocessing.R 
#But the feature barcode matrix file is from soupx cleaned matrix 
#The soupx cleaning process is run on the original filtered barcode matrix (by cellranger) along with its umap results 

#save.image("results/SoupX_Smith_SVFmix.RData")
####Get RNA assay and HTO assay as separate Seurat object####
getwd()
# [1] "/Users/yujunfeng/SVF_myscseq/SoupX_SVFmix_process"

library(tidyverse)
library(Seurat)
library(scater) 
library(sctransform) 
library(ggplot2) 
library(dittoSeq) 
library(SoupX)
?Read10X

# Load data that is cleaned by soupx 
scData<-Read10X(data.dir="soupX_SVFmix_filt/")
#133 MB

dim(scData)
# [1] 41954  2051
head(colnames(scData))

?CreateSeuratObject
# Create Seurat object 
SVFmix<-CreateSeuratObject(
  counts=scData,
  names.field=1,
  names.delim = "-",
  min.cells = 10)
# min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. 
dim(SVFmix)
# if set min.cells=10 [1] 27295  2051  
# if set min.cells=20 [1] 25547  2051 
# I keep using 10, since there could be some minor immune cell types 

# about 15K genes were filtered out after setting min.cell=20) 
DefaultAssay(SVFmix)
# [1] "RNA" 

SVFmix[["nCount_RNA"]] %>% head()
levels(SVFmix$orig.ident)
levels(SVFmix_barcode$orig.ident)
#barcode assay is not attached here 
rm(scData)
gc() #to reclaim memory 

#### Check library QC####
# look at total RNA count and Feature RNA count for RNA and HTO assay 
#beacuse I used min.cells=10, both nCount and nFeature increased a little bit 

summary(SVFmix[["nCount_RNA"]])
#nCount_RNA    
#Min.   :   634  
#1st Qu.:  2790  
#Median : 15443  
#Mean   : 58291  
#3rd Qu.: 90528  
#Max.   :567514  
summary(SVFmix[["nFeature_RNA"]])
#nFeature_RNA  
#Min.   :  569  
#31st Qu.: 1708  
#3Median : 4231  
#Mean   : 5159  
#3rd Qu.: 8780  
#Max.   :18617 

#### Calculate MT percentage####

#mitochondria genes start with mt- that is shown in row names 
#^ means "start with" 
MT.genes.human <- grep("^GRCh38-MT-", rownames(SVFmix), value = TRUE)
length(MT.genes.human)
#"13"
MT.genes.mouse <- grep("^GRCm39-mt-", rownames(SVFmix), value = TRUE)
length(MT.genes.mouse)
#"13"
MT.genes<-c(MT.genes.human, MT.genes.mouse)

#calculate percentage of UMIs in MT.genes
?PercentageFeatureSet
SVFmix[["percent.mt"]] <- PercentageFeatureSet(SVFmix, 
                                               features = MT.genes)
head(SVFmix[[]])
# this is using the matrix data without running soupX 
#orig.ident nCount_RNA nFeature_RNA percent.mt
#AATCACGAGATCGCCC-1 SeuratProject      45882         7429   2.508609
#AATCACGAGCAATTAG-1 SeuratProject      28804         5071   3.232190
#AATCACGAGCATTGAA-1 SeuratProject      19631         4136   1.151240
#AATCACGAGCGAATGC-1 SeuratProject      42468         6896   2.352359

# this is using the matrix data after running soupX 

#orig.ident nCount_RNA nFeature_RNA percent.mt
#AATCACGAGATCGCCC-1 SeuratProject      45420         7193   2.514311
#AATCACGAGCAATTAG-1 SeuratProject      28557         4920   3.246139
#AATCACGAGCATTGAA-1 SeuratProject      19553         4115   1.145604
#AATCACGAGCGAATGC-1 SeuratProject      42198         6784   2.338973
#not change much in percent.mt 

summary(SVFmix[["percent.mt"]])
#before soupX: percent.mt      
#Min.   : 0.05079  
#1st Qu.: 2.16764  
#Median : 2.83303  
#Mean   : 4.23114  
#3rd Qu.: 3.91613  
#Max.   :78.63085

#percent.mt     
#Min.   :0.03005  
#1st Qu.:2.05421  
#Median :2.63632  
#Mean   :2.71285  
#3rd Qu.:3.34630  
#Max.   :5.02344 

#open XQuartz
x11(width = 12, height = 6)
levels(SVFmix$orig.ident)<-"SVF_HeLa_Mix"
head(SVFmix[[]])

#image on non-log scale
VlnPlot(SVFmix, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.5)



# image on log scale 
x11(width = 12, height = 6)
p<-VlnPlot(SVFmix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.5, log = TRUE)
ggsave("~/SVF_myscseq/SoupX_SVFmix_process/results/myplot.jpeg", plot = p, width = 12, height = 6, dpi = 300)

plot1 <- FeatureScatter(SVFmix, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt") 
plot2 <- FeatureScatter(SVFmix, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(SVFmix, 
                        feature1 = "percent.mt", 
                        feature2 = "nFeature_RNA")

x11()
x11(width = 15, height = 5)
p<-plot1 + plot2 + plot3 
ggsave("results/QC_preMT_Scatter.jpeg",plot=p ,width = 15, height = 5, dpi=300 ) 

#### Cell Filtering#### 

## MT threshold##
#based on MT% 
#here used a 3 mad method

BiocManager::install("scater")
browseVignettes("scater")
library(scater)

temp1 <- scater::isOutlier(SVFmix$percent.mt, 
                           nmads = 3, type = "higher")
table(temp1)

#temp1
#FALSE 
#2051 
#because the soupx used filtered cell list, so they all passed the MT filtering this time


## 
#calculate number of human and mouse gene for each cell

gene.human.temp<-grep("^GRCh",rownames(SVFmix),value=TRUE)
gene.mouse.temp<-grep("^GRCm",rownames(SVFmix),value=TRUE)

SVFmix[["UMI_human"]] <- PercentageFeatureSet(SVFmix, 
                                              features = gene.human.temp)*SVFmix$nCount_RNA/100

SVFmix[["UMI_mouse"]] <- PercentageFeatureSet(SVFmix, 
                                              features = gene.mouse.temp)*SVFmix$nCount_RNA/100
plot1 <- FeatureScatter(SVFmix, 
                        feature1 = "UMI_human", 
                        feature2 = "UMI_mouse") +
  scale_x_log10() + scale_y_log10()
ggsave("results/UMI_before_filter.jpeg",plot=plot1 ,width = 10, height = 10, dpi=300 ) 

plot2 <- FeatureScatter(SVFmix, 
                        feature1 = "UMI_human", 
                        feature2 = "UMI_mouse")  +
  scale_x_log10() + scale_y_log10() +
  geom_hline(yintercept=1000) +
  geom_vline(xintercept=5000) 
plot2 
ggsave("results/cutoff-2.jpeg",plot=plot2 ,width = 10, height = 10, dpi=300 )  

SVFmix_filter<-subset(SVFmix, subset =UMI_human > 5000 |UMI_mouse > 2000) 
dim(SVFmix_filter) 
#[1] 27295  1487  

####Normalization#### 
library(sctransform)
library(ggplot2)
DefaultAssay(SVFmix)
#"RNA"
?SCTransform
SVFmix_filter <- SCTransform(SVFmix_filter, 
                             method = "glmGamPoi",
                             vars.to.regress = "percent.mt", 
                             return.only.var.genes = FALSE, 
                             vst.flavor = "v2")
#now default assay is set to SCT 

DefaultAssay(SVFmix_filter)
#"SCT" 

#### Run PCA####
SVFmix_filter <- RunPCA(SVFmix_filter, verbose = FALSE)
# Use an ElbotPlot to see how many 
# PCs explain non-negligible
# variation in the 3000 genes
x11()
p<-ElbowPlot(SVFmix_filter, 50)
ggsave("results/2024-Feb-rerun/elbowplot.jpeg",plot=p ,width = 10, height = 10, dpi=300 )  
#looks like it plateaued at 40
ggsave("results/PNGs/elbowplot.jpeg") 

#### Run UMAP and tSNE#### 

# specify the same number of PC dims:

SVFmix_filter <- RunUMAP(SVFmix_filter, 
                         dims = 1:50, 
                         verbose = FALSE)

#### Visualize UMAP and tSNE plots ####

x11()
DimPlot(SVFmix_filter, label = FALSE, 
        reduction = "umap") 
ggsave("results/2024-Feb-rerun/UMAP.jpeg")

#### Use nearest neighbor graphs to call clusters ####

SVFmix_filter <- FindNeighbors(SVFmix_filter, dims = 1:50, verbose = FALSE)
# first start with defaults:
SVFmix_filter <- FindClusters(SVFmix_filter, verbose = TRUE)

head(SVFmix_filter[[]])
#                     orig.ident nCount_RNA nFeature_RNA percent.mt UMI_human UMI_mouse nCount_SCT nFeature_SCT
# AATCACGAGATCGCCC-1 SVF_HeLa_Mix      45420         7193   2.514311       937     44483      28919         6918
# AATCACGAGCAATTAG-1 SVF_HeLa_Mix      28557         4920   3.246139      1010     27547      27468         4814
# AATCACGAGCATTGAA-1 SVF_HeLa_Mix      19553         4115   1.145604      1332     18221      25638         4056
# AATCACGAGCGAATGC-1 SVF_HeLa_Mix      42198         6784   2.338973      1055     41143      28611         6535
# AATCACGAGTCCCAAT-1 SVF_HeLa_Mix      20255         4593   1.342878      1597     18658      25536         4526
# AATCACGAGTGGCCTC-1 SVF_HeLa_Mix      15166         4641   2.499011      1186     13980      24200         4551
#                    SCT_snn_res.0.8 seurat_clusters
# AATCACGAGATCGCCC-1               2               2
# AATCACGAGCAATTAG-1               1               1
# AATCACGAGCATTGAA-1               6               6
# AATCACGAGCGAATGC-1               0               0
# AATCACGAGTCCCAAT-1               2               2
# AATCACGAGTGGCCTC-1               1               1

table(SVFmix_filter$seurat_clusters)
# 0   1   2   3   4   5   6   7   8   9  10 
# 357 205 153 139 132 129 124 118  81  33  16 




