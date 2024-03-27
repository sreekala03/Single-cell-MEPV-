####Clustering of ME-PV cell types#########

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Nebulosa)
library(harmony)
library(readr)
library(RColorBrewer)
dir.create("output")

######LOAD SAMPLES & CREATE SEURAT OBJECT#######

sample_list<-lapply(paste0("s",1:9),function(s){
  sample<-CreateSeuratObject(Read10X(file.path("/data/Raw",s)),project = "MBH")
  sample[["sample"]]<-s
  sample<-PercentageFeatureSet(sample,pattern = "^mt-",col.name = "percent.mt")
  sample<-subset(sample, subset = nCount_RNA < 20000 & nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 15)
  counts<- GetAssayData(sample, assay="RNA")
  Rpl.genes <- rownames(counts) %>% stringr::str_subset(string = ., pattern = "^Rp[sl]")
  counts<-counts[-(which(rownames(counts) %in% c('Ehd2', 'Espl1', 'Jarid1d', 'Pnpla4',  'Rps4y1', 'Xist', 'Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d', 'Fos', 'Fosb', 'Gstp1', 'Egr1', 'Jun', 'Junb', 'Jund','Erh', 'Slc25a5', 'Pgk1', 'Eno1', 'Npas4', 'Tubb2a', 'Emc4', 'Scg5', Rpl.genes, 'Gm42418'))),]
  sample<-subset(sample, features=rownames(counts))
  return(sample)
})

###SCT transform on individual samples####
sample_list<-lapply(sample_list, SCTransform, vars.to.regress = c('percent.mt'),  return.only.var.genes=F)
saveRDS(sample_list, "output/sample_list.rds")

#Select Integration features from each sample
var_features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000, fvf.nfeatures = 3000)

#Merge samples
MEPV_cells <-Reduce(merge,sample_list)
saveRDS(MEPV_cells, "output/MEPV_cells.rds")

####Run PCA and UMAP####

MEPV_cells <- RunPCA(MEPV_cells,features = var_features)
ElbowPlot(MEPV_cells)
#Check genes between 47 and 50 dims
DimHeatmap(MEPV_cells, dims = 47:50, cells = 500, balanced = T)
MEPV_cells <- RunUMAP(object = MEPV_cells, dims = 1:50)
DimPlot(MEPV_cells,group.by="sample")+ggtitle("merged")


######## Adding metadata and reordering samples #########

MEPV_cells[["Sex"]]<-ifelse(MEPV_cells$sample%in%c("s1","s2","s4"),"Male","Female")
MEPV_cells[["Diet"]]<-sapply(MEPV_cells$sample,function(s)ifelse(s%in%c("s1","s3","s7","s8"),"SD",
                                                                 ifelse(s%in%c("s4","s9"),"HFDR",
                                                                        ifelse(s%in%c("s5"),"HFDS+D","HFDS"))))

MEPV_cells[["Sample.name"]]<-sapply(MEPV_cells$sample,function(s)ifelse(s%in%c("s1"),"Male.SD",
                                                                        ifelse(s%in%c("s2"),"Male.HFDS",
                                                                               ifelse(s%in%c("s3"),"Fem.SD.Dioest",
                                                                                      ifelse(s%in%c("s4"),"Male.HFDR",
                                                                                             ifelse(s%in%c("s5"),"Fem.HFDS+D",
                                                                                                    ifelse(s%in%c("s6"),"Fem.HFDS",
                                                                                                           ifelse(s%in%c("s7"),"Fem.SD.Proest",
                                                                                                                  ifelse(s%in%c("s8"),"Fem.SD.Oest",
                                                                                                                         ifelse(s%in%c("s9"),"Fem.HFDR"))))))))))
MEPV_cells[["Sample.name_2"]]<-sapply(MEPV_cells$sample,function(s)ifelse(s%in%c("s1"),"Male.SD",
                                                                          ifelse(s%in%c("s2"),"Male.HFDS",
                                                                                 ifelse(s%in%c("s3"),"Fem.SD.Dioest",
                                                                                        ifelse(s%in%c("s4"),"Male.HFDR",
                                                                                               ifelse(s%in%c("s5", "s6"),"Fem.HFDS",
                                                                                                      ifelse(s%in%c("s7"),"Fem.SD.Proest",
                                                                                                             ifelse(s%in%c("s8"),"Fem.SD.Oest",
                                                                                                                    ifelse(s%in%c("s9"),"Fem.HFDR")))))))))


######### INTERGRATION USING HARMONY ################

library(harmony)
MEPV_cells<-RunHarmony(samples,group.by.vars = "Sample.name",assay.use = "SCT", plot_convergence=TRUE)
# Harmony converged after 9 iterations
MEPV_cells <- RunUMAP(object = samples, dims = 1:50, reduction = "harmony", reduction.name = "humap",reduction.key = "hUMAP_", metric = 'euclidean')
p1 <- DimPlot(samples,reduction="humap",group.by="sample")+ggtitle("Harmony integrated 1:50")
p1

########## CLUSTERING ###########

MEPV_cells <- FindNeighbors(object = MEPV_cells, dims = 1:50,reduction = "harmony") 
MEPV_cells <- FindClusters(object = MEPV_cells, resolution = 1) 
MEPV_cells <- RunUMAP(object = MEPV_cells, dims = 1:50,reduction = "harmony", reduction.name = "humap",reduction.key = "hUMAP_", metric="euclidean")

DimPlot(MEPV_cells, reduction = "humap", label=TRUE)+ggtitle("Harmony Integrated dataset, res=1, dim 1:50")

######Cluster Markers 

Idents(MEPV_cells) <- "seurat_clusters"

MEPV_cells <- PrepSCTFindMarkers(MEPV_cells, assay = "SCT", verbose = TRUE)
Markers_int<-FindAllMarkers(MEPV_cells, min.pct=0.25, only.pos = TRUE, logfc.threshold = 0.4)

############# CELL TYPE ANNOTATION ###############

Idents(object = MEPV_cells) <- "seurat_clusters"

#Level 1 cell type annotation 

Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(4, 21, 29)))       <- "VLMC"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(11)))              <- "Plvap Endothelial"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(1)))               <- "Endothelial cells"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(14, 24)))          <- "Pericytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(23)))              <- "VSMC"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(22, 31)))          <- "Ccl5+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(10, 30)))          <- "Microglia" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(17,32)))           <- "CAMs" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(15)))              <- "Progenitors"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(20)))              <- "Differentiating"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(3,5, 33, 34, 26))) <- "Mature"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(27)))              <- "Lhb.Npy.Rax+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(7,2,9,13)))        <- "Tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(35)))              <- "Npy+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(8,16)))            <- "Cell membrane projections"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(19)))              <- "Avp+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(25)))              <- "Oxt+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(0,12)))            <- "Astrocytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(18, 28)))          <- "Pars.Tuberalis"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(6)))               <- "Ependymocytes"

MEPV_cells$Level1 <- Idents(object = MEPV_cells)

plot1 = DimPlot(MEPV_cells, reduction = "humap") & theme(legend.text = element_text(size = 6)) & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot1, id = "ident", size = 4, repel = T)

#No of cells in each cluster each condition
Idents (object = MEPV_cells) <-"Level1"
cell_num_v1 = table(Idents(MEPV_cells), MEPV_cells$Sample.name)


###Level 2 cell type annotation 

Idents(object = MEPV_cells) <- "seurat_clusters"

Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(4)))               <- "VLMC.1"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(21)))              <- "VLMC.2"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(29)))              <- "Dural fibroblasts"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(11)))              <- "Plvap Endothelial"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(1)))               <- "Endothelial cells"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(14, 24)))          <- "Pericytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(23)))              <- "VSMC"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(22, 31)))          <- "Immune cells"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(10, 30)))          <- "Microglia" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(17,32)))           <- "CAMs" 
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(15)))              <- "Progenitors"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(20)))              <- "Differentiating"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(33,26,3,5,34)))     <- "Mature"
#Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(27)))              <- "Lhb.Npy.Rax+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(7)))               <- "DMH tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(2)))               <- "VMH/dmARH tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(9)))               <- "vmARH tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(13)))              <- "ME tanycytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(35)))              <- "Npy+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(8,16)))            <- "Cell membrane projections"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(19)))              <- "Avp+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(25)))              <- "Oxt+"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(0,12)))            <- "Astrocytes"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(18)))              <- "Pars.Tuberalis"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(28)))              <- "Tight junction epithelial"
Idents(object = MEPV_cells, cells = WhichCells(MEPV_cells, ident = c(6)))               <- "Ependymocytes"

MEPV_cells$Level2 <- Idents(object = MEPV_cells)
plot1 = DimPlot(MEPV_cells, reduction = "humap") & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot1, id = "ident", size = 4, repel = F) 

################Remove Lhb.Npy.Rax+ cells which is be a dissection artefact################

Idents(object = MEPV_cells) <- "Level2"
MEPV_cells <-subset(MEPV_cells, idents=c("Lhb.Npy.Rax+"), invert = T)
MEPV_cells
# An object of class Seurat 
# 49838 features across 46922 samples within 2 assays 
# Active assay: SCT (18914 features, 0 variable features)
# 1 other assay present: RNA
# 4 dimensional reductions calculated: pca, umap, harmony, humap

#No of cells in each cluster each condition
Idents (object = MEPV_cells) <-"Level2"
cell_num_v2 = table(Idents(MEPV_cells), MEPV_cells$Sample.name)

##Rearrange Idents order 

MEPV_cells$Level2 <- factor(MEPV_cells$Level2,levels=c("Astrocytes", "Ependymocytes", "DMH tanycytes", "VMH/dmARH tanycytes", "vmARH tanycytes", 
                                                       "ME tanycytes", "Npy+", "Cell membrane projections", "Avp+", "Oxt+", 
                                                       "Microglia", "CAMs", "Immune cells", "Endothelial cells","Plvap Endothelial","Pericytes","VSMC",
                                                       "VLMC.1", "VLMC.2", "Dural fibroblasts", "Pars.Tuberalis", "Tight junction epithelial",
                                                       "Progenitors", "Differentiating", "Mature"))

###Assign colors to clusters####

col.pal <- list()
col.pal$celltype <- c("Astrocytes"="#A6CEE3", "DMH tanycytes"="#1F78B4", "VMH/dmARH tanycytes"="#B2DF8A", 
                      "vmARH tanycytes"="#33A02C", "ME tanycytes"="red",
                      "Npy+"="deeppink","Oxt+"="#FF7F00","Cell membrane projections"= "#CAB2D6", 
                      "Avp+"="#6A3D9A", "Immune cells"="#D95F02", "VLMC.1"="#abc4ff", 
                      "VLMC.2"="#B15928","Dural fibroblasts"="#CBD52E","Microglia"="#7570B3", 
                      "CAMs"="#E7298A", "Pars.Tuberalis"="#66A61E", "Tight junction epithelial"="#E6AB02",
                      "Plvap Endothelial"="#A6761D","Endothelial cells"="lightsalmon", "Pericytes"="yellow2","VSMC"="tomato1",
                      "Mature"="#0466c8", "Progenitors"="#38b000", "Differentiating"="ivory4", "Ependymocytes"="darkgoldenrod1")




pdf("output/Fig1c_Dimplot_MEPV.pdf", width = 12.5, height = 8)
DimPlot(MEPV_cells, group.by = "seurat_clusters", cols = col.pal$celltype, reduction = "humap" , label.box = T) & NoLegend() & NoAxes()
dev.off()

###Plots markers 

Idents(MEPV_cells) <- "Level2"
Markers_ano <-FindAllMarkers(MEPV_cells, min.pct=0.25, only.pos = TRUE, logfc.threshold = 0.6)

write.csv(Markers_ano, "MEPV/output/TableS1_MEPV_markers.csv")

Markers_ano %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
top5 <- top5[!duplicated(top5$gene),]

pdf("output/ExtFig2_Dotplot_markers.pdf", width = 9, height = 18)
DotPlot(object = MEPV_cells, 
        features = (top5$gene), 
        group.by = "Level2",
        assay = "SCT",
        scale = T,
        col.max = 2, 
        col.min = -2, cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size=9))+
  labs(title = "", x = "", y = "") +
  guides(colour = guide_colorbar(title = "Scaled average expression", 
                                 order = 1)) 
dev.off()

saveRDS(MEPV_cells, "MEPV_cells.rds")

######Extended Figure 6 ###########

pdf("output/ExtFig6_MEPV_Irf7.pdf", width = 9, height = 4.5)
FeaturePlot_scCustom(seurat_object = MEPV_cells, features = c("Irf7"),
                     num_columns = 4, reduction = "humap", alpha_exp = 0.75, na_color = "lightgray", na_cutoff = 0.1, split.by = "Sample.name_2") & NoAxes()
dev.off()

#####Extended Figure 7 #####

FeaturePlot_scCustom(seurat_object = MEPV_cells, features = c("Tlr2", "Tlr3", "Tlr4"),
                     num_columns = 1, reduction = "humap", alpha_exp = 0.75, na_color = "lightgray", na_cutoff = 0.5) 

###### Extended Figure 8 ######
pdf("output/Feature_MEPV_Cxcl10.pdf", width = 9, height = 4.5)
FeaturePlot_scCustom(seurat_object = MEPV_cells, features = c("Cxcl10"),
                     num_columns = 4, reduction = "humap", alpha_exp = 0.75, na_color = "lightgray", na_cutoff = 0.1, split.by = "Sample.name_2") & NoAxes()
dev.off()

#####Extended Figure 10######

pdf("output/Ifng_tbet_MEPV.pdf", width = 7.23, height = 5.85)
FeaturePlot_scCustom(seurat_object = MEPV_cells, features = c("Ifng", "Tbx21", "Ifngr1", "Ifngr2"),
                     num_columns = 2, reduction = "humap", alpha_exp = 0.75, na_color = "lightgray", na_cutoff = 1)
dev.off()

#####Extended Figure 13b######
pdf("output/Feature_MEPV_Hcar2.pdf", width = 5.1, height = 4.1)
FeaturePlot_scCustom(seurat_object = MEPV_cells, features = c("Hcar2"),
                     num_columns = 1, reduction = "humap", alpha_exp = 0.75, na_color = "lightgray", na_cutoff = 0.1) 

dev.off()

