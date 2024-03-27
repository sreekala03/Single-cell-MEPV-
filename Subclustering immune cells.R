######Subclustering immune cells#######


options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 3300 * 1024^2)
set.seed(1234)
library(dplyr)
library(Seurat)
library(ggplot2)
library(Nebulosa)
library(harmony)
library(readr)
library(RColorBrewer)
library(dittoSeq)
library(scCustomize)

##Load main object
MEPV_cells <- readRDS("MEPV_cells.rds")

###Subset Immune cells
Idents(object = MEPV_cells) <- "Level2"
Immune_cells <-subset(MEPV_cells, idents=c("Immune cells"))

Immune_cells <- FindNeighbors(object = Immune_cells, dims = 1:7, reduction = "harmony") 
Immune_cells <- FindClusters(object = Immune_cells, resolution = 0.3) 
Immune_cells <- RunUMAP(object = Immune_cells, dims = 1:7,reduction = "harmony", reduction.name = "humap",reduction.key = "hUMAP_", metric="euclidean")

Idents(Immune_cells) <- "seurat_clusters"

######Figure 4a#####
pdf("output/DimPLot_Immune.pdf", width = 4.2, height = 3.1)
DimPlot(Immune_cells,reduction = "humap", cols = c("0"="#FFC300", "1"="#1ABC9C", "2"="#FF5733"), label=F)
dev.off()

##Find Markers
Immune_cells <- PrepSCTFindMarkers(Immune_cells, assay = "SCT", verbose = TRUE)
Immune.Markers<-FindAllMarkers(Immune_cells, only.pos = TRUE)

Idents(object = Immune_cells) <- "seurat_clusters"


####Extended Figure 9 ######
pdf("output/Dotplot_Immune.pdf", width = 5.5, height = 6.3)
DotPlot(object = Immune_cells, 
        features = c("Ly6a", "Cd3g", "Ifi27l2a", "Txn1", "Ifitm3", "Gata3", "Selenof", "Ckb", "Txnip", "Aplp2", "Tln1", "Calm2", "Ptma", "Hspa1a", "H2-Q7",
                     "Cd3e", "Cd3d", "Cd8a", "Cd4", "Il7r", "Klre1", "Klrd1", "Klrk1", "Klrb1b", "Klrb1a", "Nkg7", "Fcer1g"),
        assay = "SCT",
        scale = T,
        cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size=9))+
  labs(title = "", x = "", y = "") +
  guides(colour = guide_colorbar(title = "Scaled average expression",
                                 order = 1)) + coord_flip()
dev.off()

####Figure 4b######
pdf("output/Immune_dot_estrous.pdf", width = 4.13, height =4.67)
DotPlot(object = Immune_cells_estr, 
        features = c(c("Nkg7", "Klrk1","Klrd1", "Klre1","Klrc2", "Klrb1a", "Klrb1b", "Cd8a", "Cd4", "Il7r","Gata3", "Cd3e")),
        assay = "SCT",
        scale = T,
        cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size=9))+
  labs(title = "", x = "", y = "") +
  guides(colour = guide_colorbar(title = "Scaled average expression",
                                 order = 1)) + coord_flip()
dev.off()