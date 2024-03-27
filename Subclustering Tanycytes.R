#########Tanycyte Subclustering#######

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
library(dittoSeq)
library(scCustomize)

##Load main object
MEPV_cells <- readRDS("MEPV_cells.rds")

Idents(object = Tanycytes) <- "Level2"
DimPlot(Tanycytes, reduction = "humap", label=F)

####Subset Tanycytes####
Idents(object = MEPV_cells) <- "seurat_clusters"
Tanycytes <-subset(MEPV_cells, idents=c(7, 2, 9, 13))

Tanycytes <- FindNeighbors(object = Tanycytes, dims = 1:28, reduction = "harmony") 
Tanycytes <- FindClusters(object = Tanycytes, resolution = 0.6) 
Tanycytes <- RunUMAP(object = Tanycytes, dims = 1:28,reduction = "harmony", reduction.name = "humap",reduction.key = "hUMAP_", metric="euclidean")

pdf(file = "output/DimPlot_Tany.pdf",   
    width = 6.2, 
    height = 5.4)
DimPlot(Tanycytes, reduction = "humap", label=TRUE)
dev.off()

##Find Markers
Tanycytes <- PrepSCTFindMarkers(Tanycytes, assay = "SCT", verbose = TRUE)
Tanycyte.Markers<-FindAllMarkers(Tanycytes, only.pos = TRUE)

#########Remove cluster 2 with Plp1 and cytoplasmic RNA ####
Idents(Tanycytes) <-  "seurat_clusters"
Tanycytes <-subset(Tanycytes, idents=c(2), invert = T)

Tanycytes <- FindNeighbors(object = Tanycytes, dims = 1:28, reduction = "harmony") 
Tanycytes <- FindClusters(object = Tanycytes, resolution = 0.6) 
Tanycytes <- RunUMAP(object = Tanycytes, dims = 1:28,reduction = "harmony", reduction.name = "humap",reduction.key = "hUMAP_", metric="euclidean")

DimPlot(Tanycytes, reduction = "humap", label=TRUE)

###Name tanycyte cell states##
Idents(Tanycytes) <-  "seurat_clusters"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(5)))              <- "Tany.1"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(3)))              <- "Tany.2"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(6)))              <- "Tany.3"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(1)))              <- "Tany.4"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(0)))              <- "Tany.5"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(4)))              <- "Tany.6"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(2)))              <- "Tany.7"
Idents(object = Tanycytes, cells = WhichCells(Tanycytes, ident = c(7)))              <- "Tany.8"

Tanycytes$Tany.cell.state <- Idents(object = Tanycytes)
plot1 = DimPlot(Tanycytes, reduction = "humap") & NoLegend() #& NoAxes() for plots without axes
LabelClusters(plot1, id = "ident", size = 4, repel = F) 

pdf("output/Fig1d_Dimplot_Tanycyte.pdf", width = 5.7, height = 4.9)
dittoDimPlot(Tanycytes, var = "Tany.cell.state",reduction.use = 'humap', opacity = 0.5)
dev.off()

####Dot plot of Tanycyte Markers#####

Idents(Tanycytes) <- "Tany.cell.state"
Tanycyte.sub.Markers<-FindAllMarkers(Tanycytes, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25)
write.csv(Tanycyte.sub.Markers, "output/TableS2_Tanycyte_markers.csv")

levels(Tanycytes) <- c('Tany.1', 'Tany.2', 'Tany.3', 'Tany.4','Tany.5','Tany.6','Tany.7', 'Tany.8')

pdf("output/ExtFig3_Dotplot_Tany_markers.pdf", width = 6.4, height = 10)
DotPlot(object = Tanycytes, 
        features = rev(c("Cd81", "Ntsr2", "Slc6a11", "Clu", "Gpr37l1", "Cspg5","Slc1a2", 
                         "Cpe", "Rspo3", "Tgfb2","Flt1", "Nr2f1",	"Timp4", "S100b","Gfap", "Cd24a",
                         "Rarres2", "Dbi", "Tmem212", "Gja1", "Fam183b", "Pltp", "Pcsk2", 
                         "Crym","Itm2b", "Id4", "Prdx4", "S100a6", "6330403K07Rik",
                         "Prdx6", "Thrsp","Ephb1","Vcan","Pdzph1", "Ptn", "Frzb", "Penk", "Mia", "Dlk1", 
                         "Col25a1", "Six6", "Scn7a", "Cldn10", "Mest", "Rgcc", "Fndc3c1", "Adm", "Sfrp2", "Oasl2",
                         "Ifitm3", "Ifit1", "Ifit3", "Igtp","Isg15", "Irf7", "Stat2", "Cxcl10")),
        assay = "SCT",
        scale = T,
        cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size=9))+
  labs(title = "", x = "", y = "") +
  guides(colour = guide_colorbar(title = "Scaled average expression", 
                                 order = 1)) +coord_flip()
dev.off()


saveRDS(Tanycytes, "Tanycytes.rds")

#####Density Plot of Markers ##########

pdf("output/ExtFig3_density_Cd81.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Cd81")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Gpr37l1.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Gpr37l1")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Cspg5.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Cspg5")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Tgfb2.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Tgfb2")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Cd24a.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Cd24a")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Pcsk2.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Pcsk2")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Ephb1.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Ephb1")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Vcan.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Vcan")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Dlk1.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Dlk1")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Sfrp2.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Sfrp2")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Thrsp.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Thrsp")) & NoAxes() & NoLegend()
dev.off()

pdf("output/ExtFig3_density_Stat2.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Stat2")) & NoAxes() & NoLegend()
dev.off()

###Density plot of inflammatory markers

pdf("MEPV/output/Fig1e_density_B2m.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("B2m")) & NoAxes() & NoLegend()
dev.off()

pdf("MEPV/output/Fig1e_density_Ncoa7.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Ncoa7")) & NoAxes() & NoLegend()
dev.off()

pdf("MEPV/output/Fig1e_density_Ly6e.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Ly6e")) & NoAxes() & NoLegend()
dev.off()

pdf("MEPV/output/Fig1a_density_Ifit3.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Ifit3")) & NoAxes() & NoLegend()
dev.off()

pdf("MEPV/output/Fig1a_density_Trim25.pdf", width = 4, height = 4)
plot_density(Tanycytes, c("Trim25")) & NoAxes() & NoLegend()
dev.off()


##########Functional Enrichment of Tany.8#######

library(fgsea)

# Load GSEA gene sets: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
#M2_curated_all <- fgsea::gmtPathways("MEPV/m2.all.v2023.2.Mm.symbols.gmt")
# Reactome <- fgsea::gmtPathways("MEPV/m2.cp.reactome.v2023.2.Mm.symbols.gmt")
# GO_BP <- fgsea::gmtPathways("MEPV/m5.go.bp.v2023.2.Mm.symbols.gmt")

Hallmark <- fgsea::gmtPathways("MEPV/mh.all.v2023.2.Mm.symbols.gmt")

Idents(Tanycytes) <- "Tany.cell.state"

Tany.sub.Marker <-FindAllMarkers(Tanycytes, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25)

Tany.8 = subset(Tany.sub.Marker, cluster == 'Tany.8')
Tany.8 <- Tany.8 %>% arrange(desc(avg_log2FC))
fold_changes_8 <- Tany.8$avg_log2FC
names(fold_changes_8) <- Tany.8$gene

gsea_Tany.8_H <- fgsea(pathways = Hallmark,
                       stats = fold_changes_8,
                       eps = 0.0,
                       minSize=10,
                       maxSize=500)

# head(gsea_Tany.8_R[order(pval), ])

df_gsea_Tany.8_H = gsea_Tany.8_H$leadingEdge

# make a table plot for a bunch of selected pathways:
topPathwaysUp8 <- subset(gsea_Tany.8_H, ES > 0 & pval < 0.05)

topPathwaysUp8_list = topPathwaysUp8$leadingEdge

write.csv(topPathwaysUp8_list[[3]], "output/genes_interferon_gamma_h.csv")
write.csv(topPathwaysUp8_list[[2]], "output/genes_interferon_alpha_h.csv")
write.csv(topPathwaysUp8_list[[1]], "output/genes_inflammatory response_h.csv")

### Figure 1F Bar plot###########
theme_set(theme_bw())
pdf("output/Gsea_Reactome_Tany8.pdf", width = 8.5, height = 2.7)
ggplot(topPathwaysUp8, aes (x = pathway, y = size, fill = padj)) + geom_bar (stat = "identity") + coord_flip()+ 
  scale_fill_gradient(low="red", high="pink", breaks = c(0.01,0.05,0.1))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 11),
        axis.text = element_text(size = 11, colour="black")) +labs(x = "", y = "Gene Count")
dev.off()


#####Figure 1G Dot plot of leading edge genes ########

pdf("output/Dotplot_genes_hallmark_interferon.pdf", width = 5, height = 7.5)
DotPlot(object = Tanycytes, 
        features = c("B2m", "Bst2", "Cmpk2", "Cxcl10", "Eif2ak2", "Gbp3", "Herc6", "Ifi27", "Ifi35", "Ifi44", "Ifih1",
                     "Ifit3", "Ifitm3", "Irf7", "Isg15", "Lgals3bp", "Ly6e", "Psmb8", "Psmb9", "Psme1", "Rsad2", 
                     "Rtp4", "Samd9l", "Tap1", "Usp18", "Stat1", "Xaf1", "Socs3", "Psmb10", "Samhd1", "Pnp", "Rnf213",
                     "Parp9", "Parp14", "Ube2l6", "Ddx60", "Irf9", "Adm", "Btg2", "Fzd5", "Klf6", "Tlr3", "Gpc3", "Irf1", "Csf1",
                     "Tapbp", "Nmi"),
        assay = "SCT",
        scale = T,
        cols = "RdYlBu") +
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.02) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5,
                                   size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size=9))+
  labs(title = "", x = "", y = "") +
  guides(colour = guide_colorbar(title = "Scaled average expression", 
                                 order = 1)) +coord_flip()
dev.off()


####Figure 1H Tanycyte cell state Proportions by conditions#####
#reorder sample conditions

Tanycytes$Sample.name_2 <- factor(Tanycytes$Sample.name_2,levels=c("Male.SD", "Male.HFDS", "Male.HFDR", "Fem.SD.Dioest","Fem.SD.Proest", "Fem.SD.Oest", "Fem.HFDS", "Fem.HFDR" ))

Idents(object = Tanycytes) <- "Tany.cell.state"
cells_by_cluster = as.matrix(table(Idents(Tanycytes), Tanycytes$Sample.name_2))

# 1. convert the data as a table
cells_by_cluster <- as.table(as.matrix(cells_by_cluster))
#Compute Chi-square residuals
chisq <- chisq.test(cells_by_cluster)
chisq

# Observed counts
chisq$observed

# Expected counts
round(chisq$expected,0)

round(chisq$residuals, 0)

library(corrplot)
pdf(file = "output/Fig1h_corrplot_Tany_cond.pdf",   # The directory you want to save the file in
    width = 6.5, # The width of the plot in inches
    height = 5.35) # The height of the plot in inches
corrplot(chisq$residuals, method = 'square', is.cor = FALSE, tl.col = 'black', cl.cex = 0.8, tl.cex = 0.8,col = rev(brewer.pal(n=8, name="RdGy")))
dev.off()

####Figure 2a Highlight Tany.8 #####

pdf("output/Fig2a_Highlight_Tany8.pdf", width = 6, height = 4.8)
DimPlot(Tanycytes, reduction = "humap", cells.highlight=WhichCells(Tanycytes, idents = c("Tany.8")), sizes.highlight = 1.5)
dev.off()

####Figure 2D####
Idents(MEPV_cells) <- "Sample.name_2"
pdf("Fig2d_Feature_Tany_Irf7.pdf", width = 6, height = 4.6)
FeaturePlot_scCustom(seurat_object = Tanycytes, features = c("Irf7"),
                     num_columns = 4, reduction = "humap", alpha_exp = 0.75, na_color = "lightgray", na_cutoff = 0.3, split.by = "Sample.name_2") & NoAxes()
dev.off()

######Extended Figure 4 ######

pdf("output/ExtFig4_Tany_Irfs.pdf", width = 6, height = 4.6)
FeaturePlot_scCustom(seurat_object = Tanycytes, features = c("Irf7", "Irf2", "Irf9"),
                     num_columns = 4, reduction = "humap", alpha_exp = 0.75, na_color = "lightgray", na_cutoff = 0.3, split.by = "Sample.name_2") 
dev.off()

