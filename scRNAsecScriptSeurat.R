library(Seurat)
library(tidyverse)
library(magrittr)
library(purrr)
#######   Setting up the files to create A Seurat object to work on ###############


# step1 list sample directories ----------------------------------------------
dir.ls <- list.dirs(path = '/data/SeqFac/Lavigne/marianna/fastqsTalian/results',
            full.names = T,
            recursive = F)
dir.ls %<>% map( ~ paste0(.x, "/outs/filtered_feature_bc_matrix"))
names(dir.ls) <- c('sample21L0058', 'sample21L0062', 'sample21L0066', 'sample21L0070')

# step2 check whether dir exist -------------------------------------------
dir.ls %>% map( ~ dir.exists(.x))

# step3 create seurat per samples -----------------------------------------
obj.ls <- dir.ls %>% map( ~ Read10X(.x)) %>% map( ~ CreateSeuratObject(.x, min.cells = 3))

############## Merge samples into WildType and KnockOut##############################

WildType <- merge( x = obj.ls[[1]],y = obj.ls[[2]],add.cell.ids = c("WT_A", "WT_B"))
KnockOut <- merge( x = obj.ls[[3]],y = obj.ls[[4]],add.cell.ids = c("KO_A", "KO_B"))

str(WildType)
str(KnockOut)

metaWT <- WildType@meta.data 
metaKO <- KnockOut@meta.data
dim(metaWT)
dim(metaKO)
head(metaWT)
head(metaWT)
KnockOut[["origin"]] <- "KnockOut"
WildType[["origin"]] <- "WildType"



############################ Quality Control #######################################
####################################################################################


pdf(file= "QC_plotsWT.pdf" , width =  12, height = 16)
VlnPlot(WildType, features = c("nFeature_RNA","nCount_RNA"),ncol = 2,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
FeatureScatter(WildType, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(WildType, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()


pdf(file= "QC_plotsKO.pdf" , width = 12 , height = 16)
VlnPlot(KnockOut, features = c("nFeature_RNA","nCount_RNA"), ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
FeatureScatter(KnockOut, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(KnockOut, feature1 = "percent.rb", feature2 = "percent.mt")
dev.off()
##############The plots above clearly show that high MT percentage strongly correlates
# with low UMI counts, and usually is interpreted as dead cells. 
#High ribosomal protein content, however, strongly anti-correlates with MT, 
#and seems to contain biological signal. There’s also a strong correlation between 
#the doublet score and number of expressed genes.
#Let’s set QC column in metadata and define it in an informative way.
 ############################################################################################
 #######################################################################
table(rownames(WildType) %in% rownames(KnockOut)) 



############# Standard Process ###########################################
######### Normalize , HVG , Scaling ##############################################
############ Doublet Finder , filtering Doublets ################################

mouse_sample_list <- list()
mouse_sample_list[["WildType"]] <- WildType
mouse_sample_list[["KnockOut"]] <- KnockOut
QC
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
nExp_poi <- round(0.075*nrow(WildType@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
WildType_filtered <- doubletFinder_v3(WildType, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi)
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
nExp_poiKO <- round(0.075*nrow(KnockOut@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
KnockOut_filtered <- doubletFinder_v3(KnockOut, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi)
DF.name = colnames(WildType_filtered@meta.data)[grepl("DF.classification", colnames(WildType_filtered@meta.data))]
DF.name = colnames(KnockOut_filtered@meta.data)[grepl("DF.classification", colnames(KnockOut_filtered@meta.data))]
WildType_filtered = WildType_filtered[, WildType_filtered@meta.data[, DF.name] == "Singlet"]
KnockOut_filtered = KnockOut_filtered[, KnockOut_filtered@meta.data[, DF.name] == "Singlet"]
WildType_filtered <- subset(WildType_filtered, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20 & nCount_RNA < 25000)
KnockOut_filtered <- subset(KnockOut_filtered, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20 & nCount_RNA < 25000)






######################## Continuing With The **MERGED** Data Sets ################

mouse_seurat1=merge(WildType_filtered, y=KnockOut_filtered, add.cell.ids=c("WT","KO"), project="mouse_seurat")
all.genes <- rownames(mouse_seurat1)
 
# Violin Plot
pdf("violin_plot_mouse_beforenormseurat1.pdf", width=18, height=10)
VlnPlot(mouse_seurat1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

# Violin Plot with merged replicates
pdf("violin_plot_groupByBeforeNormmouse_seurat1.pdf", width=18, height=8)
VlnPlot(mouse_seurat1, features = c("nFeature_RNA", "nCount_RNA"), group.by = "origin", ncol = 2)
plot2 <- FeatureScatter(mouse_seurat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "origin")
plot2
dev.off()

# Feature Scatter Plot
#plot1 <- FeatureScatter(mouse_seurat1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_seurat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("feauture_scatterbeforenormmouseseurat1.pdf", width=15, height=8)
plot2
dev.off()

# Feature Scatter Plot
plot2 <- FeatureScatter(mouse_seurat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "origin")
pdf("QC_DIAGR.pdf", width=15, height=8)
plot2
dev.off()

mouse_seurat1 <- subset(mouse_seurat1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20 & nCount_RNA < 25000)


mouse_seurat1[["percent.mt"]] <- PercentageFeatureSet(mouse_seurat1, pattern = "^mt-")
mouse_seurat1 <- subset(mouse_seurat1, subset = nFeature_RNA > 400 & nFeature_RNA < 3500 & percent.mt < 5 & nCount_RNA < 20000 & nCount_RNA > 2000 )

mouse_seurat1 <- NormalizeData(mouse_seurat1, verbose = F)
mouse_seurat1 <- FindVariableFeatures(mouse_seurat1, selection.method = "vst", nfeatures = 2000, verbose = F)
mouse_seurat1 <- ScaleData(mouse_seurat1, verbose = F, features= all.genes)
mouse_seurat1 <- RunPCA(mouse_seurat1, npcs = 30,  features = VariableFeatures(object = mouse_seurat1))
mouse_seurat1 <- JackStraw(mouse_seurat1, num.replicate = 100)
mouse_seurat1 <- ScoreJackStraw(mouse_seurat1, dims = 1:20)


# Violin Plot
pdf("violin_plot_mouse_after_normseurat1.pdf", width=18, height=10)
VlnPlot(mouse_seurat1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

# Violin Plot with merged replicates
pdf("violin_plot_groupBy_after_Normmouse_seurat1.pdf", width=18, height=10)
VlnPlot(mouse_seurat1, features = c("nFeature_RNA", "nCount_RNA"), group.by = "origin", ncol = 2)
dev.off()

# Feature Scatter Plot
#plot1 <- FeatureScatter(mouse_seurat1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse_seurat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("feauture_scatter_after_normmouseseurat1.pdf", width=15, height=8)
plot2
dev.off()

# Feature Scatter Plot
plot2 <- FeatureScatter(mouse_seurat1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "origin")
pdf("feauture_scatter_merged_after_normmouseseurat1.pdf", width=15, height=8)
plot2
dev.off()



pdf("JackstrawPlot_mouseseurat1.pdf", width=15, height=8)
JackStrawPlot(mouse_seurat1, dims = 1:15)
dev.off()

pdf("ElbowPlot_mouseseurat1.pdf", width=15, height=8)
ElbowPlot(mouse_seurat1)
dev.off()

######################Continuing with clustering #########################################


mouse_seurat1 <- FindNeighbors(mouse_seurat1, dims = 1:8, graph.name = "test", k.param = 10, verbose = F)
mouse_seurat1 <- FindClusters(mouse_seurat1, verbose = F, graph.name = "test", resolution=0.13)
count_table <- table(mouse_seurat1@meta.data$seurat_clusters, mouse_seurat1@meta.data$origin)
count_table

################### Creating the first & second pdf with plots ############################

pdf(file = "PCA_Plots.pdf" , width = 4, height = 4)
DimHeatmap(mouse_seurat1, dims = 1:8, nfeatures = 20, cells = 500, balanced = T, split_by = "origin" )
DimPlot(mouse_seurat1, reduction = "pca", split_by = "origin")
dev.off()

mouse_seurat1 <- RunUMAP(mouse_seurat1, reduction = "pca", dims = 1:8, verbose = F)

pdf(file = "UMAP_Plotsmouseseurat1.pdf" , width = 6, height = 8)
DimPlot(mouse_seurat1, reduction = "umap", split.by = "origin") + NoLegend()
dev.off()
############################### we do the marker genes ####################################

all.markers1 <- FindAllMarkers(mouse_seurat1, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
write.table(all.markers, file = "all_markersRes0.13newsubclusters.csv", sep = ",", row.names = TRUE)
all.genes1 <- rownames(mouse_seurat1)
write.table(all.genes, file = "all_genes.tsv", sep = "\t", row.names = TRUE)
write.table(all.genes, file = "all_genes.csv", sep = ",", row.names = TRUE)

###################### Renaming Identities ##############################
mouse_seurat1[["oldold.ident"]] <- Idents(object = mouse_seurat1)
mouse_seurat1 <- RenameIdents(object = mouse_seurat1, "0" = "stem", "1" = "progenitor", "2" = "goblet", "3" = "enterocyte","4" = "enteroendocrine","5" = "Tuft","6" = "Paneth")
mouse_seurat1[["new.ident"]] <- Idents(object = mouse_seurat1)
saveRDS(mouse_seurat1, "BecauseStupid.rds")


pdf(file = "UMAP_PlotsRenamed.pdf" , width = 8 , height = 12)
DimPlot(mouse_seurat1, reduction = "umap", split.by = "origin") + plot_annotation(title = "sample cells after integration (Seurat 3)")
dev.off()

all.markers1 <- FindAllMarkers(mouse_seurat1, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
write.table(all.markers, file = "all_markersRes0.15.csv", sep = ",", row.names = TRUE)
markersstemKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="stem")
write.table(markersstemKOvsWT, file = 'markersstemKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersprogenitorKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="progenitor")
write.table(markersprogenitorKOvsWT, file = 'markersprogenitorKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersgobletKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="goblet")
write.table(markersgobletKOvsWT, file = 'markersgobletKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersenterocyteKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="enterocyte")
write.table(markersenterocyteKOvsWT, file = 'markersenterocyteKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersenteroendocrineKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="enteroendocrine")
write.table(markersenteroendocrineKOvsWT, file = 'markersenteroendocrineKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersTuftKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="Tuft")
write.table(markersTuftKOvsWT, file = 'markersTuftKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersPanethKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="Paneth")
write.table(markersPanethKOvsWT, file = 'markersPanethKOvsWT.tsv', sep = "\t", row.names = TRUE)

c('Xcl1','Olfm4','Gzma','Mptx1','Spink4')
c('Lgals3','Kcnq1ot1','Rgs16','Fcgbp','S100a11')
c('Klrd1','Hmgb2','Fabp1','Gip','Nap1l1')
c('Dek','Psmb9','Zg16','Klk1','Rbp2')
c('Pycard','Lsm2','Glud1','Dock10','Chaf1b')
c('Msn','Chgb','Fabp2','Trp53inp1')

markersstem <- FindMarkers(mouse_seurat1, ident.1 ="stem")
markersprogenitor <- FindMarkers(mouse_seurat1,ident.1 = "progenitor")
markersenterocyte <- FindMarkers(mouse_seurat1,ident.1 = "enterocyte")
write.table(markersstem, file = 'markersstemRNA.tsv', sep = "\t", row.names = TRUE)
write.table(markersprogenitor, file = 'markersprogenitorRNA.tsv', sep = "\t", row.names = TRUE)
write.table(markersenterocyte, file = 'markersenterocyteRNA.tsv', sep = "\t", row.names = TRUE)

write.table(markersstem, file = 'markersstemKO_ONLY.tsv', sep = "\t", row.names = TRUE)
markersprogenitorKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut', group.by = 'origin', subset.ident ="stem")
write.table(markersprogenitorKOvsWT, file = 'markersprogenitorKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersgobletKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="goblet")
write.table(markersgobletKOvsWT, file = 'markersgobletKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersenterocyteKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="enterocyte")
write.table(markersenterocyteKOvsWT, file = 'markersenterocyteKOvsWT.tsv', sep = "\t", row.names = TRUE)
markersenteroendocrineKOvsWT<- FindMarkers(mouse_seurat1,ident.1 = 'KnockOut',ident.2 = 'WildType', group.by = 'origin', subset.ident ="enteroendocrine")


mouse_seurat1 <- AddMetaData(mouse_seurat1, metadata = paste0(mouse_seurat1$new.ident, "_", mouse_seurat1$origin), col.name = "celltypeBycondition")


markersstemWT <- FindMarkers(mouse_seurat1, ident.1 ="stem_WildType")
markersprogenitorWT <- FindMarkers(mouse_seurat1,ident.1 = "progenitor_WildType")
markersenterocyteWT <- FindMarkers(mouse_seurat1,ident.1 = "enterocyte_WildType")
write.table(markersstemWT, file = 'markersstemRNA_ONLYWT.tsv', sep = "\t", row.names = TRUE)
write.table(markersprogenitorWT, file = 'markersprogenitorRNA_ONLYWT.tsv', sep = "\t", row.names = TRUE)
write.table(markersenterocyteWT, file = 'markersenterocyteRNA_ONLYWT.tsv', sep = "\t", row.names = TRUE)

markersstemKO <- FindMarkers(mouse_seurat1, ident.1 ="stem_knockOut")
markersprogenitorKO <- FindMarkers(mouse_seurat1,ident.1 = "progenitor_knockOut")
markersenterocyteKO <- FindMarkers(mouse_seurat1,ident.1 = "enterocyte_knockOut")
write.table(markersstemKO, file = 'markersstemRNA_ONLYKO.tsv', sep = "\t", row.names = TRUE)
write.table(markersprogenitorKO, file = 'markersprogenitorRNA_ONLYKO.tsv', sep = "\t", row.names = TRUE)
write.table(markersenterocyteKO, file = 'markersenterocyteRNA_ONLYKO.tsv', sep = "\t", row.names = TRUE)


########### dotplots per condition per cluster for top marker genes ###############

all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

pdf(file = "dotplotTOP10splitbyConditionALLnew.pdf", width = 35 , height = 20) 
DotPlot(object = mouse_seurat1, features = top10$gene,  cols = c("orange", "red"), split.by = "origin") &coord_flip()
dev.off()

all.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5

pdf(file = "dotplotTOP5splitbyConditionALLnew.pdf", width = 35 , height = 20) 
DotPlot(object = mouse_seurat1, features = top5$gene,  cols = c("orange", "red") ,split.by = "origin") &coord_flip() 
dev.off()


############################### Cell cycle Genes ###############################

cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

mouse_seurat <- CellCycleScoring(mouse_seurat, s.features = s.genes, g2m.features = g2m.genes)
table(mouse_seurat[[]]$Phase)

########## Plotting Technical Differences about mitochondria , ribosomic proteins & cell phase score
pdf(file = "CellCycleGenes_Plots.pdf" , width = 4, height = 4)
VlnPlot(mouse_seurat,features = c("nCount_RNA","nFeature_RNA"), group.by = "origin") &  theme(plot.title = element_text(size=10))
FeaturePlot(mouse_seurat,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T, group.by = "origin") & theme(plot.title = element_text(size=10))
VlnPlot(mouse_seurat,features = c("S.Score","G2M.Score") , group.by = "origin") & theme(plot.title = element_text(size=10))
dev.off()




##############################Barplot per centage cells per cluster

x = c("stem", "progenitor", "goblet", "enterocyte","enteroendocrine",  "Paneth", "Tuft"))
knockOut= c(31.15153251217416, 23.80406760240619, 17.23002005156116, 14.172156975078774, 5.4712116871956455, 5.012890289315383, 3.158120882268691)
WildType= c(36.41828069840664, 28.27974429920809, 15.084438507775976, 11.411124892662913, 2.8527812231657284, 3.301211716439271, 2.6524186623413795)
matrices <- cbind(knockOut, WildType)
matrices


barplot(matrices, main="Numbers of Cell per Cluster per Sample", 
        xlab = "number of cells", ylab = "samples",
        col = c("pink", "red", "cyan", "darkblue", "yellow", "purple", "green", "orange", "black"), 
        legend.text = x, beside = FALSE)


matrices

######################Find subclusters under one cluster############

for (cluster_to_subcluster in c("stem", "progenitor", "enterocyte")) {
  mouse_seurat1 <- FindSubCluster(mouse_seurat1, cluster_to_subcluster, "test", subcluster.name = sprintf("%s_subclusters", cluster_to_subcluster), resolution = 0.2, algorithm = 1)
}


pdf(file = ("stem_subclusters.pdf"), height = 8, width = 9)
DimPlot(mouse_seurat1, reduction = "umap", group.by = "stem_subclusters" , label = TRUE, label.size = 4, split.by = "origin")
dev.off()



pdf(file = ("progenitor_subclusters.pdf"), height = 8, width = 9)
DimPlot(mouse_seurat1, reduction = "umap", group.by = "progenitor_subclusters" , label = TRUE, label.size = 4, split.by = "origin")
dev.off()


pdf(file = ("enterocyte_subclusters.pdf"), height = 8, width = 9)
DimPlot(mouse_seurat1, reduction = "umap", group.by = "enterocyte_subclusters" , label = TRUE, label.size = 4, split.by = "origin")
dev.off()



for (i in c(count_table, count_table2, count_table3)) {
  sum(i[, "KnockOut"])
  sum(i[, "WildType"])


}




################# Dotplots Markers for each Cluster ##############################
#################### To check our biasted genes and SCSA algorithm ###############



# DotPlots per Cluster

pdf("Enterocyte_Immature Markers.pdf", width=15, height=8)
p =  DotPlot(mouse_seurat1,features = c("Reg3g", "Gsdmc4", "Prss32", "Krt8"),cols = c("blue","orange"))
p + ggtitle("Enterocyte Immature Markers") + RotatedAxis()
dev.off()

pdf("Enterocyte_Mature Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Elf3","Sis","Fabp1","Hnf4aos","Hnf4a", "Hnf4g", "Tmigd1","Fabp6","Slc51b","Slc51a","Mep1a","Fam151a","Naaladl1","Slc34a2","Plb1","Nudt4","Dpep1","Pmp22","Xpnpep2","Muc3","Neu1","Clec2h","Phgr1","Prss30","Aldob","Alpi","Apoa1","Apoa4","Lct"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Mature Markers")
dev.off()

pdf("Enterocyte Progenitor Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Ccnb1","Cdc20","Cenpa","Cdkn3","Cdc25c","Ccnb2","Kif22","Ube2c","Sapcd2","Rbp7","Ccna2","Aurka","Cdkn2d","Kif23","Nek2","Birc5","Plk1","Tacc3","Melk","Cps1"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Enterocyte Progenitor Markers")
dev.off()

pdf("Goblet Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Muc2","Spdef","Foxa1","Agr2","Spink4","Fcgbp","Tff3","Zg16","Clca1","Ccl6","Klk1","Tpsg1","Ccl9","Txndc5","Tspan13","Atoh1","Lrrc26","Clca3a1","Klf4"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Goblet Markers")
dev.off()

pdf("Paneth Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Lyz1","Mmp7","Dll4","Sox9","Gfi1","Gm14851","Defa21","Defa22","Defa17","Defa24","Defa3","Mptx2","Ang4"),cols = c("blue","orange")) + RotatedAxis()
# Gene Gm15284 & Defa-rs1 excluded as there are not in the dataset 
p + ggtitle("Paneth Markers")
dev.off()

pdf("Enteroendocrine Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Gfi1","Neurog3","Neurod1","Chga","Chgb","Isl1","Arx","Pax6","Foxa2","Sst","Gck","Gcg","Tph1","Pyy","Gfra3","Cpe","Tac1","Fam183b", "Hmgn3","Cck","Fev","Gch1","Pcsk1n", "Bex2","Vwa5b2","Nkx2-2","Marcksl1","Neurod2","Insm1"),cols = c("blue","orange")) + RotatedAxis()
# Gene Ngfrap1 excluded. It can not be found in the dataset
p + ggtitle("Enteroendocrine Markers")
dev.off()

pdf("Stem Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Lgr5","Ascl2","Olfm4","Prom1","Axin2","Fzd2","Fzd7","Lrp5","Lrp6","Notch1","Hes1","Smo","Yap1","Igfbp4","Bex1","Gkn3","Slc12a2"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Stem Markers")
dev.off()

pdf("Tuft Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Dclk1","Ptprc","Avil","Lrmp","Alox5ap","Rgs13","Sh2d6","Ltc4s","Hck","Cd24a","Trpm5","Kctd12","Aldh2","Il13ra1","Gng13","Tmem176a","Skap2","Ptpn6","Ly6g6f","Fyb","Adh1","Gfi1b","Il25"),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Tuft Markers")
dev.off()

pdf("Necroptosis Markers.pdf", width=15, height=8)
p = DotPlot(mouse_seurat1,features = c("Zbp1", "Ripk1", "Ripk3", "Cxcl1", "Ccl20", "Tnf", "Csf1" ),cols = c("blue","orange")) + RotatedAxis()
p + ggtitle("Necroptosis Markers")
dev.off()



pdf(file = "ViolinGene_Plots.pdf" , width = 4, height = 4)
VlnPlot(mouse_seurat1,features = c("Rgcc",	"Cd74",	"Akr1c13"	,"Cd74"	), group.by = "origin") &  theme(plot.title = element_text(size=10))
VlnPlot(mouse_seurat1,features = c("Mmp7"	,"Tstd1"	,"Adh6a",	"Rgs2"		), group.by = "origin") &  theme(plot.title = element_text(size=10))
VlnPlot(mouse_seurat1,features = c(	"Tff3",	"Ccl5",	"Crip1"	,"Mptx2"), group.by = "origin") &  theme(plot.title = element_text(size=10))
VlnPlot(mouse_seurat1,features = c(		"Zg16",	"Reg4"	,"Muc2"	), group.by = "origin") &  theme(plot.title = element_text(size=10))
dev.off()