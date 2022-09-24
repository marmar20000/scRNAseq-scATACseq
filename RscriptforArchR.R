library(ArchR)
library(parallel)
library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
#ArchR::installExtraPackages()
########################################################################################################################################################################################################################
#setwd("~/Documents/Marianna")
#InputFiles <- getInputFiles("toservernew")
setwd(("/c/Users/maria/Documents/Πτυχιακή/scRNA+scATAC"))
#setwd("/data/SeqFac/Lavigne/NGS/DATA_IT_sc/zip_files/unzipedATAC/ATAC/Cellranger/cr_count/MUC23141/outs")
#InputFiles <- getInputFiles("fragSamples2")
#InputFiles
#length(InputFiles)
addArchRGenome("mm10")
#addArchRThreads(threads = 16)
#ArrowFiles <- createArrowFiles(
#  inputFiles = (InputFiles),
#  sampleNames = names(InputFiles),
#  minTSS = 2, 
#  minFrags  = 100, 
#  addTileMat = TRUE,
#  addGeneScoreMat = TRUE)


ArrowFiles <- c("MUC23141NEW.arrow", "MUC23142NEW.arrow")

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
########################Creating ArchR project to work on###############################
projColonCancerMerged<- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "colonCancerWholeProjectMergedNew1",
  copyArrows = TRUE) #This is recommened so that you maintain an unaltered copy for later usage.



projColonCancerMergedNew
##Output
###class: ArchRProject 
#outputDirectory: /home/user/Documents/Marianna/colonCancerATAC 
#samples(2): cc_atac_MUC23141 cc_atac_MUC23142
#sampleColData names(1): ArrowFiles
#cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment BlacklistRatio
#numberOfCells(1): 2573
#medianTSS(1): 8.119
#medianFrags(1): 3667


getAvailableMatrices(projColonCancerMerged)
head(projColonCancerMerged$cellNames)
head(projColonCancerMerged$Sample)
quantile(projColonCancerMerged$TSSEnrichment)




WT <- BiocGenerics::which(projColonCancerMerged$Sample %in% "MUC23141NEW")
cellsSampleWT <- projColonCancerMerged$cellNames[WT]
KO <- BiocGenerics::which(projColonCancerMerged$Sample %in% "cc_atac_MUC23142")
cellsSampleKO <- projColonCancerMerged$cellNames[KO]
plotEmbedding(projColonCancerMerged[cellsSampleWT, ], ....)



########### After we got the information form computing doublets scores, we can remove them######
projColonCancerMerged <- filterDoublets(projColonCancerMerged, filterRatio = 1.5)
saveArchRProject(ArchRproj = projColonCancerMerged)
###############################Computing LSI####################################


projColonCancerMerged<- addIterativeLSI(
  ArchRProj = projColonCancerMerged,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI_tile", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
    ), 
  varFeatures = 25000, 
  dimsToUse = 1:30)

###########################CLUSTERING SEURAT#######################################

projColonCancerMerged<- addClusters(
  input = projColonCancerMerged,
  reducedDims = "IterativeLSI_tile",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.25, force = TRUE)
head(projColonCancerMerged$Clusters)
table(projColonCancerMerged$Clusters)
cM <- confusionMatrix(paste0(projColonCancerMerged$Clusters), paste0(projColonCancerMerged$Sample))
cM
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
pdf(file= "clustersPerSample.pdf", width = 4, height = 4)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p
dev.off()


##################################UMAP################################################

projColonCancerMerged<- addUMAP(
  ArchRProj = projColonCancerMerged, 
  reducedDims = "IterativeLSI_tile", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine", force = TRUE)

sum(table(projColonCancerMerged$CellTypeByCondition))

pdf(file =  "Plot-UMAP-Sample-Clusters.pdf", width = 5, height = 5)
p1 <- plotEmbedding(projColonCancerMerged, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p1
p2 <- plotEmbedding(ArchRProj = projColonCancerMerged, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2
p3 <- plotEmbedding(projColonCancerMerged[cellsSampleWT, ], colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p4 <- plotEmbedding(projColonCancerMerged[cellsSampleKO, ], colorBy = "cellColData", name = "Clusters", embedding = "UMAP")



ggAlignPlots(p3,p4 type = "h")
dev.off()



############################## Integrating With RNA seq #########################################################
seRNA <- readRDS("/data/SeqFac/Lavigne/NGS/DATA_IT_sc/zip_files/unzipedATAC/ATAC/Cellranger/cr_count/MUC23141/outs/BecauseStupid.rds")
seRNA
meta <- readRDS("/data/SeqFac/Lavigne/marianna/SETDBexpRNArenamed015Meta.rds")

rse <- SummarizedExperiment(assays = SimpleList(counts = seRNA), colData = seRNA@meta.data)
table(colData(rse)$new.ident)
seRNA1 =  CreateSeuratObject(
  counts=seRNA@assays$RNA@counts,
  meta.data=seRNA@meta.data)

projColonCancerMerged <- addGeneIntegrationMatrix(
  ArchRProj = projColonCancerMerged, 
 useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI_tile",
  seRNA = seRNA1,
  addToArrow = TRUE,
  groupRNA = "new.ident",
 nameCell = "predictedCell",
 nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force=TRUE)


################Gene Scores and Marker Genes with ArchR############################################


projColonCancerMerged <- addCellColData(ArchRProj = projColonCancerMerged, data = paste0(projColonCancerMerged$predictedGroup, "_", projColonCancerMerged$Sample),
    cells = projColonCancerMerged$cellNames, name = "CellTypeByCondition", force = TRUE)


markersGS <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "GeneIntegrationMatrix", 
  groupBy = "Clusters",  
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerList <- getMarkers(markersGS, cutOff = "Pval < 0.05 & Log2FC >= 0.58")
#write.table(markerList, file= "markers_Genes_CONDITION.tsv", sep = "\t")

pdf("heatmapGeneMarkers.pdf", width = 8, height = 6)
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "Pval < 0.05 & Log2FC >= 0.58", 
  transpose = TRUE, plotLog2FC = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


############################## PEAKS CALLING ############################################################

##############################################################################################################################
############################# MARKER PEAKS #################################################################################
##########################     PER CLUSTER #############################################################################
##############################################################################################################################



idxSample <- BiocGenerics::which(projColonCancerMerged$predictedGroup %in% "stem")
idxSample2 <- BiocGenerics::which(projColonCancerMerged$predictedGroup %in% "progenitor")
idxSample3 <- BiocGenerics::which(projColonCancerMerged$predictedGroup %in% "enterocyte")
idxSample4 <- c(idxSample, idxSample2, idxSample3)
cellsSample <- projColonCancerMerged$cellNames[idxSample4]
projSub <- subsetCells(ArchRProj = projColonCancerMerged, cellNames = projColonCancerMerged$cellNames[idxSample4]
)


## create pseudo-bulk replicates
projSub <- addGroupCoverages(
  ArchRProj = projSub,
  groupBy = "predictedGroup",
  force = TRUE
)

## call peaks
pathToMacs2 <- findMacs2()

projSub <- addReproduciblePeakSet(
  ArchRProj = projSub,
  groupBy = "predictedGroup",
  pathToMacs2 =pathToMacs2,
  force = TRUE
)

## calculate peak matrix
projSub <- addPeakMatrix(projSub,
binarize = TRUE, force = TRUE)
getAvailableMatrices(projSub)



markersPeaks_Trajectory <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerpeakList <- getMarkers(markersPeaks_Trajectory, cutOff = "Pval <= 0.05 & Log2FC >= 0.25")
markerpeakList
#write.table(markerpeakList, file = "markerPeakPval005.tsv", sep = "\t")
#write.table(markerpeakList$progenitor, file = "markerPeak_progenitor_Pval005.tsv", sep = "\t")
#write.table(markerpeakList$stem, file = "markerPeak_stem_Pval005.tsv", sep = "\t")
#rite.table(markerpeakList$enterocyte, file = "markerPeak_enterocyte_Pval005.tsv", sep = "\t")
#write.table(markerpeakList$stem_cKO, file = "markerPeaklistMergedSTEM_cKO.tsv", sep = "\t")
#write.table(markerpeakList$progenitor_cKO, file = "markerPeaklistMergedPROGENITOR_cKO.tsv", sep = "\t")
#write.table(markerpeakList$goblet_cKO, file = "markerPeaklistMergedGOBLET_cKO.tsv", sep = "\t")
#write.table(markerpeakList$Paneth_cKO, file = "markerPeaklistMergedPANETH_cKO.tsv", sep = "\t")
#write.table(markerpeakList$enterocyte_cKO, file = "markerPeaklistMergedENTEROCYTE_cKO.tsv", sep = "\t")
#write.table(markerpeakList$enteroendocrine_cKO, file = "markerPeaklistMergedENTERONDOCRINE_cKO.tsv", sep = "\t")
#write.table(markerpeakList$Tuft_cKO, file = "markerPeaklistMergedTUFT_cKO.tsv", sep = "\t")
#write.table(markerpeakList$stem_WT, file = "markerPeaklistMergedSTEM_WT.tsv", sep = "\t")
#write.table(markerpeakList$progenitor_WT, file = "markerPeaklistMergedPROGENITOR_WT.tsv", sep = "\t")
#write.table(markerpeakList$goblet_WT, file = "markerPeaklistMergedGOBLET_WT.tsv", sep = "\t")
#write.table(markerpeakList$Paneth_WT, file = "markerPeaklistMergedPANETH_WT.tsv", sep = "\t")
#write.table(markerpeakList$enterocyte_WT, file = "markerPeaklistMergedENTEROCYTE_WT.tsv", sep = "\t")
#write.table(markerpeakList$enteroendocrine_WT, file = "markerPeaklistMergedENTERONDOCRINE_WT.tsv", sep = "\t")
#write.table(markerpeakList$Tuft_WT, file = "markerPeaklistMergedTUFT_WT.tsv", sep = "\t")
############### Motif Enrichment in Marker Peaks ####################################################################

projSub <- addMotifAnnotations(ArchRProj = projSub,
motifSet = "cisbp", name = "Motif", species = "Mus musculus", force = TRUE)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks_Trajectory,
  ArchRProj = projSub, 
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
)
enrichMotifs


#pdf(file = "Motifs-Enriched-MarkerFDR03.pdf", width = 8, height = 6)
#heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
#ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", 
#annotation_legend_side = "bot")
#dev.off()


projSub <- addMotifAnnotations(ArchRProj = projSub,
motifSet = "cisbp", name = "Motif", species = "Mus musculus", force = TRUE)
 

markersPeaks_Trajectory_perCondition <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellTypeByCondition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


enrichMotifsTrajectory <- peakAnnoEnrichment(
  seMarker = markersPeaks_Trajectory_perCondition,
  ArchRProj = projSub, 
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC >= 0.60"
)
enrichMotifsTrajectory


#############################################################################
#######
##############################################################################################################
##################################################################################
##############            PER CLUSTER PER CONDITION ###########################################################
######################       AGAIN PEAK CALLING   ########################################################
###############################################################################
#############################################################################


## create pseudo-bulk replicates
projSub <- addGroupCoverages(
  ArchRProj = projSub,
  groupBy = "CellTypeByCondition",
  force = TRUE
)

## call peaks
pathToMacs2 <- findMacs2()

projSub <- addReproduciblePeakSet(
  ArchRProj = projSub,
  groupBy = "CellTypeByCondition",
  pathToMacs2 =pathToMacs2,
  force = TRUE
)

## calculate peak matrix
projSub <- addPeakMatrix(projSub,
binarize = TRUE, force = TRUE)
getAvailableMatrices(projSub)



markersPeaks_PERCLUSTER_PERCONDITION <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellTypeByCondition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerListPerClusterPerCondition <- getMarkers(markersPeaks_PERCLUSTER_PERCONDITION, cutOff = "Pval <= 0.05 & Log2FC >= 0.58", returnGR = TRUE)
markerListPerClusterPerCondition



############################################################################################################################
#####################  Pairwise Testing Between Groups####################################################################
markerTeststemWT <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByCondition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "stem_MUC23141NEW",
  bgdGroups = "stem_MUC23142NEW")


markerTestProgenitorWT <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByCondition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "progenitor_MUC23141NEW",
  bgdGroups = "progenitor_MUC23142NEW")

markerTestEnterocyteWT <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByCondition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "enterocyte_MUC23141NEW",
  bgdGroups = "enterocyte_MUC23142NEW")

markerTeststemKO <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByCondition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "stem_MUC23142NEW",
  bgdGroups = "stem_MUC23141NEW")


markerTestProgenitorKO <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByCondition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "progenitor_MUC23142NEW",
  bgdGroups = "progenitor_MUC23141NEW")

markerTestEnterocyteKO <- getMarkerFeatures(
  ArchRProj = projSub, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByCondition",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "enterocyte_MUC23142NEW",
  bgdGroups = "enterocyte_MUC23141NEW")


projSub <- addMotifAnnotations(ArchRProj = projSub, motifSet = "cisbp", name = "Motif", force = TRUE)


motifsSTEM_WT <- peakAnnoEnrichment(
    seMarker = markerTeststemWT,
    ArchRProj = projSub,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
  )


motifsSTEM_KO <- peakAnnoEnrichment(
    seMarker = markerTeststemKO,
    ArchRProj = projSub,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
  )


motifsPROGENITOR_WT <- peakAnnoEnrichment(
    seMarker = markerTestProgenitorWT,
    ArchRProj = projSub,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
  )


motifsPROGENITOR_KO <- peakAnnoEnrichment(
    seMarker = markerTestProgenitorKO,
    ArchRProj = projSub,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
  )


motifsENTEROCYTE_WT <- peakAnnoEnrichment(
    seMarker = markerTestEnterocyteWT,
    ArchRProj = projSub,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
  )


motifsENTEROCYTE_KO <- peakAnnoEnrichment(
    seMarker = markerTestEnterocyteKO,
    ArchRProj = projSub,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
  )



summtoplot = cbind(enrichMotifs, motifsSTEM_WT, motifsSTEM_KO, motifsPROGENITOR_WT, motifsPROGENITOR_KO,
motifsENTEROCYTE_KO, motifsENTEROCYTE_WT)


#summtoplot2 = cbind(enrichMotifsTrajectory,trajectory_perCondition)#, motifsSTEM_WT, motifsSTEM_KO, motifsPROGENITOR_WT, motifsPROGENITOR_KO )



pdf(file = "allIwannasee100.pdf", width = 20, height = 15)
heatmapEM <- plotEnrichHeatmap(summtoplot,clusterCols = FALSE, n = 50, returnMatrix = FALSE,  transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()



motifPositions <- getPositions(projSub)
projSub <- addGroupCoverages(ArchRProj = projSub, groupBy = "CellTypeByCondition")
seFoot <- getFootprints(
  ArchRProj = projSub, 
  positions = motifPositions[markerMotifs], 
  groupBy = "CellTypeByCondition"
)

footprints = plotFootprints(
  seFoot = seFoot,
  ArchRProj = projSub, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5, 
  height = 10, width = 6
)



if("Motif" %ni% names(projSub@peakAnnotation)){
    projSub <- addMotifAnnotations(ArchRProj = projSub, motifSet = "cisbp", name = "Motif")
}

projSub <- addMotifAnnotations(ArchRProj = projSub,
motifSet = "cisbp", name = "Motif", species = "Mus musculus", force = TRUE)


projSub <- addBgdPeaks(projSub)


projSub <- addDeviationsMatrix(
  ArchRProj = projSub, 
  peakAnnotation = "Motif",
  force = TRUE, matrixName = "Motives"
)

plotVarDev <- getVarDeviations(projSub, name = "MotifMatrix", plot = TRUE)

plotVarDev

#plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projSub, addDOC = FALSE)

motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- getFeatures(projSub, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

p <- plotGroups(ArchRProj = projSub, 
  groupBy = "Clusters2", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projSub)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

p <- plotEmbedding(
    ArchRProj = projSub, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projSub)
)

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))



markerRNA <- getFeatures(projSub, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA

p <- plotEmbedding(
    ArchRProj = projSub, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projSub)
)


p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

markerRNA <- getFeatures(projSub, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
markerRNA

p <- plotEmbedding(
    ArchRProj = projHeme5, 
    colorBy = "GeneIntegrationMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    continuousSet = "blueYellow",
    imputeWeights = getImputeWeights(projHeme5)
)

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))


