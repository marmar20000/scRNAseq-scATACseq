
inputFiles <- c("MUC23141/MUC23141.fragments.tsv.gz", "MUC23142/MUC23142.fragments.tsv.gz")
names(inputFiles) <- c("WT","KO")

inputFiles
length(inputFiles)
addArchRGenome("mm10")
addArchRThreads(threads = 8)
ArrowFiles <- createArrowFiles(
  inputFiles = (inputFiles),
  sampleNames = names(inputFiles),
  minTSS = 4, 
  minFrags  = 2000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE)


ArrowFiles

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
########################Creating ArchR project to work on###############################

proj_subset<- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "sc_atac_filter_in_2000_for_peak_calling",
  copyArrows = TRUE) 

proj_subset

######################## Plotting Sample Statistics from an ArchRProject ###############################
p1 <- plotGroups(
  ArchRProj = proj_subset, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges",
  baseSize = 10
)


p2 <- plotGroups(
  ArchRProj = proj_subset, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE,
)

p3 <- plotGroups(
  ArchRProj = proj_subset, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
  baseSize = 10
)

p4 <- plotGroups(
  ArchRProj = proj_subset, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 10,
  addBoxPlot = TRUE
)

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_subset, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj_subset)
p2 <- plotTSSEnrichment(ArchRProj = proj_subset)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj_subset, addDOC = FALSE, width = 5, height = 5)
proj_subset <- saveArchRProject(ArchRProj = proj_subset, outputDirectory = "Save-LAST", load = TRUE)

########### After we got the information form computing doublets scores, we can remove them######
proj_subset <- filterDoublets(proj_subset, filterRatio = 1.5)


########### ########### ########### ########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 
proj_subset <- addIterativeLSI(
  ArchRProj = proj_subset,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)
proj_subset <- addHarmony(
  ArchRProj = proj_subset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)
proj_subset <- addClusters(
  input = proj_subset,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
cM <- confusionMatrix(paste0(proj_subset$Clusters), paste0(proj_subset$Sample))

proj_subset <- addUMAP(
  ArchRProj = proj_subset, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
p1 <- plotEmbedding(ArchRProj = proj_subset, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj_subset, addDOC = FALSE, width = 5, height = 5)


plot_UMAP = plotEmbedding(
  ArchRProj = proj_subset,
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "CellTypeByCondition",
  size = 1,
  sampleCells = NULL,
  highlightCells = getCellNames(ArchRProj = proj_subset[which(proj_subset@cellColData$Sample == "WT")]),
  baseSize = 10,
  plotAs = "points")

plotPDF(plot_UMAP, name = "Plot-UMAP-WT.pdf", ArchRProj = proj_subset, addDOC = FALSE, width = 5, height = 5)


plot_UMAP = plotEmbedding(
  ArchRProj = proj_subset,
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "CellTypeByCondition",
  size = 1,
  sampleCells = NULL,
  highlightCells = getCellNames(ArchRProj = proj_subset[which(proj_subset@cellColData$Sample == "KO")]),
  baseSize = 10,
  plotAs = "points")

plotPDF(plot_UMAP, name = "Plot-UMAP-KO.pdf", ArchRProj = proj_subset, addDOC = FALSE, width = 5, height = 5)



########### ########### ########### ########### ########### ########### ########### ########### 
########### ###########   Sc RNA Integration    ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 


seRNA <- readRDS("Last_Project/Lgr5Cre_MERGED_Renamed.rds")
seRNA



seRNA1 =  CreateSeuratObject(
  counts=seRNA@assays$RNA@counts,
  meta.data=seRNA@meta.data)



proj_subset <- addGeneIntegrationMatrix(
  ArchRProj = proj_subset, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA1,
  addToArrow = TRUE,
  groupRNA = "CellType",
  nameCell = "predictedCell",
  nameGroup = "Clusters_of_RNA",
  nameScore = "predictedScore",
  force=TRUE)


proj_subset <- addCellColData(ArchRProj = proj_subset, data = paste0(proj_subset$Clusters_of_RNA, "_", proj_subset$Sample),
                                        cells = proj_subset$cellNames, name = "CellTypeByCondition", force = TRUE)

(table(proj_subset$CellTypeByCondition))
saveArchRProject(ArchRProj = proj_subset )



########### ########### ########### ########### ########### ########### ########### ########### 
########### ###########  PEAK CALLING  ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 



high_frag_cells <- rownames(proj_subset@cellColData[proj_subset$nFrags >= 2000, ])


selected_cells <- rownames(proj_subset1@cellColData[
  proj_subset$CellTypeByCondition %in% c("Stem I_KO", "Stem I_WT"), ])


proj_subset <- addGroupCoverages(
  ArchRProj = proj_subset,
  groupBy = "CellTypeByCondition",
  force = TRUE)

pathToMacs2 <- findMacs2()
# 
proj_subset <- addReproduciblePeakSet(
  ArchRProj = proj_subset,
  groupBy = "CellTypeByCondition",
  peakMethod = "Macs2", pathToMacs2 = pathToMacs2
)


# Add peak matrix
proj_subset <- addPeakMatrix(proj_subset)

# Save updated project
saveArchRProject(proj_subset)
getAvailableMatrices(proj_subset)









getGroupBW(
  ArchRProj = proj_subset,
  groupBy = "CellTypeByCondition",
  normMethod = "ReadsInTSS",
  tileSize = 10,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)






########### ########### ########### ########### ########### ########### ########### ########### 
########### ###########  MOTIVES    ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 



proj_subset2 <- addMotifAnnotations(ArchRProj = proj_subset2, motifSet = "cisbp", name = "Motif", force = TRUE)
proj_subset2 <- addBgdPeaks(proj_subset2)

proj_subset2 <- addDeviationsMatrix(
  ArchRProj = proj_subset2, 
  peakAnnotation = "Motif",
  force = TRUE)




########### ########### ########### ########### ########### ########### ########### ########### 
########### ###########  POSITIVE REGULTORS    ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 

seGroupMotif <- getGroupSE(ArchRProj = proj_subset1, useMatrix = "MotifMatrix", groupBy = "CellTypeByCondition")

corGIM_MM <- correlateMatrices(
  ArchRProj = proj_subset1,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI")


rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.05 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.1))] <- "YES"


positive_regulators = (corGIM_MM[corGIM_MM$TFRegulator=="YES",2])
positive_regulators

write.csv(corGIM_MM, 'motif_gene_expr_corellation_matrix_filtered_over_3000_please_stem_only.csv' )





cellMeta <- getCellColData(proj_subset2)
unique(cellMeta$CellTypeByCondition)
table(proj_subset2$CellTypeByCondition)



seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
seZ

# Extract the assay matrix
z_scores <- assay(seZ, "MotifMatrix")
# Ensure rownames and colnames are properly set
rownames(z_scores) <- rowData(seZ)$name  
colnames(z_scores) <- colnames(seZ)   
write.csv(z_scores, "z_scores_of_deviant_TF_0.5_filtered_above_3000_stem_only_second.csv")
z_scores_subset <- z_scores[rownames(z_scores) %in% positive_regulators, ]
z_scores_subset
write.csv(z_scores_subset, "z_scores_of_deviant_TF_0.5_padj_0_05_stem_only_positive_regulators_above_3000_macs2.csv")



# Extract the assay matrix
deviations <- assay(se_deviations, "MotifMatrix")
# Ensure rownames and colnames are properly set
rownames(deviations) <- rowData(se_deviations)$name  
colnames(deviations) <- colnames(se_deviations)   
write.csv(deviations, "deviations_of_deviant_TF_0.5_filtered_above_3000.csv")
deviations_subset <- deviations[rownames(deviations) %in% positive_regulators, ]
deviations_subset
write.csv(deviations_subset, "deviations_of_deviant_TF_0.5_padj_0_05_stem_only_positive_regulators_above_3000_macs2.csv")


markerMotifs = c('Atoh1_94', 'Sox9_725', 'Onecut2_827', 'Mafk_102', 'Ppara_670',
                 'Nr3c1_673', 'Bhlha15_87', 'Maf_793', 'Foxa2_307', 'Foxa3_831', 'Mafb_136', 'Arid3a_7', 'Esrra_674')

markerMotifs <- getFeatures(proj_subset1, select = paste(markerMotifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
keepCells <- rownames(cellMeta)[
  (grepl("_WT$", cellMeta$CellTypeByCondition)  & 
     cellMeta$CellTypeByCondition != "Unclassified_WT" )  # Optional: exclude this too
]


selected_cells <- rownames(proj_subset1@cellColData[
  proj_subset1$CellTypeByCondition %in% c("Stem I_KO", "Stem I_WT"), ])

for (i in markerMotifs) {
  
  p <- plotGroups(ArchRProj = proj_subset1[selected_cells,], 
                  groupBy = "CellTypeByCondition", 
                  colorBy = "MotifMatrix", 
                  name = i, plotAs = "ridges",   imputeWeights = getImputeWeights(proj_subset1[selected_cells,])) 
  
  # Define the desired order of clusters
  desiredOrder <- c(
    "Stem I_WT",  "Stem I_KO", "Stem II_WT", "Stem III_WT", "Progenitor I_WT", "Progenitor II_WT",
    "Ent.Immature_WT", "Ent.Mature_WT", "Paneth_WT",  "Goblet III_WT",
    "Tuft_WT",   "Enteroendocrine I_WT" )
  
  # Filter rows to keep only desired clusters
  df <- p[[1]]
  df_filtered <- df[df$x %in% desiredOrder, ]
  
  
  p <- ggGroup(
    x = as.character(df_filtered$x),
    y = df_filtered$y,
    groupOrder = rev(desiredOrder),
    xlabel = "CellTypeByCondition",
    ylabel = "Z-scores ",
    baseSize = 10,
    addBoxPlot = TRUE,
    plotAs = "ridges",   ridgeScale = 3,   alpha = 2, 
    pal = paletteDiscrete(values = desiredOrder, set = "stallion"), title = i)
  pdf_file <- paste0("plot_macs2_wt", i)

  plotPDF(p, name = pdf_file, width = 5, height = 8, ArchRProj = proj_subset1, addDOC = TRUE)
}


########### ########### ########### ########### ########### ########### ########### ########### 
########### ########### PEAK TO GENES ANALYSIS  ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 
cellMeta = getCellColData(project)
keepCells <- rownames(cellMeta)[
  (cellMeta$CellTypeByCondition == "Stem I_KO") |
    (cellMeta$CellTypeByCondition == "Stem I_WT") 
]


project <- addPeak2GeneLinks(project, reducedDims = "IterativeLSI",   k = 40,  
                                   useMatrix = "GeneIntegrationMatrix", cellsToUse  = keepCells,   addPermutedPval = TRUE, 
                                   predictionCutoff = 0.3)

p2g <- getPeak2GeneLinks(
  ArchRProj = project,
  corCutOff = 0.3,
  resolution = 1000,
  returnLoops = TRUE)
#p2g = p2g$Peak2GeneLinks
p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
p2g

palette <- c('Ent.Immature_KO'     =  "#D51F26"    ,  'Ent.Immature_WT' =  "#5A2955"    ,       'Ent.Mature_KO'     =   "#245359"   ,   'Ent.Mature_WT' =  "#9B8EC4" ,   'Enteroendocrine I_KO'  =   "#8AC972" , 
             
             'Enteroendocrine I_WT' =    "#FDDD03" ,  'Enteroendocrine II_KO' =   "#2B7F4A" ,    'Goblet I_KO'   =    "#C0545B"  ,          'Goblet I_WT'  =   "#6264A0"   ,
             
             'Goblet III_KO'   =   "#F69421"   ,        'Goblet III_WT'      =   "#B4B883" ,         'Paneth_KO'   =     "#DCABCF"   , 'Paneth_WT'  =  "#C39C6F" ,          'Progenitor I_KO' =     "#8ED2D0"  ,
             
             'Progenitor I_WT'    =  "#BFCADF"     ,    'Progenitor II_KO'    =   "#DB7D8D" ,    'Progenitor II_WT'  =  "#9C82BA" ,      'Stem I_KO'     =   "#3EB3A7"  ,      'Stem I_WT'  =   "#C16FAC" , 
             
             'Stem II_KO'     =     "#BF5D58"   ,         'Stem II_WT'   =  "#3E5D8D"  ,          'Stem III_KO'    =  "#216069"   ,         'Stem III_WT'      =    "#711E21" ,           'Tuft_KO'  = "#9A7456" ,
             
             'Tuft_WT'    =              "#B36B45"                ,     'Unclassified_KO'   =  "#AA875A" ,     'Unclassified_WT'  =  "#3D3D3D")

PLT = plotPeak2GeneHeatmap(
  ArchRProj = project,
  corCutOff = 0.3,
  FDRCutOff =1e-04,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  limitsATAC = c(-2, 2),
  limitsRNA = c(-2, 2),
 groupBy = 'CellTypeByCondition',
  k = 6, returnMatrices=FALSE,  palGroup = palette)


plotPDF('peak2genes_stem1_macs2_above_2000_cor_03.pdf', PLT)


matrix = plotPeak2GeneHeatmap(
  ArchRProj = project,
  corCutOff = 0.3,
  FDRCutOff =1e-04,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25, groupBy = 'CellTypeByCondition',
  k = 6, returnMatrices=TRUE)



matrix
rna_clusters <- matrix$RNA$kmeansId  
atac_clusters <- matrix$ATAC$kmeansId 
length(atac_clusters)
length(rna_clusters)
# Extract Peak-to-Gene links
p2g_links <- matrix$Peak2GeneLinks
length(p2g_links)

# Ensure gene and peak names are accessible
p2g_links$geneCluster <- rna_clusters 
p2g_links$peakCluster <- atac_clusters 

p2g_links

write.csv(p2g_links, 'p2glinks_stem_i_macs2_above_2000_cor_0_35.csv')




########### ########### ########### ########### ########### ########### ########### ########### 
########### ###########     MARKER MOTIVES      ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 


MOTIF_MARKERS_STEM_I_KO_WT = getMarkerFeatures(
  ArchRProj = proj_subset1,
  groupBy = "CellTypeByCondition",
  useGroups = 'Stem I_KO',
  bgdGroups = 'Stem I_WT',
  useMatrix = "MotifMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  normBy = NULL,
  testMethod = "wilcoxon",
  maxCells = 500,
  scaleTo = 10^4,
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  closest = FALSE,
  verbose = TRUE,
  logFile = createLogFile("getMarkerFeatures")
)


#markerList <- getMarkers(MOTIF_MARKERS_STEM_I_KO_WT, cutOff = "MeanDiff > 0 & Pval < 0.05")
#markerList$"Stem I_KO"$name
#[1] "NP_604"      "Ovol1_166"   "Zbtb3_235"   "Gm5294_338"  "Foxl1_364"
#[6] "Sox13_749"   "Trp63_853"   "Neurod2_69"  "Olig2_72"    "Arid5b_8"
#[11] "Neurog2_49"  "Neurog3_75"  "Hoxb3_546"   "Hoxd3_583"   "Sox30_736"
#[16] "Olig1_81"    "Atoh1_94"    "Mafg_123"    "Elf1_285"    "Evx1_410"
#[21] "Pdx1_472"    "Gbx2_494"    "Sox9_725"    "T_771"       "Neurod1_784"
#[26] "Neurod6_786" "Neurod4_787"





########### ########### ########### ########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 
########### ########### ########### ########### ########### ########### ########### ########### 

  
  cell_types <- c("Stem I_KO", "Stem I_WT")
  
  for (ct in cell_types) {
  
  message("Processing: ", ct)
  
  # Filter cells for current group
  cellMeta <- getCellColData(project)
  keepCells <- rownames(cellMeta)[cellMeta$CellTypeByCondition == ct]
  
  # Add peak2gene links only for current group
  project <- addPeak2GeneLinks(
    ArchRProj = project,
    reducedDims = "IterativeLSI",
    k = 40,
    useMatrix = "GeneIntegrationMatrix",
    cellsToUse = keepCells,
    addPermutedPval = TRUE,
    predictionCutoff = 0.3
  )
  
  # Get peak2gene links
  p2g <- getPeak2GeneLinks(
    ArchRProj = project,
    corCutOff = 0.55,
    resolution = 1000,
    returnLoops = FALSE
  )
  
  # Add gene and peak names
  p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
  p2g$peakName <- (metadata(p2g)$peakSet %>%
                     {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
  
  # Plot and save heatmap
  PLT <- plotPeak2GeneHeatmap(
    ArchRProj = project,
    corCutOff = 0.45,
    FDRCutOff = 1e-04,
    varCutOffATAC = 0.25,
    varCutOffRNA = 0.25,
    limitsATAC = c(-2, 2),
    limitsRNA = c(-2, 2),
    groupBy = 'CellTypeByCondition',
    k = 6,
    returnMatrices = FALSE,
    palGroup = palette
  )
  
  plotPDF(paste0("peak2genes_", gsub(" ", "_", tolower(ct)), ".pdf"), PLT)
  
  # Get matrices and clusters
  matrix <- plotPeak2GeneHeatmap(
    ArchRProj = project,
    corCutOff = 0.45,
    FDRCutOff = 1e-04,
    varCutOffATAC = 0.25,
    varCutOffRNA = 0.25,
    groupBy = 'CellTypeByCondition',
    k = 6,
    returnMatrices = TRUE
  )
  
  rna_clusters <- matrix$RNA$kmeansId
  atac_clusters <- matrix$ATAC$kmeansId
  p2g_links <- matrix$Peak2GeneLinks
  
  # Add cluster labels to links
  p2g_links$geneCluster <- rna_clusters
  p2g_links$peakCluster <- atac_clusters
  #p2g_links$Cluster <- atac_clusters
  
  # Export to CSV
  write.csv(p2g_links, paste0("p2glinks_", gsub(" ", "_", tolower(ct)), "_macs2_above_2000_cor_0_45.txt"))
  
  message("Finished: ", ct, "\n")
  }
  
  
  

p <- plotBrowserTrack(
  ArchRProj = project, 
  groupBy = "CellTypeByCondition", 
  useGroups = c('Stem I_KO', 'Stem I_WT'),
  geneSymbol = 'Itln1', 
  upstream = 100000,
  downstream = 100000,
  loops = getPeak2GeneLinks(project)
)
  


