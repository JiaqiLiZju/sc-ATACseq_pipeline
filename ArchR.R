# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# # devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
# BiocManager::install("motifmatchr")
# BiocManager::install("chromVAR")
# devtools::install_local("~/Desktop/ArchR-master.zip")
# library(ArchR)
# ArchR::installExtraPackages()

setwd("/media/ggj/Files/scATAC/Sci1-ATAC/ArchR/")
library(ArchR)
set.seed(1)

addArchRThreads(threads = 16) 

addArchRGenome("mm9")

inputFiles <- "../mouse.sorted.bam"
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = "mouse.sorted.bam",
  filterTSS = 0, #Dont set this too high because you can always increase later
  filterFrags = 0,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)
proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)
p$CD14

#Rearrange for grid plotting
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

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$CD14)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

proj <- saveArchRProject(ArchRProj = proj)
Sys.Date()
sessionInfo()
