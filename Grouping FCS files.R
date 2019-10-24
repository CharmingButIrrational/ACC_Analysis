#Aggregate flow cytometry files
library(flowCore)
library(FlowSOM)
library(matrixStats)
library(ConsensusClusterPlus)

#Import data Group 1
group1 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 1.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Import data Group 2
group2 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 2.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Import data Group 3
group3 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 3.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Import data Group 4
group4 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 4.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Import data Group 5
group5 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 5 minus 20.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))

#Remove extra columns from data
marker_cols_G1 <- c(1,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,27,28,29,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69)
marker_cols_G2 <- c(1,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,27,28,29,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69)
marker_cols_G3 <- c(1,2,3,4,5,6,8,11,12,13,14,15,16,17,18,19,20,24,25,26,27,29,30,31,32,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,52,53,54,56,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
marker_cols_G4 <- c(1,2,3,4,5,6,8,11,12,13,14,15,16,17,18,19,20,24,25,26,27,29,30,31,32,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,52,53,54,56,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
marker_cols_G5 <- c(1,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,27,28,29,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69)


group1_sub <- group1[,c(marker_cols_G1)]
group2_sub <- group2[,c(marker_cols_G2)]
group3_sub <- group3[,c(marker_cols_G3)]
group4_sub <- group4[,c(marker_cols_G4)]
group5_sub <- group5[,c(marker_cols_G5)]

dim(group1)
dim(group1_sub)
dim(group2)
dim(group2_sub)
dim(group3)
dim(group3_sub)
dim(group4)
dim(group4_sub)
dim(group5)
dim(group5_sub)

#Transforming the data
asinh_scale <- 5
group1_sub[, 1:60] <- asinh(group1_sub[, 1:60] / asinh_scale)
group2_sub[, 1:60] <- asinh(group2_sub[, 1:60] / asinh_scale)
group3_sub[, 1:60] <- asinh(group3_sub[, 1:60] / asinh_scale)
group4_sub[, 1:60] <- asinh(group4_sub[, 1:60] / asinh_scale)
group5_sub[, 1:60] <- asinh(group5_sub[, 1:60] / asinh_scale)

#Scaling the data from 0-1
group1_sub <- apply(group1_sub, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
group2_sub <- apply(group2_sub, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
group3_sub <- apply(group3_sub, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
group4_sub <- apply(group4_sub, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
group5_sub <- apply(group5_sub, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

dim(group1_sub)
dim(group2_sub)
dim(group3_sub)
dim(group4_sub)
dim(group5_sub)

#Checking all the samples have the same columns
length(unique(lapply(list(group1_sub, group2_sub, group3_sub, group4_sub, group5_sub), dimnames))) == 1 
    #If FALSE then at least 1 set of dimnames is different 

#Use rbind to join the matrixes together
ComData <- rbind(group1_sub, group2_sub, group3_sub, group4_sub, group5_sub)

dim(ComData)

#Select the columns used for clustering
ClusterCols <- c(9,10,11,12,15,16,18,21,22,26,29,31,33,34,42,47,50,52,53,56,57,58,59,60)

#Build flowframe object for combind data
ComFlowSOM <- flowCore::flowFrame(ComData)

#Build flowframe objects for the individual groups
group1_FlowSOM <- flowCore::flowFrame(group1_sub)
group2_FlowSOM <- flowCore::flowFrame(group2_sub)
group3_FlowSOM <- flowCore::flowFrame(group3_sub)
group4_FlowSOM <- flowCore::flowFrame(group4_sub)
group5_FlowSOM <- flowCore::flowFrame(group5_sub)

#Run FlowSOM ComData
set.seed(234)
out <- FlowSOM::ReadInput(ComFlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = ClusterCols)
out <- FlowSOM::BuildMST(out)

#Visualisation of the combined groups 
PlotStars(out)

#Mapping the indivual groups onto the combine SOM structure 
Group1_remap <- NewData(out, group1_FlowSOM)
Group2_remap <- NewData(out, group2_FlowSOM)
Group3_remap <- NewData(out, group3_FlowSOM)
Group4_remap <- NewData(out, group4_FlowSOM)
Group5_remap <- NewData(out, group5_FlowSOM)

##############################################################
#################    Meta clustering    ######################
##############################################################
#Meta clustering
metaClusteringOut <- metaClustering_consensus(out$map$codes,k=10)
metaClusteringG1 <- metaClustering_consensus(Group1_remap$map$codes,k=10)
metaClusteringG2 <- metaClustering_consensus(Group2_remap$map$codes,k=10)
metaClusteringG3 <- metaClustering_consensus(Group3_remap$map$codes,k=10)
metaClusteringG4 <- metaClustering_consensus(Group4_remap$map$codes,k=10)
metaClusteringG5 <- metaClustering_consensus(Group5_remap$map$codes,k=10)

#Plot the metaclustering
PlotPies(out, cellTypes = out$map$mapping[,1], backgroundValues = as.factor(metaClusteringOut))
PlotPies(out, cellTypes = out$map$mapping[,1], backgroundValues = as.factor(metaClusteringOut), view = "grid")


##############################################################

#Plot group 1 markers
PlotStars(Group1_remap)

PlotMarker(Group1_remap,"Bi209Di")
PlotMarker(Group1_remap,"Dy161Di")
PlotMarker(Group1_remap,"Dy162Di")
PlotMarker(Group1_remap,"Dy163Di")
PlotMarker(Group1_remap,"Dy164Di")
PlotMarker(Group1_remap,"Er166Di")
PlotMarker(Group1_remap,"Er167Di")
PlotMarker(Group1_remap,"Er168Di")
PlotMarker(Group1_remap,"Er170Di")
PlotMarker(Group1_remap,"Eu151Di")
PlotMarker(Group1_remap,"Eu153Di")
PlotMarker(Group1_remap,"Gd155Di")
PlotMarker(Group1_remap,"Gd156Di")
PlotMarker(Group1_remap,"Gd160Di")
PlotMarker(Group1_remap,"Ho165Di")
PlotMarker(Group1_remap,"Lu175Di")
PlotMarker(Group1_remap,"Nd143Di")
PlotMarker(Group1_remap,"Nd145Di")
PlotMarker(Group1_remap,"Nd146Di")
PlotMarker(Group1_remap,"Nd148Di")
PlotMarker(Group1_remap,"Nd150Di")
PlotMarker(Group1_remap,"Pr141Di")
PlotMarker(Group1_remap,"Sm147Di")
PlotMarker(Group1_remap,"Sm149Di")
PlotMarker(Group1_remap,"Sm152Di")
PlotMarker(Group1_remap,"Sm154Di")
PlotMarker(Group1_remap,"Tb159Di")
PlotMarker(Group1_remap,"Tm169Di")
PlotMarker(Group1_remap,"Yb171Di")
PlotMarker(Group1_remap,"Yb172Di")
PlotMarker(Group1_remap,"Yb173Di")
PlotMarker(Group1_remap,"Yb174Di")
PlotMarker(Group1_remap,"Yb176Di")

#Plot group 2 markers
PlotStars(Group2_remap)

PlotMarker(Group2_remap,"Bi209Di")
PlotMarker(Group2_remap,"Dy161Di")
PlotMarker(Group2_remap,"Dy162Di")
PlotMarker(Group2_remap,"Dy163Di")
PlotMarker(Group2_remap,"Dy164Di")
PlotMarker(Group2_remap,"Er166Di")
PlotMarker(Group2_remap,"Er167Di")
PlotMarker(Group2_remap,"Er168Di")
PlotMarker(Group2_remap,"Er170Di")
PlotMarker(Group2_remap,"Eu151Di")
PlotMarker(Group2_remap,"Eu153Di")
PlotMarker(Group2_remap,"Gd155Di")
PlotMarker(Group2_remap,"Gd156Di")
PlotMarker(Group2_remap,"Gd160Di")
PlotMarker(Group2_remap,"Ho165Di")
PlotMarker(Group2_remap,"Lu175Di")
PlotMarker(Group2_remap,"Nd143Di")
PlotMarker(Group2_remap,"Nd145Di")
PlotMarker(Group2_remap,"Nd146Di")
PlotMarker(Group2_remap,"Nd148Di")
PlotMarker(Group2_remap,"Nd150Di")
PlotMarker(Group2_remap,"Pr141Di")
PlotMarker(Group2_remap,"Sm147Di")
PlotMarker(Group2_remap,"Sm149Di")
PlotMarker(Group2_remap,"Sm152Di")
PlotMarker(Group2_remap,"Sm154Di")
PlotMarker(Group2_remap,"Tb159Di")
PlotMarker(Group2_remap,"Tm169Di")
PlotMarker(Group2_remap,"Yb171Di")
PlotMarker(Group2_remap,"Yb172Di")
PlotMarker(Group2_remap,"Yb173Di")
PlotMarker(Group2_remap,"Yb174Di")
PlotMarker(Group2_remap,"Yb176Di")

#Plot group 3 markers
PlotStars(Group3_remap)

PlotMarker(Group3_remap,"Bi209Di")
PlotMarker(Group3_remap,"Dy161Di")
PlotMarker(Group3_remap,"Dy162Di")
PlotMarker(Group3_remap,"Dy163Di")
PlotMarker(Group3_remap,"Dy164Di")
PlotMarker(Group3_remap,"Er166Di")
PlotMarker(Group3_remap,"Er167Di")
PlotMarker(Group3_remap,"Er168Di")
PlotMarker(Group3_remap,"Er170Di")
PlotMarker(Group3_remap,"Eu151Di")
PlotMarker(Group3_remap,"Eu153Di")
PlotMarker(Group3_remap,"Gd155Di")
PlotMarker(Group3_remap,"Gd156Di")
PlotMarker(Group3_remap,"Gd160Di")
PlotMarker(Group3_remap,"Ho165Di")
PlotMarker(Group3_remap,"Lu175Di")
PlotMarker(Group3_remap,"Nd143Di")
PlotMarker(Group3_remap,"Nd145Di")
PlotMarker(Group3_remap,"Nd146Di")
PlotMarker(Group3_remap,"Nd148Di")
PlotMarker(Group3_remap,"Nd150Di")
PlotMarker(Group3_remap,"Pr141Di")
PlotMarker(Group3_remap,"Sm147Di")
PlotMarker(Group3_remap,"Sm149Di")
PlotMarker(Group3_remap,"Sm152Di")
PlotMarker(Group3_remap,"Sm154Di")
PlotMarker(Group3_remap,"Tb159Di")
PlotMarker(Group3_remap,"Tm169Di")
PlotMarker(Group3_remap,"Yb171Di")
PlotMarker(Group3_remap,"Yb172Di")
PlotMarker(Group3_remap,"Yb173Di")
PlotMarker(Group3_remap,"Yb174Di")
PlotMarker(Group3_remap,"Yb176Di")

#Plot group 4 markers
PlotStars(Group4_remap)

PlotMarker(Group4_remap,"Bi209Di")
PlotMarker(Group4_remap,"Dy161Di")
PlotMarker(Group4_remap,"Dy162Di")
PlotMarker(Group4_remap,"Dy163Di")
PlotMarker(Group4_remap,"Dy164Di")
PlotMarker(Group4_remap,"Er166Di")
PlotMarker(Group4_remap,"Er167Di")
PlotMarker(Group4_remap,"Er168Di")
PlotMarker(Group4_remap,"Er170Di")
PlotMarker(Group4_remap,"Eu151Di")
PlotMarker(Group4_remap,"Eu153Di")
PlotMarker(Group4_remap,"Gd155Di")
PlotMarker(Group4_remap,"Gd156Di")
PlotMarker(Group4_remap,"Gd160Di")
PlotMarker(Group4_remap,"Ho165Di")
PlotMarker(Group4_remap,"Lu175Di")
PlotMarker(Group4_remap,"Nd143Di")
PlotMarker(Group4_remap,"Nd145Di")
PlotMarker(Group4_remap,"Nd146Di")
PlotMarker(Group4_remap,"Nd148Di")
PlotMarker(Group4_remap,"Nd150Di")
PlotMarker(Group4_remap,"Pr141Di")
PlotMarker(Group4_remap,"Sm147Di")
PlotMarker(Group4_remap,"Sm149Di")
PlotMarker(Group4_remap,"Sm152Di")
PlotMarker(Group4_remap,"Sm154Di")
PlotMarker(Group4_remap,"Tb159Di")
PlotMarker(Group4_remap,"Tm169Di")
PlotMarker(Group4_remap,"Yb171Di")
PlotMarker(Group4_remap,"Yb172Di")
PlotMarker(Group4_remap,"Yb173Di")
PlotMarker(Group4_remap,"Yb174Di")
PlotMarker(Group4_remap,"Yb176Di")

#Plot group 5 markers
PlotStars(Group5_remap)

PlotMarker(Group5_remap,"Bi209Di")
PlotMarker(Group5_remap,"Dy161Di")
PlotMarker(Group5_remap,"Dy162Di")
PlotMarker(Group5_remap,"Dy163Di")
PlotMarker(Group5_remap,"Dy164Di")
PlotMarker(Group5_remap,"Er166Di")
PlotMarker(Group5_remap,"Er167Di")
PlotMarker(Group5_remap,"Er168Di")
PlotMarker(Group5_remap,"Er170Di")
PlotMarker(Group5_remap,"Eu151Di")
PlotMarker(Group5_remap,"Eu153Di")
PlotMarker(Group5_remap,"Gd155Di")
PlotMarker(Group5_remap,"Gd156Di")
PlotMarker(Group5_remap,"Gd160Di")
PlotMarker(Group5_remap,"Ho165Di")
PlotMarker(Group5_remap,"Lu175Di")
PlotMarker(Group5_remap,"Nd143Di")
PlotMarker(Group5_remap,"Nd145Di")
PlotMarker(Group5_remap,"Nd146Di")
PlotMarker(Group5_remap,"Nd148Di")
PlotMarker(Group5_remap,"Nd150Di")
PlotMarker(Group5_remap,"Pr141Di")
PlotMarker(Group5_remap,"Sm147Di")
PlotMarker(Group5_remap,"Sm149Di")
PlotMarker(Group5_remap,"Sm152Di")
PlotMarker(Group5_remap,"Sm154Di")
PlotMarker(Group5_remap,"Tb159Di")
PlotMarker(Group5_remap,"Tm169Di")
PlotMarker(Group5_remap,"Yb171Di")
PlotMarker(Group5_remap,"Yb172Di")
PlotMarker(Group5_remap,"Yb173Di")
PlotMarker(Group5_remap,"Yb174Di")
PlotMarker(Group5_remap,"Yb176Di")


###################################################################################################
#Attempt to change MST plot
#layout.reingold.tilford gives an interesting tree
BuildMST <- function(fsom, silent=FALSE, tSNE=FALSE){
  
  fsom$MST <- list()
  if(!silent) message("Building MST\n")
  
  adjacency <- stats::dist(fsom$map$codes, method = "euclidean")
  fullGraph <- igraph::graph.adjacency(as.matrix(adjacency), 
                                       mode = "undirected", 
                                       weighted = TRUE)
  fsom$MST$graph <- igraph::minimum.spanning.tree(fullGraph)
  fsom$MST$l <- igraph::layout.kamada.kawai(fsom$MST$graph)    
  
  if(tSNE){
    fsom$MST$l2 <- tsne(fsom$map$codes)   
    #library(RDRToolbox)
    #fsom$MST$l2 <- Isomap(fsom$map$codes,dims=2,k=3)[[1]]
  }
  
  UpdateNodeSize(fsom)
}

#Run FlowSOM ComData
set.seed(234)
out2 <- FlowSOM::ReadInput(ComFlowSOM, transform = FALSE, scale = FALSE)
out2 <- FlowSOM::BuildSOM(out2, colsToUse = ClusterCols)
out2 <- BuildMST(out2)

PlotStars(out2)


######################
#custom attempt
BuildMST <- function(fsom, silent=FALSE, tSNE=FALSE){
  
  fsom$MST <- list()
  if(!silent) message("Building MST\n")
  
  adjacency <- stats::dist(fsom$map$codes, method = "euclidean")
  fullGraph <- igraph::graph.adjacency(as.matrix(adjacency), 
                                       mode = "undirected", 
                                       weighted = TRUE)
  fsom$MST$graph <- igraph::minimum.spanning.tree(fullGraph)
  fsom$MST$l <- igraph::layout.kamada.kawai(fsom$MST$graph)    
  UpdateNodeSize(fsom)
}





set.seed(234)
out2 <- FlowSOM::ReadInput(ComFlowSOM, transform = FALSE, scale = FALSE)
out2 <- FlowSOM::BuildSOM(out2, colsToUse = ClusterCols)
out2 <- 
  
  PlotStars(out2)

###################################################################################################
