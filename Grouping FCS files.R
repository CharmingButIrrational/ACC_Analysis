#Aggregate flow cytometry files
library(flowCore)
library(FlowSOM)
library(matrixStats)

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
group5 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 5.fcs",
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

PlotStars(Group1_remap)
PlotStars(Group2_remap)
PlotStars(Group3_remap)
PlotStars(Group4_remap)
PlotStars(Group5_remap)
