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

#Labelling each of the samples 
group1_sub <- cbind(group1_sub, Sample = 1 )
group2_sub <- cbind(group2_sub, Sample = 2 )
group3_sub <- cbind(group3_sub, Sample = 3 )
group4_sub <- cbind(group4_sub, Sample = 4 )
group5_sub <- cbind(group5_sub, Sample = 5 )

dim(group1_sub)
dim(group2_sub)
dim(group3_sub)
dim(group4_sub)
dim(group5_sub)

#Checking all the samples have the same columns
setdiff(colnames(group1_sub), colnames(group2_sub))
setdiff(colnames(group1_sub), colnames(group3_sub))
setdiff(colnames(group1_sub), colnames(group4_sub))
setdiff(colnames(group1_sub), colnames(group5_sub))


