#FlowMap package development script

library(flowMap)

#Import data Group 1
group1 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 1.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))

#Group 1 columns 
marker_cols_G1 <- c(1,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18,22,23,24,25,27,28,29,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69)
plot_cols_G1 <- c(11,12,13,14,17,18,22,25,27,32,35,37,39,40,49,56,59,61,62,65,66,67,68,69)
phenotype_G1 <- c(10,11,12,13,14,15,16,17,18,22,23,25,27,32,35,37,38,39,40,49,56,57,58,59,61,62,65,66,67,68,69)

#Transform the data
asinh_scale <- 5
group1[, marker_cols_G1] <- asinh(group1[, marker_cols_G1] / asinh_scale)

#Create flowFrame object
group1_FlowSOM <- flowCore::flowFrame(group1)

#Run FlowSOM
set.seed(1234)
out_G1 <- FlowSOM::ReadInput(group1_FlowSOM, transform = FALSE, scale = FALSE)
out_G1 <- FlowSOM::BuildSOM(out_G1, colsToUse = plot_cols_G1)
out_G1 <- FlowSOM::BuildMST(out_G1)

#Import data Group 2
group2 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 2.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))

#Group 2 columns 
marker_cols_G2 <- c(1,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,27,28,29,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69)
plot_cols_G2 <- c(11,12,13,14,17,18,22,25,27,32,35,37,39,40,49,56,59,61,62,65,66,67,68,69)
phenotype_G2 <- c(10,11,12,13,14,15,16,17,18,22,23,25,27,32,35,37,38,39,40,49,56,57,58,59,61,62,65,66,67,68,69)

#Transform the data
group2[, marker_cols_G2] <- asinh(group2[, marker_cols_G2] / asinh_scale)

#Create flowFrame object
group2_FlowSOM <- flowCore::flowFrame(group2)

#Run FlowSOM
set.seed(1234)
out_G2 <- FlowSOM::ReadInput(group2_FlowSOM, transform = FALSE, scale = FALSE)
out_G2 <- FlowSOM::BuildSOM(out_G2, colsToUse = plot_cols_G2)
out_G2 <- FlowSOM::BuildMST(out_G2)


###########################################################
################  Compare cell populations  ###############
###########################################################

#Select samples from groups 
#Pick the group ids for comparison 
Sam1 <- out_G1[out_G1$id==1,]
Sam2 <- out_G2[out_G2$id==1,]

#Combine the samples
matrix(data, nrow = rows, ncol = cols) <- rbind(Sam1, Sam2)

#Select number of events to sample
sampleSize = 100

#Sample events from the two populations+
nn1 = round(sampleSize*table(mat$id)[1]/nrow(mat))
nn2 = round(sampleSize*table(mat$id)[2]/nrow(mat))
submat = rbind(mat1[sample(nrow(mat1),nn1),],mat2[sample(nrow(mat2),nn2),])
colnames(submat)[5] = "sam"

#Plot MST of the 100 events
g1 = makeFRMST(submat)
par(mar=c(0,0,0,0))
plot(g1$g,vertex.label.cex=0.01,
     layout=layout.fruchterman.reingold(g1$g))

#############################################################
########### Mapping cell population across FCMs #############
#############################################################

#Generate a matrix of Friedman-Rafsky statistics
res1 = getFRest(out_G1,out_G2,sampleMethod="proportional",
                sampleSize=100,
                ndraws=100,estStat="median",
                ncores=NULL)

library(gplots)
par(mar=c(0,0,0,0))
heatmapCols <- colorRampPalette(c("red","yellow","white","blue"))(50)
#Plot the F-R statistics comparing Sample 1 & 2
heatmap.2(res1@ww,trace="none",col=heatmapCols,symm=FALSE,dendrogram="none",
          Rowv=FALSE,Colv=FALSE,xlab="Sample 2",ylab="Sample 1")
#Alternative heatmap of p-values
heatmap.2(res1@pNorm,trace="none",col=heatmapCols,symm=FALSE,dendrogram="none",
          Rowv=FALSE,Colv=FALSE,xlab="Sample 2",ylab="Sample 1")