library(flowCore)
library(FlowSOM)
library(ggplot2)
library(Rtsne)
library(scatterplot3d)

#Load data
group1 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 1.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))

#select marker columns to use for clustering
marker_cols <- c(1,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,27,28,29,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69)

#Transform data
asinh_scale <- 5
group1[, marker_cols] <- asinh(group1[, marker_cols] / asinh_scale)

#create flowFrame object 
group1_FlowSOM <- flowCore::flowFrame(group1)

plot_cols <- c(11,12,13,14,17,18,22,25,27,32,35,37,39,40,49,56,59,61,62,65,66,67,68,69)
#run FlowSOM
set.seed(1234)

#initial step prior to meta-clustering
out <- FlowSOM::ReadInput(group1_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = plot_cols)
out <- FlowSOM::BuildMST(out)

#BuildingMST using selected markers
example <- FlowSOM::ReadInput(group1_FlowSOM, transform = FALSE, scale = FALSE)
example <- FlowSOM::BuildSOM(example, colsToUse = marker_cols)
example <- FlowSOM::BuildMST(example)

#Number the nodes
PlotNumbers(UpdateNodeSize(out,reset=TRUE))

#Plot the SOM
PlotStars(out)  #May be difficult to view the central area

PlotStars(example)

#Alternate methods of visualisation 
PlotStars(out, view = "grid")
#PlotNumbers(UpdateNodeSize(out,reset = TRUE), view = "grid")
PlotStars(out, view = "tSNE")

#Investigation of a particular marker
#Print the markers
print(colnames(out$map$medianValues))

#Plot the different markers
PlotMarker(out,"Dy162Di")
PlotMarker(out,"Dy163Di")
PlotMarker(out,"Dy164Di")
PlotMarker(out,"Er166Di")
PlotMarker(out,"Er170Di")
PlotMarker(out,"Eu151Di")
PlotMarker(out,"Gd155Di")
PlotMarker(out,"Gd160Di")
PlotMarker(out,"Ho165Di")
PlotMarker(out,"Lu175Di")
PlotMarker(out,"Nd143Di")
PlotMarker(out,"Nd145Di")
PlotMarker(out,"Nd148Di")
PlotMarker(out,"Nd150Di")
PlotMarker(out,"Pr141Di")
PlotMarker(out,"Sm147Di")
PlotMarker(out,"Sm154Di")
PlotMarker(out,"Tb159Di")
PlotMarker(out,"Tm169Di")
PlotMarker(out,"Yb171Di")
PlotMarker(out,"Yb172Di")
PlotMarker(out,"Yb173Di")
PlotMarker(out,"Yb174Di")
PlotMarker(out,"Yb176Di")

#Investigate a particular node(s)
PlotNumbers(UpdateNodeSize(out,reset=TRUE))
#Select clusters and markers e.g Podoplanin and Drp1 
PlotClusters2D(out,"Gd160Di","Yb171Di",c(99,1,34))

#Extract cellular labels before meta-clusteringt
labels_pre <- out$map$mapping[, 1]

#Running meta clustering
k <- 7  #Expected number of classes
seed <- 46456
meta <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = k, seed = seed)
meta <- meta[[k]]$consensusClass


#Cluster labels from output object
labels <- meta[labels_pre]

#Save cluster labels
res <- data.frame(cluster = labels)
write.table(res, file = "C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group1_cluster_labels_FlowSOM.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

#Subsampling 
n_sub <- 2000
set.seed(1234)
ix <- sample(1:length(labels), n_sub)

#Prepare data for Rtsne (matrix format)
data_Rtsne <- group1[ix, plot_cols]
data_Rtsne <- as.matrix(data_Rtsne)

#Remove any near-duplicate rows (required by Rtsne)
dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]

#Run Rtsne (Barnes-Hut-SNE algorithm)
set.seed(6784)
out_Rtsne <- Rtsne(data_Rtsne,  dims = 3, perplexity = 10, pca = FALSE, verbose = TRUE)
                                        #Perplexity generally is between 5-50

#Load cluster label
file_labels <- "C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group1_cluster_labels_FlowSOM.txt"
data_labels <- read.table(file_labels, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
labels <- data_labels[, "cluster"]

#Select points used by Rtsne
labels_plot <- labels[ix][!dups]

#Prepare Rtsne output data for plot
data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")
data_plot[, "cluster"] <- as.factor(labels_plot)

#Plot 2D t-SNE projection
ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 2) + 
  coord_fixed(ratio = 1) + 
  ggtitle("t-SNE projection with FlowSOM clustering") + 
  theme_bw()

#Plot 3D t-SNE 
scatterplot3d(x=out_Rtsne$Y[,1],y=out_Rtsne$Y[,2],z=out_Rtsne$Y[,3],
              color = labels_plot)



########################################################################
#CITRUS clustering

library("citrus")
#Run launcher, navigate to folder containing fcs files 
citrus.launchUI()

#Load data
group1 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 1.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Load data
group2 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 2.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Load data
group3 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 3.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Load data
group4 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 4.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
#Load data
group5 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 5.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
dim(group1)
dim(group2)
dim(group3)
dim(group4)
dim(group5)

vec.1 <- colnames(group1)
vec.2 <- colnames(group2)
vec.3 <- colnames(group3)
vec.4 <- colnames(group4)
vec.5 <- colnames(group5)

Reduce(intersect, list(vec.1, vec.2, vec.3, vec.4, vec.5))

