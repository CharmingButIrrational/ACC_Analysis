library(flowCore)
library(FlowSOM)
library(ggplot2)
library(Rtsne)

#import data

group1 <- flowCore::exprs(flowCore::read.FCS("C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group 1.fcs",
                                             transformation = FALSE,
                                             truncate_max_range = FALSE))
head(group1)
dim(group1)

#select marker columns to use for clustering

marker_cols <- c(1,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,27,28,29,30,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,48,49,50,51,52,54,56,57,58,59,60,61,62,63,64,65,66,67,68,69)

plot_cols <- c(11,12,13,14,17,18,22,25,27,32,35,37,39,40,49,56,59,61,62,65,66,67,68,69)

#apply arcsinh transformation
#with standard factor 5 for cyTOF data

asinh_scale <- 5
group1[, marker_cols] <- asinh(group1[, marker_cols] / asinh_scale)

summary(group1)

#create flowFrame object (required for FlowSOM) from data matrix

group1_FlowSOM <- flowCore::flowFrame(group1)

#run FlowSOM
#set seed for reproducibility
set.seed(1234)

#initial step prior to meta-clustering

out <- FlowSOM::ReadInput(group1_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = plot_cols)
out <- FlowSOM::BuildMST(out)

#Visualisation
FlowSOM::PlotStars(out)
PlotStars(out, view = "tSNE")

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

#extract cellular labels (pre meta-clustering) from output object

labels_pre <- out$map$mapping[, 1]

#run meta clustering
k <- 20
seed <- 4356456
                    
out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = k, seed = seed)
out <- out[[k]]$consensusClass


#extract cluster labels from output object

labels <- out[labels_pre]


#summary of cluster sizes and number of clusters

table(labels)
length(table(labels))

# save cluster labels

res <- data.frame(cluster = labels)

write.table(res, file = "C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group1_cluster_labels_FlowSOM.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")


# subsampling (required due to runtime)

n_sub <- 2000

set.seed(1234)
ix <- sample(1:length(labels), n_sub)

# prepare data for Rtsne (matrix format required)

data_Rtsne <- group1[ix, marker_cols]
data_Rtsne <- as.matrix(data_Rtsne)

head(data_Rtsne)
dim(data_Rtsne)

# remove any near-duplicate rows (required by Rtsne)

dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]

dim(data_Rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm)

# note initial PCA is not required, since we do not have too many dimensions
# (i.e. not thousands, which may be the case in other domains)

set.seed(1234)
out_Rtsne <- Rtsne(data_Rtsne,  dims = 3, perplexity = 40, pca = FALSE, verbose = TRUE)
                                         

# load cluster labels (if not still loaded)

file_labels <- "C:/Users/oisin/Desktop/Analysis Data/cyTOF/CyTOF samples/Group1_cluster_labels_FlowSOM.txt"
data_labels <- read.table(file_labels, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
labels <- data_labels[, "cluster"]

# select points used by Rtsne

labels_plot <- labels[ix][!dups]
length(labels_plot)  ## should be same as number of rows in data_Rtsne

# prepare Rtsne output data for plot

data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

head(data_plot)
dim(data_plot)  ## should match length of labels_plot (otherwise labels will not match up correctly)

data_plot[, "cluster"] <- as.factor(labels_plot)

head(data_plot)

# plot 2-dimensional t-SNE projection

ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 2) + 
  coord_fixed(ratio = 1) + 
  ggtitle("t-SNE projection with FlowSOM clustering") + 
  theme_bw()




