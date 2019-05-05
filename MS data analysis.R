library(e1071)
library(caret)
library(data.table)

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/MS data sheet 1.csv", sep = ";", header = F,)
#Transpose data so that outcome is a column
df_transpose <- transpose(mydata)
#Covert first row to column names
colnames(df_transpose) <- as.character(unlist(df_transpose[1,]))#
#remove gene names(symbols causing errors) and sample with no outcome
df_subset <- df_transpose[-c(1,2),]
#Add a name to the status column
names(df_subset)[1]<-paste("Samples")
setwd("C:/Users/oisin/Desktop/Analysis Data/")

write.csv(df_subset, file = "Altered MS data")

df_alt <- read.csv("C:/Users/oisin/Desktop/Analysis Data/Altered MS data")
#Remove X column at the start of the df
df_alt$X <- NULL
#convert target to factor
df_alt$Samples <- as.factor(df_alt$Samples)

#LVQ 
train.lvq <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
lvq <- train(Samples ~ .,
                          data = df_alt, method = "lvq",
                          trControl = train.lvq,
                          preProcess = c("center", "scale"),
                          tuneGrid = data.frame(size = 36, k = 1:38),
                          tuneLength = 10,
                          na.action = na.pass)
importance.lvq <- varImp(lvq, scale = FALSE)
print(importance.lvq)
plot(importance.lvq, top = 20)
#SVM 
train.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
svm <- train(Samples ~ .,
                          data = df_alt, method = "svmLinear",
                          trControl = train.svm,
                          preProcess = c("center", "scale"),
                          tuneLength = 10,
                          na.action = na.pass)
importance.svm <- varImp(svm, scale = FALSE)
print(importance.svm)
plot(importance.svm, top = 20)

#GBM 
#train.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#gbm <- train(Samples ~ .,
#                          data = df_alt
#                          , method = "gbm",
#                          trControl = train.gbm,
#                          preProcess = c("center", "scale"),
#                          tuneGrid = expand.grid(interaction.depth=seq(1,6,by=1),
#                                                 n.trees=c(25,50,100,200),
#                                                 shrinkage=c(1e-3),
#n.minobsinnode = 10),
#                          tuneLength = 10,
#                          na.action = na.pass)
#importance.gbm <- varImp(gbm, scale = FALSE)
#print(importance.gbm)
#plot(importance.gbm)

#Common predictors for  progress..yes.no.
LVQ.predictors <- importance.lvq$importance
SVM.predictors <- importance.svm$importance
#GBM.predictors <- importance.gbm$importance
all(LVQ.predictors$X0 == LVQ.predictors$X1)
all(SVM.predictors$X0 == SVM.predictors$X1)
LVQ.predictors$X1 <- NULL
SVM.predictors$X1 <- NULL
colnames(LVQ.predictors)[1] <- "Overall"
colnames(SVM.predictors)[1] <- "Overall"
LVQ.predictors <- setDT(LVQ.predictors, keep.rownames = TRUE)[]
SVM.predictors <- setDT(SVM.predictors, keep.rownames = TRUE)[]
#GBM.predictors <- setDT(GBM.predictors, keep.rownames = TRUE)[]
LVQ.predictors <- LVQ.predictors[order(-LVQ.predictors$Overall),]
SVM.predictors <- SVM.predictors[order(-SVM.predictors$Overall),]
#GBM.predictors <- GBM.predictors[order(-GBM.predictors$Overall),]
LVQ.predictors.Sig <- head(LVQ.predictors, 20)
SVM.predictors.Sig <- head(SVM.predictors, 20)
#Adjust GBM selection based on printed variables
#GBM.predictors.Sig <- head(GBM.predictors, 9)
LVQ.predictors.Sig$Overall <- NULL
SVM.predictors.Sig$Overall <- NULL
#GBM.predictors.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.predictors.Sig <- unlist(LVQ.predictors.Sig)
SVM.predictors.Sig <- unlist(SVM.predictors.Sig)
#GBM.predictors.Sig <- unlist(GBM.predictors.Sig)
library(VennDiagram)
venn.data.predictors <- list(LVQ.predictors.Sig, SVM.predictors.Sig)
grid.newpage()
venn.plot.predictors <- venn.diagram(x = list(LVQ.predictors.Sig=LVQ.predictors.Sig, SVM.predictors.Sig=SVM.predictors.Sig),
                                       filename=NULL, 
                                       fill = c("red", "blue"),
                                       alpha = 0.50,
                                       col = "transparent")
grid.draw(venn.plot.predictors)
venn.intersect.predictors <- calculate.overlap(venn.data.predictors)
print(venn.intersect.predictors$a3)



library(parallel)
library(pvclust)
#This cluster method works on columns, not rows, so it fits with the transposed data I used for LVQ and SVM
#Use parallel package to support parallel computing
cores <- detectCores() - 1 
cl <- makeCluster(cores)
cluster.fit <- parPvclust(cl, df_alt, method.hclust="ward.D",
               method.dist="euclidean", )
plot(cluster.fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(cluster.fit, alpha=.95) 



library(samr)

