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
df_alt.num <- df_alt
#convert target to factor
df_alt$Samples <- as.factor(df_alt$Samples)

df_alt.num$Samples <- as.numeric(as.character(df_alt.num$Samples))

###################################################################
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
                          na.action = na.pass,
                          allowParallel = TRUE)
importance.svm <- varImp(svm, scale = FALSE)
print(importance.svm)
plot(importance.svm, top = 20)

#GBM 
train.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 6, 9, 10),
                        n.trees = (0:50)*50, 
                        shrinkage = seq(.0005, .05,.0005),
                        n.minobsinnode = 10)
set.seed(86)
gbm <- train(Samples ~ .,
                          data = df_alt
                          , method = "gbm",
                          trControl = train.gbm,
                          preProcess = c("center", "scale"),
                          tuneLength = 10,
                          tuneGrid = gbmGrid,
                          na.action = na.pass,
                          allowParallel = TRUE)
importance.gbm <- varImp(gbm, scale = FALSE)
print(importance.gbm)
plot(importance.gbm)

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

#############################################################
set.seed(233)
elastic.ms <- train(Samples ~ ., 
  data = df_alt, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneLength = 10)
print(elastic.ms)
print(elastic.ms$bestTune)
coef(elastic.ms)


#Correlation matrix
library(caret)
library(corrplot)
set.seed(767)
corr.sig <- cor(df_alt.num, df_alt.num, method = "spearman")
corr.sig
corrplot(corr.sig, type = "upper", method="number", is.corr=FALSE)

library(randomForest)
library(mlbench)
set.seed(7658)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
results <- rfe(df_alt, df_alt, rfeControl=control)
print(results)
predictors(results)
plot(results, type=c("g", "o"))

######################################################

library(parallel)
library(pvclust)
#This cluster method works on columns, not rows, so it fits with the transposed data I used for LVQ and SVM
#Use parallel package to support parallel computing
cores <- detectCores() - 1 
cl <- makeCluster(cores)
cluster.fit <- parPvclust(cl, df_alt[,-1], method.hclust="ward.D",
               method.dist="euclidean", nboot = 10000)
plot(cluster.fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(cluster.fit, alpha=.95) 


cluster.fit.2 <- parPvclust(cl, df_alt[,-1], method.hclust="average",
                          method.dist="euclidean", nboot = 10000)
plot(cluster.fit.2) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(cluster.fit.2, alpha=.99) 

#############################################################

sam.data <- mydata

sam.data.red <- sam.data[,-1]

Outcome <- as.character(as.vector(sam.data[1,]))
Outcome[1] <- NULL


library(samr)

sam.data <- SAM(sam.data.red,censoring.status=NULL,
            resp.type=c("Two class unpaired"),
            s0=NULL, 
            s0.perc=NULL, 
            nperms=10000, 
            center.arrays=TRUE, 
            testStatistic=c("standard"), 
            time.summary.type=c("slope"), 
            regression.method=c("standard"), 
            knn.neighbors=10, 
            random.seed=101332,
            logged2 = FALSE)

########################################################################
df_alt.num.mat <- as.matrix(df_alt.num)

heatmap(df_alt.num.mat)

