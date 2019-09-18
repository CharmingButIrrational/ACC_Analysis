#Pheochromocytoma analysis

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/Phenochromacytoma/SPSS.csv", sep = ";", header = T)

class <- as.factor(mydata$class)

alt.data <- mydata[,21:150]

alt.data <- cbind(class, alt.data)

library(gbm)
library(glmnet)
library(caret)
library(pROC)
library(ROCR)
library(MASS)
library(SDMTools)


smp_size <- floor(0.8 * nrow(alt.data))
set.seed(345)
train_ind <- sample(seq_len(nrow(alt.data)), size = smp_size)

train <- alt.data[train_ind, ]
test <- alt.data[-train_ind, ]

############################################################################
#Generate ROC from whole dataset
#Elastic net using caret
train.elastic.whole <- trainControl(method = "repeatedcv", number = 10, repeats = 5, savePredictions = TRUE, search = "random")
set.seed(233)
elastic.spss.whole <- train(class ~ ., 
                      data = alt.data, 
                      method = "glmnet",
                      trControl = train.elastic.whole,
                      preProc = c("center", "scale"),
                      tuneLength = 10)
importance.elastic <- varImp(elastic.spss.whole, scale = FALSE)
print(importance.elastic)
plot(importance.elastic, top = 20)



############################################################################################

#Elastic net using caret
train.elastic <- trainControl(method = "repeatedcv", number = 10, repeats = 5, savePredictions = TRUE, search = "random")
set.seed(233)
elastic.spss <- train(class ~ ., 
                      data = train, 
                      method = "glmnet",
                      trControl = train.elastic,
                      preProc = c("center", "scale"),
                      tuneLength = 10)
importance.elastic <- varImp(elastic.spss, scale = FALSE)
print(importance.elastic)
plot(importance.elastic, top = 20)

elastic.pred <- predict(elastic.spss, test)

confusion.elastic <- confusionMatrix(elastic.pred, test$class, positive = "1")
confusion.elastic

elastic.outcome <- test$class

### create the ROC curve ###
pred.elastic <- prediction(as.numeric(elastic.pred), as.numeric(elastic.outcome))
perf.elastic <- performance(pred.elastic, "tpr", "fpr")
plot(perf.elastic, col = "red", lwd = 2)
abline(0, 1, lty = 2)

### create AUC value ###
auc.perf.elastic = performance(pred.elastic, measure = "auc")
auc.perf.elastic@y.values


#SVM
train.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5, savePredictions = TRUE, search = "random")
set.seed(233)
svm.spss <- train(class ~ ., 
                  data = train, 
                  method = "svmLinear",
                  trControl = train.svm,
                  preProc = c("center", "scale"),
                  tuneLength = 10,
                  allowParallel = T)
importance.svm <- varImp(svm.spss, scale = FALSE)
print(importance.svm)
plot(importance.svm, top = 20)

svm.pred <- predict(svm.spss, test)

confusion.svm <- confusionMatrix(svm.pred, test$class)
confusion.svm

svm.outcome <- test$class

### create the ROC curve ###
pred.svm <- prediction(as.numeric(svm.pred), as.numeric(svm.outcome))
perf.svm <- performance(pred.svm, "tpr", "fpr")
plot(perf.svm, col = "red", lwd = 2)
abline(0, 1, lty = 2)

### create AUC value ###
auc.perf.svm = performance(pred.svm, measure = "auc")
auc.perf.svm@y.values


#GBM
train.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5, savePredictions = TRUE)
caretGrid <- expand.grid(interaction.depth=c(1, 3, 5), n.trees = (0:50)*50,
                         shrinkage=c(0.01, 0.001),
                         n.minobsinnode=10)
set.seed(233)
gbm.spss <- train(class ~ ., 
                  data = train, 
                  method = "gbm",
                  trControl = train.gbm,
                  preProc = c("center", "scale"),
                  tuneGrid = caretGrid,
                  tuneLength = 10)
importance.gbm <- varImp(gbm.spss, scale = T)   #9 variables of importance
print(importance.gbm)
plot(importance.gbm, top = 9)

gbm.pred <- predict(gbm.spss, test)

confusion.gbm <- confusionMatrix(gbm.pred, test$class)
confusion.gbm

gbm.outcome <- test$class

### create the ROC curve ###
pred.gbm <- prediction(as.numeric(gbm.pred), as.numeric(gbm.outcome))
perf.gbm <- performance(pred.gbm, "tpr", "fpr")
plot(perf.gbm, col = "red", lwd = 2)
abline(0, 1, lty = 2)

### create AUC value ###
auc.perf.gbm = performance(pred.gbm, measure = "auc")
auc.perf.gbm@y.values


##################################################################
#Compare models
results <- resamples(list(ELA=elastic.spss, GBM=gbm.spss, SVM=svm.spss))
summary(results)
bwplot(results)

#In this dataset GBM has greater accuracy and kappa values
##################################################################

library(data.table)

#Common predictors
elast.predictors <- importance.elastic$importance
SVM.predictors <- importance.svm$importance
GBM.predictors <- importance.gbm$importance
#Remove dupicate column
all(SVM.predictors$X0 == SVM.predictors$X1)
SVM.predictors$X1 <- NULL
colnames(SVM.predictors)[1] <- "Overall"
elastic.predictors <- setDT(elast.predictors, keep.rownames = TRUE)[]
SVM.predictors <- setDT(SVM.predictors, keep.rownames = TRUE)[]
GBM.predictors <- setDT(GBM.predictors, keep.rownames = TRUE)[]
elastic.predictors <- elastic.predictors[order(-elastic.predictors$Overall),]
SVM.predictors <- SVM.predictors[order(-SVM.predictors$Overall),]
GBM.predictors <- GBM.predictors[order(-GBM.predictors$Overall),]
elastic.predictors.Sig <- head(elastic.predictors, 20)
SVM.predictors.Sig <- head(SVM.predictors, 20)
#Adjust GBM selection based on printed variables
GBM.predictors.Sig <- head(GBM.predictors, 9)
elastic.predictors.Sig$Overall <- NULL
SVM.predictors.Sig$Overall <- NULL
GBM.predictors.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
elastic.predictors.Sig <- unlist(elastic.predictors.Sig)
SVM.predictors.Sig <- unlist(SVM.predictors.Sig)
GBM.predictors.Sig <- unlist(GBM.predictors.Sig)

library(VennDiagram)
venn.data.predictors <- list(elastic.predictors.Sig, SVM.predictors.Sig, GBM.predictors.Sig)
grid.newpage()
venn.plot.predictors <- venn.diagram(x = list(elastic.predictors.Sig=elastic.predictors.Sig, SVM.predictors.Sig=SVM.predictors.Sig, GBM.predictors.Sig=GBM.predictors.Sig),
                                     filename=NULL, 
                                     fill = c("red", "blue", "green"),
                                     alpha = 0.50,
                                     col = "transparent")
grid.draw(venn.plot.predictors)
venn.intersect.predictors <- calculate.overlap(venn.data.predictors)
print(venn.intersect.predictors$a5) #Variables selected by all models

###################################################################
