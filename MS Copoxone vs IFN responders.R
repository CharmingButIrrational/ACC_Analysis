#MS treatment dataset analysis (groups 1,3)
library(e1071)
library(gbm)
library(caret)
library(data.table)
library(pROC)
library(ROCR)
library(MASS)
library(SDMTools)

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/MS project/Data summary_für Bioinformatik.csv", sep = ";",dec = ",", header = F,)
#Transpose data so that outcome is a column
df_transpose <- transpose(mydata)
#Covert first row to column names
colnames(df_transpose) <- as.character(unlist(df_transpose[1,]))
#Remove names and extra columns
df_subset <- df_transpose[-c(1,131,132),]
#Rename the status column
names(df_subset)[1]<-paste("Therapy")
#Convert variables to numeric
df_subset[] <- lapply(df_subset, function(x) as.numeric(as.character(x)))
#Convert target to factor
df_subset$Therapy <- as.factor(df_subset$Therapy)

#Check the data for structure and duplicate colnames 
str(df_subset)
any(duplicated(names(df_subset)))
length(which(duplicated(names(df_subset))==TRUE))
#P_AQP4 141-160 col: 60 and 258
identical(df_subset[,c(60)],df_subset[,c(258)])
#Rename the columns to differenciate them
names(df_subset)[60]<-paste("P_AQP4 141-160_A")
names(df_subset)[258]<-paste("P_AQP4 141-160_B")
#Recheck columns
any(duplicated(names(df_subset)))

#Split the dataset 
Therapy1.3 <- split(df_subset, df_subset$Therapy) 
#Create df for Copaxone reponders and non-responders
response <- rbind(Therapy1.3$`1`, Therapy1.2$'3')
response$Therapy <- factor(response$Therapy)
#Split data into testing and training sets
smp_size <- floor(0.8 * nrow(response))
set.seed(345)
train_ind <- sample(seq_len(nrow(response)), size = smp_size)

train <- response[train_ind, ]
test <- response[-train_ind, ]

#Drop factor levels which don't occur
train$Therapy <- factor(train$Therapy)
test$Therapy <- factor(test$Therapy)

#Elastic net 
train.elastic <- trainControl(method = "repeatedcv", number = 10, repeats = 5, savePredictions = TRUE, search = "random")
set.seed(233)
elastic <- train(Therapy ~ ., 
                 data = train, 
                 method = "glmnet",
                 trControl = train.elastic,
                 preProc = c("center", "scale"),
                 tuneLength = 10)
importance.elastic <- varImp(elastic, scale = FALSE)
print(importance.elastic)
plot(importance.elastic, top = 12)

elastic.pred <- predict(elastic, test)

confusion.elastic <- confusionMatrix(elastic.pred, test$Therapy, positive = "1")
confusion.elastic

elastic.outcome <- test$Therapy

#Create the ROC curve 
pred.elastic <- prediction(as.numeric(elastic.pred), as.numeric(elastic.outcome))
perf.elastic <- performance(pred.elastic, "tpr", "fpr")
plot(perf.elastic, col = "red", lwd = 2)
abline(0, 1, lty = 2)

#Create AUC value
auc.perf.elastic = performance(pred.elastic, measure = "auc")
auc.perf.elastic@y.values


#SVM
train.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5, savePredictions = TRUE, search = "random")
set.seed(233)
svm <- train(Therapy ~ ., 
             data = train, 
             method = "svmLinear",
             trControl = train.svm,
             preProc = c("center", "scale"),
             tuneLength = 10,
             allowParallel = T)
importance.svm <- varImp(svm, scale = FALSE)
print(importance.svm)
plot(importance.svm, top = 20)

svm.pred <- predict(svm, test)

confusion.svm <- confusionMatrix(svm.pred, test$Therapy)
confusion.svm

svm.outcome <- test$Therapy

#Create the ROC curve
pred.svm <- prediction(as.numeric(svm.pred), as.numeric(svm.outcome))
perf.svm <- performance(pred.svm, "tpr", "fpr")
plot(perf.svm, col = "red", lwd = 2)
abline(0, 1, lty = 2)

#Create AUC value 
auc.perf.svm = performance(pred.svm, measure = "auc")
auc.perf.svm@y.values


#GBM
train.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5, savePredictions = TRUE)
caretGrid <- expand.grid(interaction.depth=c(1, 3, 5), n.trees = (0:50)*50,
                         shrinkage=c(0.01, 0.001),
                         n.minobsinnode=10)
set.seed(233)
gbm <- train(Therapy ~ ., 
             data = train, 
             method = "gbm",
             trControl = train.gbm,
             preProc = c("center", "scale"),
             tuneGrid = caretGrid,
             tuneLength = 10)
importance.gbm <- varImp(gbm, scale = T)  
print(importance.gbm)
plot(importance.gbm, top = 9)

gbm.pred <- predict(gbm, test)

confusion.gbm <- confusionMatrix(gbm.pred, test$Therapy)
confusion.gbm

gbm.outcome <- test$Therapy

#Create the ROC curve
pred.gbm <- prediction(as.numeric(gbm.pred), as.numeric(gbm.outcome))
perf.gbm <- performance(pred.gbm, "tpr", "fpr")
plot(perf.gbm, col = "red", lwd = 2)
abline(0, 1, lty = 2)

#Create AUC value 
auc.perf.gbm = performance(pred.gbm, measure = "auc")
auc.perf.gbm@y.values

################################################################
#Compare models
results <- resamples(list(ELA=elastic, SVM=svm, GBM=gbm))
summary(results)
bwplot(results)

################################################################

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
elastic.predictors.Sig <- head(elastic.predictors, 12)
SVM.predictors.Sig <- head(SVM.predictors, 20)
#Adjust GBM selection based on printed variables
GBM.predictors.Sig <- head(GBM.predictors, 20)
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
print(venn.intersect.predictors$a4) 
