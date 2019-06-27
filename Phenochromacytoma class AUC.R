#Pheochromocytoma analysis

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/SPSS.csv", sep = ";", header = T)

class <- as.factor(mydata$class)

alt.data <- mydata[,21:150]

alt.data <- cbind(class, alt.data)

library(gbm)
library(glmnet)
library(caret)
library(pROC)
library(ROCR)


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

#ROC plot for elastic net
selectedglmnet <- elastic.spss.whole$pred$pred == 2
# Plot:
plot.roc(elastic.spss.whole$pred[selectedglmnet],
         elastic.spss.whole$pred[selectedglmnet])

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

confusion.elastic <- confusionMatrix(elastic.pred, test$class)
confusion.elastic
#Plotting ROC graph
ROC.elastic <- performance(elastic.pred, "tpr", "fpr")
plot(ROC.elastic)
plot(ROC.elastic, add = TRUE, col = "green")
legend("right", legend = c("rf"), col = c("green"), lty = 1:2, cex = 0.6)

#Calculating AUC
AUC.elastic <- performance(elastic.pred, "AUC")
AUC.elastic <- as.numeric(AUC.elastic@y.values)
AUC.elastic


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