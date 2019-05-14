#Pheochromocytoma analysis

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/SPSS.csv", sep = ";", header = T)

class <- as.factor(mydata$class)

alt.data <- mydata[,21:150]

alt.data <- cbind(class, alt.data)

#Elastic net using caret
train.elastic <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(233)
elastic.spss <- train(class ~ ., 
                      data = alt.data, 
                      method = "glmnet",
                      trControl = train.elastic,
                      preProc = c("center", "scale"),
                      tuneLength = 10)
print(elastic.spss)
importance.elastic <- varImp(elastic.spss, scale = FALSE)
print(importance.elastic)
plot(importance.elastic, top = 20)


#SVM
train.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(233)
svm.spss <- train(class ~ ., 
                      data = alt.data, 
                      method = "svmLinear",
                      trControl = train.svm,
                      preProc = c("center", "scale"),
                      tuneLength = 10)
importance.svm <- varImp(svm.spss, scale = FALSE)
print(importance.svm)
plot(importance.svm, top = 20)


#GBM
train.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(233)
gbm.spss <- train(class ~ ., 
                  data = alt.data, 
                  method = "gbm",
                  trControl = train.gbm,
                  preProc = c("center", "scale"),
                  tuneLength = 10)
importance.gbm <- varImp(gbm.spss, scale = FALSE)
print(importance.gbm)
plot(importance.gbm, top = 20)
