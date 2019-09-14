#MS treatment dataset analysis 
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
Therapy1.2 <- split(df_subset, df_subset$Therapy) 
#Create df for Copaxone reponders and non-responders
response <- rbind(Therapy1.2$`1`, Therapy1.2$`2`)
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
print(venn.intersect.predictors$a4) #Intersect of GBM and elastic

#############################################################

Antigen <- response$Antigen

`A5024 hIgGM` <- response$`A5024 hIgGM`                    
`LS3 AQP4 11-30` <- response$`LS3 AQP4 11-30`             
`LS4 AQP4 21-40` <- response$`LS4 AQP4 21-40`                
`14 KIR4.1 131-150` <- response$`14 KIR4.1 131-150`           
`A4874 MBP 31-50` <- response$`A4874 MBP 31-50`           
`A4754 hAQP4 - Mayo in 0.3HAPS` <- response$`A4754 hAQP4 - Mayo in 0.3HAPS`
`LS20 AQP4 181-200` <- response$`LS20 AQP4 181-200`

var1 <- as.data.frame(cbind(Antigen, `A5024 hIgGM`))
var2 <- as.data.frame(cbind(Antigen, `LS3 AQP4 11-30`))
var3 <- as.data.frame(cbind(Antigen, `LS4 AQP4 21-40`))
var4 <- as.data.frame(cbind(Antigen, `14 KIR4.1 131-150`))
var5 <- as.data.frame(cbind(Antigen, `A4874 MBP 31-50`))
var6 <- as.data.frame(cbind(Antigen, `A4754 hAQP4 - Mayo in 0.3HAPS`))
var7 <- as.data.frame(cbind(Antigen, `LS20 AQP4 181-200`))




#`A5024 hIgGM`  
group_by(var1, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A5024 hIgGM`, na.rm = TRUE),
            median = median(`A5024 hIgGM`, na.rm = TRUE),
            sd = sd(`A5024 hIgGM`, na.rm = TRUE))

ggboxplot(var1, x = "Antigen", y = "`A5024 hIgGM`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A5024 hIgGM`", xlab = "Antigen")

Var1.test <- wilcox.test(`A5024 hIgGM` ~ Antigen)
Var1.test

#`LS3 AQP4 11-30` 
group_by(var2, Antigen) %>%
  summarise(count = n(),
            mean = mean(`LS3 AQP4 11-30`, na.rm = TRUE),
            median = median(`LS3 AQP4 11-30`, na.rm = TRUE),
            sd = sd(`LS3 AQP4 11-30`, na.rm = TRUE))

ggboxplot(var2, x = "Antigen", y = "`LS3 AQP4 11-30`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`LS3 AQP4 11-30`", xlab = "Antigen")

Var2.test <- wilcox.test(`LS3 AQP4 11-30` ~ Antigen)
Var2.test

#`LS4 AQP4 21-40` 
group_by(var3, Antigen) %>%
  summarise(count = n(),
            mean = mean(`LS4 AQP4 21-40`, na.rm = TRUE),
            median = median(`LS4 AQP4 21-40`, na.rm = TRUE),
            sd = sd(`LS4 AQP4 21-40`, na.rm = TRUE))

ggboxplot(var3, x = "Antigen", y = "`LS4 AQP4 21-40`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`LS4 AQP4 21-40`", xlab = "Antigen")

Var3.test <- wilcox.test(`LS4 AQP4 21-40` ~ Antigen)
Var3.test

#`14 KIR4.1 131-150`
group_by(var4, Antigen) %>%
  summarise(count = n(),
            mean = mean(`14 KIR4.1 131-150`, na.rm = TRUE),
            median = median(`14 KIR4.1 131-150`, na.rm = TRUE),
            sd = sd(`14 KIR4.1 131-150`, na.rm = TRUE))

ggboxplot(var4, x = "Antigen", y = "`14 KIR4.1 131-150`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`14 KIR4.1 131-150`", xlab = "Antigen")

Var4.test <- wilcox.test(`14 KIR4.1 131-150` ~ Antigen)
Var4.test

#`A4874 MBP 31-50` 
group_by(var5, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4874 MBP 31-50`, na.rm = TRUE),
            median = median(`A4874 MBP 31-50`, na.rm = TRUE),
            sd = sd(`A4874 MBP 31-50`, na.rm = TRUE))

ggboxplot(var5, x = "Antigen", y = "`A4874 MBP 31-50`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4874 MBP 31-50`", xlab = "Antigen")

Var5.test <- wilcox.test(`A4874 MBP 31-50` ~ Antigen)
Var5.test

#`A4754 hAQP4 - Mayo in 0.3HAPS` 
group_by(var6, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4754 hAQP4 - Mayo in 0.3HAPS`, na.rm = TRUE),
            median = median(`A4754 hAQP4 - Mayo in 0.3HAPS`, na.rm = TRUE),
            sd = sd(`A4754 hAQP4 - Mayo in 0.3HAPS`, na.rm = TRUE))

ggboxplot(var6, x = "Antigen", y = "`A4754 hAQP4 - Mayo in 0.3HAPS`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4754 hAQP4 - Mayo in 0.3HAPS`", xlab = "Antigen")

Var6.test <- wilcox.test(`A4754 hAQP4 - Mayo in 0.3HAPS` ~ Antigen)
Var6.test

#`LS20 AQP4 181-200` 
group_by(var7, Antigen) %>%
  summarise(count = n(),
            mean = mean(`LS20 AQP4 181-200`, na.rm = TRUE),
            median = median(`LS20 AQP4 181-200`, na.rm = TRUE),
            sd = sd(`LS20 AQP4 181-200`, na.rm = TRUE))

ggboxplot(var7, x = "Antigen", y = "`LS20 AQP4 181-200`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`LS20 AQP4 181-200`", xlab = "Antigen")

Var7.test <- wilcox.test(`LS20 AQP4 181-200` ~ Antigen)
Var7.test

#Best model GBM
`ANO2 136-151` <- response$`ANO2 136-151`              	
`LS19 AQP4 171-190` <- response$`LS19 AQP4 171-190`             
`A4875 MBP 41-60` <- response$`A4875 MBP 41-60`         	
`A4752 hAQP4 - Mayo` <- response$`A4752 hAQP4 - Mayo`             
`A4895 MOBP 1-20` <- response$`A4895 MOBP 1-20`         	
`A4953 abCrys 71-91` <- response$`A4953 abCrys 71-91`              
`A4955 abCrys 91-110` <- response$`A4955 abCrys 91-110`    	
`P_AQP4  271-290` <- response$`P_AQP4  271-290`                 
`A4937 PLP 251-270` <- response$`A4937 PLP 251-270`	
`5 KIR4.1 41-60` <- response$`5 KIR4.1 41-60`
`A4905 MOBP 101-120` <- response$`A4905 MOBP 101-120`	
`LS15 AQP4 131-150` <- response$`LS15 AQP4 131-150`
`P_AQP4 141-160_B` <- response$`P_AQP4 141-160_B`


var8 <- as.data.frame(cbind(Antigen, `ANO2 136-151`))
var9 <- as.data.frame(cbind(Antigen, `LS19 AQP4 171-190`))
var10 <- as.data.frame(cbind(Antigen, `A4875 MBP 41-60`))
var11 <- as.data.frame(cbind(Antigen, `A4752 hAQP4 - Mayo`))
var12 <- as.data.frame(cbind(Antigen, `A4895 MOBP 1-20`))
var13 <- as.data.frame(cbind(Antigen, `A4953 abCrys 71-91`))
var14 <- as.data.frame(cbind(Antigen, `A4955 abCrys 91-110`))
var15 <- as.data.frame(cbind(Antigen, `P_AQP4  271-290`))
var16 <- as.data.frame(cbind(Antigen, `A4937 PLP 251-270`))
var17 <- as.data.frame(cbind(Antigen, `5 KIR4.1 41-60`))
var18 <- as.data.frame(cbind(Antigen, `A4905 MOBP 101-120`))
var19 <- as.data.frame(cbind(Antigen, `LS15 AQP4 131-150`))
var20 <- as.data.frame(cbind(Antigen, `P_AQP4 141-160_B`))

#`ANO2 136-151` 
group_by(var8, Antigen) %>%
  summarise(count = n(),
            mean = mean(`ANO2 136-151`, na.rm = TRUE),
            median = median(`ANO2 136-151`, na.rm = TRUE),
            sd = sd(`ANO2 136-151`, na.rm = TRUE))

ggboxplot(var8, x = "Antigen", y = "`ANO2 136-151`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`ANO2 136-151`", xlab = "Antigen")

Var8.test <- wilcox.test(`ANO2 136-151` ~ Antigen)
Var8.test

#`LS19 AQP4 171-190` 
group_by(var9, Antigen) %>%
  summarise(count = n(),
            mean = mean(`LS19 AQP4 171-190`, na.rm = TRUE),
            median = median(`LS19 AQP4 171-190`, na.rm = TRUE),
            sd = sd(`LS19 AQP4 171-190`, na.rm = TRUE))

ggboxplot(var9, x = "Antigen", y = "`LS19 AQP4 171-190`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`LS19 AQP4 171-190`", xlab = "Antigen")

Var9.test <- wilcox.test(`LS19 AQP4 171-190` ~ Antigen)
Var9.test

#`A4875 MBP 41-60` 
group_by(var10, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4875 MBP 41-60`, na.rm = TRUE),
            median = median(`A4875 MBP 41-60`, na.rm = TRUE),
            sd = sd(`A4875 MBP 41-60`, na.rm = TRUE))

ggboxplot(var10, x = "Antigen", y = "`A4875 MBP 41-60`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4875 MBP 41-60`", xlab = "Antigen")

Var10.test <- wilcox.test(`A4875 MBP 41-60` ~ Antigen)
Var10.test

#`A4752 hAQP4 - Mayo` 
group_by(var11, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4752 hAQP4 - Mayo`, na.rm = TRUE),
            median = median(`A4752 hAQP4 - Mayo`, na.rm = TRUE),
            sd = sd(`A4752 hAQP4 - Mayo`, na.rm = TRUE))

ggboxplot(var11, x = "Antigen", y = "`A4752 hAQP4 - Mayo`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4752 hAQP4 - Mayo`", xlab = "Antigen")

Var11.test <- wilcox.test(`A4752 hAQP4 - Mayo` ~ Antigen)
Var11.test

#`A4895 MOBP 1-20` 
group_by(var12, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4895 MOBP 1-20`, na.rm = TRUE),
            median = median(`A4895 MOBP 1-20`, na.rm = TRUE),
            sd = sd(`A4895 MOBP 1-20`, na.rm = TRUE))

ggboxplot(var12, x = "Antigen", y = "`A4895 MOBP 1-20`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4895 MOBP 1-20`", xlab = "Antigen")

Var12.test <- wilcox.test(`A4895 MOBP 1-20` ~ Antigen)
Var12.test

#`A4953 abCrys 71-91` 
group_by(var13, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4953 abCrys 71-91`, na.rm = TRUE),
            median = median(`A4953 abCrys 71-91`, na.rm = TRUE),
            sd = sd(`A4953 abCrys 71-91`, na.rm = TRUE))

ggboxplot(var13, x = "Antigen", y = "`A4953 abCrys 71-91`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4953 abCrys 71-91`", xlab = "Antigen")

Var13.test <- wilcox.test(`A4953 abCrys 71-91` ~ Antigen)
Var13.test

#`A4955 abCrys 91-110` 
group_by(var14, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4955 abCrys 91-110`, na.rm = TRUE),
            median = median(`A4955 abCrys 91-110`, na.rm = TRUE),
            sd = sd(`A4955 abCrys 91-110`, na.rm = TRUE))

ggboxplot(var14, x = "Antigen", y = "`A4955 abCrys 91-110`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4955 abCrys 91-110`", xlab = "Antigen")

Var14.test <- wilcox.test(`A4955 abCrys 91-110` ~ Antigen)
Var14.test

#`P_AQP4  271-290` 
group_by(var15, Antigen) %>%
  summarise(count = n(),
            mean = mean(`P_AQP4  271-290`, na.rm = TRUE),
            median = median(`P_AQP4  271-290`, na.rm = TRUE),
            sd = sd(`P_AQP4  271-290`, na.rm = TRUE))

ggboxplot(var15, x = "Antigen", y = "`P_AQP4  271-290`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`P_AQP4  271-290`", xlab = "Antigen")

Var15.test <- wilcox.test(`P_AQP4  271-290` ~ Antigen)
Var15.test

#`A4937 PLP 251-270` 
group_by(var16, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4937 PLP 251-270`, na.rm = TRUE),
            median = median(`A4937 PLP 251-270`, na.rm = TRUE),
            sd = sd(`A4937 PLP 251-270`, na.rm = TRUE))

ggboxplot(var16, x = "Antigen", y = "`A4937 PLP 251-270`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4937 PLP 251-270`", xlab = "Antigen")

Var16.test <- wilcox.test(`A4937 PLP 251-270` ~ Antigen)
Var16.test

#`5 KIR4.1 41-60` 
group_by(var17, Antigen) %>%
  summarise(count = n(),
            mean = mean(`5 KIR4.1 41-60`, na.rm = TRUE),
            median = median(`5 KIR4.1 41-60`, na.rm = TRUE),
            sd = sd(`5 KIR4.1 41-60`, na.rm = TRUE))

ggboxplot(var17, x = "Antigen", y = "`5 KIR4.1 41-60`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`5 KIR4.1 41-60`", xlab = "Antigen")

Var17.test <- wilcox.test(`5 KIR4.1 41-60` ~ Antigen)
Var17.test

#`A4905 MOBP 101-120` 
group_by(var18, Antigen) %>%
  summarise(count = n(),
            mean = mean(`A4905 MOBP 101-120`, na.rm = TRUE),
            median = median(`A4905 MOBP 101-120`, na.rm = TRUE),
            sd = sd(`A4905 MOBP 101-120`, na.rm = TRUE))

ggboxplot(var18, x = "Antigen", y = "`A4905 MOBP 101-120`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`A4905 MOBP 101-120`", xlab = "Antigen")

Var18.test <- wilcox.test(`A4905 MOBP 101-120` ~ Antigen)
Var18.test

#`LS15 AQP4 131-150` 
group_by(var19, Antigen) %>%
  summarise(count = n(),
            mean = mean(`LS15 AQP4 131-150`, na.rm = TRUE),
            median = median(`LS15 AQP4 131-150`, na.rm = TRUE),
            sd = sd(`LS15 AQP4 131-150`, na.rm = TRUE))

ggboxplot(var19, x = "Antigen", y = "`LS15 AQP4 131-150`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`LS15 AQP4 131-150`", xlab = "Antigen")

Var19.test <- wilcox.test(`LS15 AQP4 131-150` ~ Antigen)
Var19.test

#`P_AQP4 141-160_B`
group_by(var20, Antigen) %>%
  summarise(count = n(),
            mean = mean(`P_AQP4 141-160_B`, na.rm = TRUE),
            median = median(`P_AQP4 141-160_B`, na.rm = TRUE),
            sd = sd(`P_AQP4 141-160_B`, na.rm = TRUE))

ggboxplot(var20, x = "Antigen", y = "`P_AQP4 141-160_B`", 
          color = "Antigen", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "`P_AQP4 141-160_B`", xlab = "Antigen")

Var20.test <- wilcox.test(`P_AQP4 141-160_B` ~ Antigen)
Var20.test
