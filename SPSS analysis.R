#Pheochromocytoma analysis

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/SPSS.csv", sep = ";", header = T)

class <- as.factor(mydata$class)

alt.data <- mydata[,21:150]

alt.data <- cbind(class, alt.data)

library(gbm)
library(glmnet)
library(caret)


#Elastic net using caret
train.elastic <- trainControl(method = "repeatedcv", number = 10, repeats = 5, search = "random")
set.seed(233)
elastic.spss <- train(class ~ ., 
                      data = alt.data, 
                      method = "glmnet",
                      trControl = train.elastic,
                      preProc = c("center", "scale"),
                      tuneLength = 10)
importance.elastic <- varImp(elastic.spss, scale = FALSE)
print(importance.elastic)
plot(importance.elastic, top = 20)


#SVM
train.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5, search = "random")
set.seed(233)
svm.spss <- train(class ~ ., 
                      data = alt.data, 
                      method = "svmLinear",
                      trControl = train.svm,
                      preProc = c("center", "scale"),
                      tuneLength = 10,
                      allowParallel = T)
importance.svm <- varImp(svm.spss, scale = FALSE)
print(importance.svm)
plot(importance.svm, top = 20)


#GBM
train.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
caretGrid <- expand.grid(interaction.depth=c(1, 3, 5), n.trees = (0:50)*50,
                         shrinkage=c(0.01, 0.001),
                         n.minobsinnode=10)
set.seed(233)
gbm.spss <- train(class ~ ., 
                  data = alt.data, 
                  method = "gbm",
                  trControl = train.gbm,
                  preProc = c("center", "scale"),
                  tuneGrid = caretGrid,
                  tuneLength = 10)
importance.gbm <- varImp(gbm.spss, scale = FALSE)   #9 variables of importance
print(importance.gbm)
plot(importance.gbm, top = 20)

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
#Analysis of selected variables
library(dplyr)
library(ggpubr)

#Features selected by elastic net
H1 <- alt.data$H1
Ala <- alt.data$Ala

var1 <- as.data.frame(cbind(class, H1))
var2 <- as.data.frame(cbind(class, Ala))

class <- as.numeric(alt.data$class)

#Variable 1
group_by(var1, class) %>%
  summarise(count = n(),
            mean = mean(H1, na.rm = TRUE),
            sd = sd(H1, na.rm = TRUE))

ggboxplot(var1, x = "class", y = "H1", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "H1", xlab = "Class")

Var1.test <- wilcox.test(class, H1, alternative = "two.sided")
Var1.test

#Variable 2
group_by(var2, class) %>%
  summarise(count = n(),
            mean = mean(Ala, na.rm = TRUE),
            sd = sd(Ala, na.rm = TRUE))

ggboxplot(var2, x = "class", y = "Ala", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "Ala", xlab = "Class")

Var2.test <- wilcox.test(class, Ala, alternative = "two.sided")
Var2.test

##################################################################
#Variables selected from best model

Gly <- alt.data$Gly           
PC.aa.C34.1 <- alt.data$PC.aa.C34.1         
Lys <- alt.data$Lys               
Val <- alt.data$Val                
Gln <- alt.data$Gln               
PC.aa.C36.2 <- alt.data$PC.aa.C36.2  
PC.aa.C34.2 <- alt.data$PC.aa.C34.2    

var3 <- as.data.frame(cbind(class, Gly))
var4 <- as.data.frame(cbind(class, PC.aa.C34.1))
var5 <- as.data.frame(cbind(class, Lys))
var6 <- as.data.frame(cbind(class, Val))
var7 <- as.data.frame(cbind(class, Gln))
var8 <- as.data.frame(cbind(class, PC.aa.C36.2))
var9 <- as.data.frame(cbind(class, PC.aa.C34.2))

#Variable 3
group_by(var3, class) %>%
  summarise(count = n(),
            mean = mean(Gly, na.rm = TRUE),
            sd = sd(Gly, na.rm = TRUE))

ggboxplot(var1, x = "class", y = "Gly", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "Gly", xlab = "Class")

Var3.test <- wilcox.test(class, Gly, alternative = "two.sided")
Var3.test

#Variable 4
group_by(var4, class) %>%
  summarise(count = n(),
            mean = mean(PC.aa.C34.1, na.rm = TRUE),
            sd = sd(PC.aa.C34.1, na.rm = TRUE))

ggboxplot(var4, x = "class", y = "PC.aa.C34.1", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "PC.aa.C34.1", xlab = "Class")

Var4.test <- wilcox.test(class, PC.aa.C34.1, alternative = "two.sided")
Var4.test

#Variable 5
group_by(var5, class) %>%
  summarise(count = n(),
            mean = mean(Lys, na.rm = TRUE),
            sd = sd(Lys, na.rm = TRUE))

ggboxplot(var1, x = "class", y = "Lys", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "Lys", xlab = "Class")

Var5.test <- wilcox.test(class, Lys, alternative = "two.sided")
Var5.test

#Variable 6
group_by(var6, class) %>%
  summarise(count = n(),
            mean = mean(Val, na.rm = TRUE),
            sd = sd(Val, na.rm = TRUE))

ggboxplot(var6, x = "class", y = "Val", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "Val", xlab = "Class")

Var6.test <- wilcox.test(class, Val, alternative = "two.sided")
Var6.test

#Variable 7
group_by(var7, class) %>%
  summarise(count = n(),
            mean = mean(Gln, na.rm = TRUE),
            sd = sd(Gln, na.rm = TRUE))

ggboxplot(var7, x = "class", y = "Gln", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "Gln", xlab = "Class")

Var7.test <- wilcox.test(class, Gln, alternative = "two.sided")
Var7.test

#Variable 8
group_by(var8, class) %>%
  summarise(count = n(),
            mean = mean(PC.aa.C36.2, na.rm = TRUE),
            sd = sd(PC.aa.C36.2, na.rm = TRUE))

ggboxplot(var8, x = "class", y = "PC.aa.C36.2", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "PC.aa.C36.2", xlab = "Class")

Var8.test <- wilcox.test(class, PC.aa.C36.2, alternative = "two.sided")
Var8.test

#Variable 9
group_by(var9, class) %>%
  summarise(count = n(),
            mean = mean(PC.aa.C34.2, na.rm = TRUE),
            sd = sd(PC.aa.C34.2, na.rm = TRUE))

ggboxplot(var9, x = "class", y = "PC.aa.C34.2", 
          color = "class", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "PC.aa.C34.2", xlab = "Class")

Var9.test <- wilcox.test(class, PC.aa.C34.2, alternative = "two.sided")
Var9.test
