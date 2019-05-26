#Pheochromocytoma analysis (malignancy)

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/SPSS.csv", sep = ";", header = T)

malignancy <- as.factor(mydata$malignancy)

alt.data <- mydata[,21:150]

alt.data <- cbind(phenotype, alt.data)

library(gbm)
library(glmnet)
library(caret)


#Elastic net using caret
train.elastic <- trainControl(method = "repeatedcv", number = 10, repeats = 5, search = "random")
set.seed(233)
elastic.spss <- train(malignancy ~ ., 
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
svm.spss <- train(malignancy ~ ., 
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
caretGrid <- expand.grid(interaction.depth=c(1, 3, 5), n.trees = (0:500),
                         shrinkage=c(0.01, 0.001),
                         n.minobsinnode=10)
set.seed(233)
gbm.spss <- train(malignancy ~ ., 
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
#Venn diagram script doesn't work with multivariate 
###################################################################
#Analysis of selected variables
library(dplyr)
library(ggpubr)

malignancy <- as.numeric(mydata$malignancy)

H1 <- alt.data$H1                  
Gly <- alt.data$Gly                 
Ala <- alt.data$Ala                 
PC.aa.C34.2 <- alt.data$PC.aa.C34.2      
Val <- alt.data$Val                
Lys <- alt.data$Lys                 
PC.aa.C34.1 <- alt.data$PC.aa.C34.1       
PC.aa.C36.2 <- alt.data$PC.aa.C36.2        
Gln <- alt.data$Gln 

var1 <- as.data.frame(cbind(malignancy, H1))
var2 <- as.data.frame(cbind(malignancy, Gly))
var3 <- as.data.frame(cbind(malignancy, Ala))
var4 <- as.data.frame(cbind(malignancy, PC.aa.C34.2))
var5 <- as.data.frame(cbind(malignancy, Val))
var6 <- as.data.frame(cbind(malignancy, Lys))
var7 <- as.data.frame(cbind(malignancy, PC.aa.C34.1))
var8 <- as.data.frame(cbind(malignancy, PC.aa.C36.2))
var9 <- as.data.frame(cbind(malignancy, Gln))

#Variable 1
group_by(var1, malignancy) %>%
  summarise(count = n(),
            mean = mean(H1, na.rm = TRUE),
            median = median(H1, na.rm = TRUE),
            sd = sd(H1, na.rm = TRUE))

ggboxplot(var1, x = "malignancy", y = "H1", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "H1", xlab = "malignancy")

Var1.test <- pairwise.wilcox.test(var1$H1, var1$malignancy, p.adjust.method = "BH")
Var1.test

#Variable 2
group_by(var2, malignancy) %>%
  summarise(count = n(),
            mean = mean(Gly, na.rm = TRUE),
            median = median(Gly, na.rm = TRUE),
            sd = sd(Gly, na.rm = TRUE))

ggboxplot(var2, x = "malignancy", y = "Gly", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "Gly", xlab = "malignancy")

Var2.test <- pairwise.wilcox.test(var2$Gly, var2$malignancy, p.adjust.method = "BH")
Var2.test

#Variable 3
group_by(var3, malignancy) %>%
  summarise(count = n(),
            mean = mean(Ala, na.rm = TRUE),
            median = median(Ala, na.rm = TRUE),
            sd = sd(Ala, na.rm = TRUE))

ggboxplot(var3, x = "malignancy", y = "Ala", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "Ala", xlab = "malignancy")

Var3.test <- pairwise.wilcox.test(var3$Ala, var3$malignancy, p.adjust.method = "BH")
Var3.test

#Variable 4
group_by(var4, malignancy) %>%
  summarise(count = n(),
            mean = mean(PC.aa.C34.2, na.rm = TRUE),
            median = median(PC.aa.C34.2, na.rm = TRUE),
            sd = sd(PC.aa.C34.2, na.rm = TRUE))

ggboxplot(var4, x = "malignancy", y = "PC.aa.C34.2", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "PC.aa.C34.2", xlab = "malignancy")

Var4.test <- pairwise.wilcox.test(var4$PC.aa.C34.2, var4$malignancy, p.adjust.method = "BH")
Var4.test

#Variable 5
group_by(var5, malignancy) %>%
  summarise(count = n(),
            mean = mean(Val, na.rm = TRUE),
            median = median(Val, na.rm = TRUE),
            sd = sd(Val, na.rm = TRUE))

ggboxplot(var5, x = "malignancy", y = "Val", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "Val", xlab = "malignancy")

Var5.test <- pairwise.wilcox.test(var5$Val, var5$malignancy, p.adjust.method = "BH")
Var5.test

#Variable 6
group_by(var6, malignancy) %>%
  summarise(count = n(),
            mean = mean(Lys, na.rm = TRUE),
            median = median(Lys, na.rm = TRUE),
            sd = sd(Lys, na.rm = TRUE))

ggboxplot(var6, x = "malignancy", y = "Lys", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "Lys", xlab = "malignancy")

Var6.test <- pairwise.wilcox.test(var6$Lys, var6$malignancy, p.adjust.method = "BH")
Var6.test

#Variable 7
group_by(var7, malignancy) %>%
  summarise(count = n(),
            mean = mean(PC.aa.C34.1, na.rm = TRUE),
            median = median(PC.aa.C34.1, na.rm = TRUE),
            sd = sd(PC.aa.C34.1, na.rm = TRUE))

ggboxplot(var7, x = "malignancy", y = "PC.aa.C34.1", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "PC.aa.C34.1", xlab = "malignancy")

Var7.test <- pairwise.wilcox.test(var7$PC.aa.C34.1, var7$malignancy, p.adjust.method = "BH")
Var7.test

#Variable 8
group_by(var8, malignancy) %>%
  summarise(count = n(),
            mean = mean(PC.aa.C36.2, na.rm = TRUE),
            median = median(PC.aa.C36.2, na.rm = TRUE),
            sd = sd(PC.aa.C36.2, na.rm = TRUE))

ggboxplot(var8, x = "malignancy", y = "PC.aa.C36.2", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "PC.aa.C36.2", xlab = "malignancy")

Var8.test <- pairwise.wilcox.test(var8$PC.aa.C36.2, var8$malignancy, p.adjust.method = "BH")
Var8.test

#Variable 9
group_by(var9, malignancy) %>%
  summarise(count = n(),
            mean = mean(Gln, na.rm = TRUE),
            median = median(Gln, na.rm = TRUE),
            sd = sd(Gln, na.rm = TRUE))

ggboxplot(var9, x = "malignancy", y = "Gln", 
          color = "malignancy", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2", "3"),
          ylab = "Gln", xlab = "malignancy")

Var9.test <- pairwise.wilcox.test(var9$Gln, var9$malignancy, p.adjust.method = "BH")
Var9.test