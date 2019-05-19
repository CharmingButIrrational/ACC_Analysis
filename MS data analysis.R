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
#Elastic net 
train.elastic <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(233)
elastic <- train(Samples ~ ., 
                      data = df_alt, 
                      method = "glmnet",
                      trControl = train.elastic,
                      preProc = c("center", "scale"),
                      tuneLength = 10)
importance.elastic <- varImp(elastic, scale = FALSE)
print(importance.elastic)
plot(importance.elastic, top = 20)

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

###Dataset too small
#GBM 
#train.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 6, 9, 10),
#                        n.trees = (0:50)*50, 
#                        shrinkage = seq(.0005, .05,.0005),
#                        n.minobsinnode = 10)
#set.seed(86)
#gbm <- train(Samples ~ ., data = df_alt,
#                          method = "gbm",
#                          trControl = train.gbm,
#                          preProcess = c("center", "scale"),
#                          tuneLength = 10,
#                          tuneGrid = gbmGrid,
#                          na.action = na.pass)
#importance.gbm <- varImp(gbm, scale = FALSE)
#print(importance.gbm)
#plot(importance.gbm)

################################################################
#Compare models
results <- resamples(list(ELA=elastic, SVM=svm, LVQ=lvq))
summary(results)
bwplot(results)

################################################


#Common predictors for  progress..yes.no.
LVQ.predictors <- importance.lvq$importance
SVM.predictors <- importance.svm$importance
elastic.predictors <- importance.elastic$importance
all(LVQ.predictors$X0 == LVQ.predictors$X1)
all(SVM.predictors$X0 == SVM.predictors$X1)
LVQ.predictors$X1 <- NULL
SVM.predictors$X1 <- NULL
colnames(LVQ.predictors)[1] <- "Overall"
colnames(SVM.predictors)[1] <- "Overall"
LVQ.predictors <- setDT(LVQ.predictors, keep.rownames = TRUE)[]
SVM.predictors <- setDT(SVM.predictors, keep.rownames = TRUE)[]
elastic.predictors <- setDT(elastic.predictors, keep.rownames = TRUE)[]
LVQ.predictors <- LVQ.predictors[order(-LVQ.predictors$Overall),]
SVM.predictors <- SVM.predictors[order(-SVM.predictors$Overall),]
elastic.predictors <- elastic.predictors[order(-elastic.predictors$Overall),]
LVQ.predictors.Sig <- head(LVQ.predictors, 20)
SVM.predictors.Sig <- head(SVM.predictors, 20)
#Adjust elastic selection based on printed variables
elastic.predictors.Sig <- head(elastic.predictors, 4)
LVQ.predictors.Sig$Overall <- NULL
SVM.predictors.Sig$Overall <- NULL
elastic.predictors.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.predictors.Sig <- unlist(LVQ.predictors.Sig)
SVM.predictors.Sig <- unlist(SVM.predictors.Sig)
elastic.predictors.Sig <- unlist(elastic.predictors.Sig)

library(VennDiagram)
venn.data.predictors <- list(LVQ.predictors.Sig, SVM.predictors.Sig, elastic.predictors.Sig)
grid.newpage()
venn.plot.predictors <- venn.diagram(x = list(LVQ.predictors.Sig=LVQ.predictors.Sig, SVM.predictors.Sig=SVM.predictors.Sig, elastic.predictors.Sig=elastic.predictors.Sig),
                                       filename=NULL, 
                                       fill = c("red", "blue", "green"),
                                       alpha = 0.50,
                                       col = "transparent")
grid.draw(venn.plot.predictors)
venn.intersect.predictors <- calculate.overlap(venn.data.predictors)
print(venn.intersect.predictors$a5)

#############################################################

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
#Analysis of selected variables
library(dplyr)
library(ggpubr)

Samples <- as.numeric(df_alt$Samples)


#Features selected by all methods
A4976.OSP.131.150 <- df_alt$A4976.OSP.131.150
A4858.MOG.N.term <- df_alt$A4858.MOG.N.term
A4910.MOBP.161.180 <- df_alt$A4910.MOBP.161.180

var1 <- as.data.frame(cbind(Samples, A4976.OSP.131.150))
var2 <- as.data.frame(cbind(Samples, A4858.MOG.N.term))
var3 <- as.data.frame(cbind(Samples, A4910.MOBP.161.180))

#Variable 1 analysis
group_by(var1, Samples) %>%
        summarise(count = n(),
            mean = mean(A4976.OSP.131.150, na.rm = TRUE),
            sd = sd(A4976.OSP.131.150, na.rm = TRUE))
ggboxplot(var1, x = "Samples", y = "A4976.OSP.131.150", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4976.OSP.131.150", xlab = "Samples")

Var1.test <- wilcox.test(Samples, A4976.OSP.131.150, alternative = "two.sided")
Var1.test

#Variable 2 analysis
group_by(var2, Samples) %>%
  summarise(count = n(),
            mean = mean(A4858.MOG.N.term, na.rm = TRUE),
            sd = sd(A4858.MOG.N.term, na.rm = TRUE))
ggboxplot(var2, x = "Samples", y = "A4858.MOG.N.term", 
          color = "Samples", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "A4858.MOG.N.term", xlab = "Samples")

Var2.test <- wilcox.test(Samples, A4858.MOG.N.term, alternative = "two.sided")
Var2.test

#Variable 3 analysis
group_by(var3, Samples) %>%
  summarise(count = n(),
            mean = mean(A4910.MOBP.161.180, na.rm = TRUE),
            sd = sd(A4910.MOBP.161.180, na.rm = TRUE))
ggboxplot(var3, x = "Samples", y = "A4910.MOBP.161.180", 
          color = "Samples", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("1", "2"),
          ylab = "A4910.MOBP.161.180", xlab = "Samples")

Var3.test <- wilcox.test(Samples, A4910.MOBP.161.180, alternative = "two.sided")
Var3.test

#Variables selected by best preforming method (LVQ)

A4943.CNPase.343.362.RM <- df_alt$A4943.CNPase.343.362.RM
A4963.OSP.1.20 <- df_alt$A4963.OSP.1.20
A4874.MBP.31.50 <- df_alt$A4874.MBP.31.50            
A4945.CNPase.369.388.RM <- df_alt$A4945.CNPase.369.388.RM
A4962.abCrys.161.176  <- df_alt$A4962.abCrys.161.176
A4913.PLP.11.30 <- df_alt$A4913.PLP.11.30             
A4977.OSP.141.160 <- df_alt$A4977.OSP.141.160         
A4545.MBP <- df_alt$A4545.MBP                     
A4607.Tubulin <- df_alt$A4607.Tubulin               
A4930.PLP.181.200 <- df_alt$A4930.PLP.181.200        
A4981.OSP.181.200 <- df_alt$A4981.OSP.181.200        
A4641.Tubulin <- df_alt$A4641.Tubulin               
A4908.MOBP.141.160 <- df_alt$A4908.MOBP.141.160  
A4561.NaV.1.6L.113 <- df_alt$A4561.NaV.1.6L.113        
LS8.AQP4.61.80 <- df_alt$LS8.AQP4.61.80               
A4887.MBP.150.171 <- df_alt$A4887.MBP.150.171         
A4947.abCrys.11.30 <- df_alt$A4947.abCrys.11.30   

var4 <- as.data.frame(cbind(Samples, A4943.CNPase.343.362.RM))
var5 <- as.data.frame(cbind(Samples, A4963.OSP.1.20))
var6 <- as.data.frame(cbind(Samples, A4874.MBP.31.50))
var7 <- as.data.frame(cbind(Samples, A4945.CNPase.369.388.RM))
var8 <- as.data.frame(cbind(Samples, A4962.abCrys.161.176))
var9 <- as.data.frame(cbind(Samples, A4913.PLP.11.30))
var10 <- as.data.frame(cbind(Samples, A4977.OSP.141.160))
var11 <- as.data.frame(cbind(Samples, A4545.MBP))
var12 <- as.data.frame(cbind(Samples, A4607.Tubulin))
var13 <- as.data.frame(cbind(Samples, A4930.PLP.181.200))
var14 <- as.data.frame(cbind(Samples, A4981.OSP.181.200))
var15 <- as.data.frame(cbind(Samples, A4641.Tubulin))
var16 <- as.data.frame(cbind(Samples, A4908.MOBP.141.160))
var17 <- as.data.frame(cbind(Samples, A4561.NaV.1.6L.113))
var18 <- as.data.frame(cbind(Samples, LS8.AQP4.61.80))
var19 <- as.data.frame(cbind(Samples, A4887.MBP.150.171))
var20 <- as.data.frame(cbind(Samples, A4947.abCrys.11.30))

#Variable 4 analysis
group_by(var4, Samples) %>%
  summarise(count = n(),
            mean = mean(A4943.CNPase.343.362.RM, na.rm = TRUE),
            sd = sd(A4943.CNPase.343.362.RM, na.rm = TRUE))
ggboxplot(var4, x = "Samples", y = "A4943.CNPase.343.362.RM", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4943.CNPase.343.362.RM", xlab = "Samples")

Var4.test <- wilcox.test(Samples, A4943.CNPase.343.362.RM, alternative = "two.sided")
Var4.test

#Variable 5 analysis
group_by(var5, Samples) %>%
  summarise(count = n(),
            mean = mean(A4963.OSP.1.20, na.rm = TRUE),
            sd = sd(A4963.OSP.1.20, na.rm = TRUE))
ggboxplot(var5, x = "Samples", y = "A4963.OSP.1.20", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4963.OSP.1.20", xlab = "Samples")

Var5.test <- wilcox.test(Samples, A4963.OSP.1.20, alternative = "two.sided")
Var5.test

#Variable 6 analysis
group_by(var6, Samples) %>%
  summarise(count = n(),
            mean = mean(A4874.MBP.31.50, na.rm = TRUE),
            sd = sd(A4874.MBP.31.50, na.rm = TRUE))
ggboxplot(var6, x = "Samples", y = "A4874.MBP.31.50", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4874.MBP.31.50", xlab = "Samples")

Var6.test <- wilcox.test(Samples, A4874.MBP.31.50, alternative = "two.sided")
Var6.test

#Variable 7 analysis
group_by(var7, Samples) %>%
  summarise(count = n(),
            mean = mean(A4945.CNPase.369.388.RM, na.rm = TRUE),
            sd = sd(A4945.CNPase.369.388.RM, na.rm = TRUE))
ggboxplot(var7, x = "Samples", y = "A4945.CNPase.369.388.RM", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4945.CNPase.369.388.RM", xlab = "Samples")

Var7.test <- wilcox.test(Samples, A4945.CNPase.369.388.RM, alternative = "two.sided")
Var7.test

#Variable 8 analysis
group_by(var8, Samples) %>%
  summarise(count = n(),
            mean = mean(A4962.abCrys.161.176, na.rm = TRUE),
            sd = sd(A4962.abCrys.161.176, na.rm = TRUE))
ggboxplot(var8, x = "Samples", y = "A4962.abCrys.161.176", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4962.abCrys.161.176", xlab = "Samples")

Var8.test <- wilcox.test(Samples, A4962.abCrys.161.176, alternative = "two.sided")
Var8.test

#Variable 9 analysis
group_by(var9, Samples) %>%
  summarise(count = n(),
            mean = mean(A4913.PLP.11.30, na.rm = TRUE),
            sd = sd(A4913.PLP.11.30, na.rm = TRUE))
ggboxplot(var9, x = "Samples", y = "A4913.PLP.11.30", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4913.PLP.11.30", xlab = "Samples")

Var9.test <- wilcox.test(Samples, A4913.PLP.11.30, alternative = "two.sided")
Var9.test

#Variable 10 analysis
group_by(var10, Samples) %>%
  summarise(count = n(),
            mean = mean(A4977.OSP.141.160, na.rm = TRUE),
            sd = sd(A4977.OSP.141.160, na.rm = TRUE))
ggboxplot(var10, x = "Samples", y = "A4977.OSP.141.160", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4977.OSP.141.160", xlab = "Samples")

Var10.test <- wilcox.test(Samples, A4977.OSP.141.160, alternative = "two.sided")
Var10.test

#Variable 11 analysis
group_by(var11, Samples) %>%
  summarise(count = n(),
            mean = mean(A4545.MBP, na.rm = TRUE),
            sd = sd(A4545.MBP, na.rm = TRUE))
ggboxplot(var11, x = "Samples", y = "A4545.MBP", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4545.MBP", xlab = "Samples")

Var11.test <- wilcox.test(Samples, A4545.MBP, alternative = "two.sided")
Var11.test

#Variable 12 analysis
group_by(var12, Samples) %>%
  summarise(count = n(),
            mean = mean(A4607.Tubulin, na.rm = TRUE),
            sd = sd(A4607.Tubulin, na.rm = TRUE))
ggboxplot(var12, x = "Samples", y = "A4607.Tubulin", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4607.Tubulin", xlab = "Samples")

Var12.test <- wilcox.test(Samples, A4607.Tubulin, alternative = "two.sided")
Var12.test

#Variable 13 analysis
group_by(var13, Samples) %>%
  summarise(count = n(),
            mean = mean(A4930.PLP.181.200, na.rm = TRUE),
            sd = sd(A4930.PLP.181.200, na.rm = TRUE))
ggboxplot(var13, x = "Samples", y = "A4930.PLP.181.200", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4930.PLP.181.200", xlab = "Samples")

Var13.test <- wilcox.test(Samples, A4930.PLP.181.200, alternative = "two.sided")
Var13.test

#Variable 14 analysis
group_by(var14, Samples) %>%
  summarise(count = n(),
            mean = mean(A4981.OSP.181.200, na.rm = TRUE),
            sd = sd(A4981.OSP.181.200, na.rm = TRUE))
ggboxplot(var14, x = "Samples", y = "A4981.OSP.181.200", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4981.OSP.181.200", xlab = "Samples")

Var14.test <- wilcox.test(Samples, A4981.OSP.181.200, alternative = "two.sided")
Var14.test

#Variable 15 analysis
group_by(var15, Samples) %>%
  summarise(count = n(),
            mean = mean(A4641.Tubulin, na.rm = TRUE),
            sd = sd(A4641.Tubulin, na.rm = TRUE))
ggboxplot(var15, x = "Samples", y = "A4641.Tubulin", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4641.Tubulin", xlab = "Samples")

Var15.test <- wilcox.test(Samples, A4641.Tubulin, alternative = "two.sided")
Var15.test

#Variable 16 analysis
group_by(var16, Samples) %>%
  summarise(count = n(),
            mean = mean(A4908.MOBP.141.160, na.rm = TRUE),
            sd = sd(A4908.MOBP.141.160, na.rm = TRUE))
ggboxplot(var16, x = "Samples", y = "A4908.MOBP.141.160", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4908.MOBP.141.160", xlab = "Samples")

Var16.test <- wilcox.test(Samples, A4908.MOBP.141.160, alternative = "two.sided")
Var16.test

#Variable 17 analysis
group_by(var17, Samples) %>%
  summarise(count = n(),
            mean = mean(A4561.NaV.1.6L.113, na.rm = TRUE),
            sd = sd(A4561.NaV.1.6L.113, na.rm = TRUE))
ggboxplot(var17, x = "Samples", y = "A4561.NaV.1.6L.113", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4561.NaV.1.6L.113", xlab = "Samples")

Var17.test <- wilcox.test(Samples, A4561.NaV.1.6L.113, alternative = "two.sided")
Var17.test

#Variable 18 analysis
group_by(var18, Samples) %>%
  summarise(count = n(),
            mean = mean(LS8.AQP4.61.80, na.rm = TRUE),
            sd = sd(LS8.AQP4.61.80, na.rm = TRUE))
ggboxplot(var18, x = "Samples", y = "LS8.AQP4.61.80", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "LS8.AQP4.61.80", xlab = "Samples")

Var18.test <- wilcox.test(Samples, LS8.AQP4.61.80, alternative = "two.sided")
Var18.test

#Variable 19 analysis
group_by(var19, Samples) %>%
  summarise(count = n(),
            mean = mean(A4887.MBP.150.171, na.rm = TRUE),
            sd = sd(A4887.MBP.150.171, na.rm = TRUE))
ggboxplot(var19, x = "Samples", y = "A4887.MBP.150.171", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4887.MBP.150.171", xlab = "Samples")

Var19.test <- wilcox.test(Samples, A4887.MBP.150.171, alternative = "two.sided")
Var19.test

#Variable 20 analysis
group_by(var20, Samples) %>%
  summarise(count = n(),
            mean = mean(A4947.abCrys.11.30, na.rm = TRUE),
            sd = sd(A4947.abCrys.11.30, na.rm = TRUE))
ggboxplot(var20, x = "Samples", y = "A4947.abCrys.11.30", 
          color = "Samples",
          order = c("1", "2"),
          ylab = "A4947.abCrys.11.30", xlab = "Samples")

Var20.test <- wilcox.test(Samples, A4947.abCrys.11.30, alternative = "two.sided")
Var20.test

