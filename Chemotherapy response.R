#Load data
mydata <- read.csv("C:/Users/User/Desktop/Analysis Data/SepMutData_Oisin")

mydata$X <- NULL #column numbers appeared as a seperate column
alt.data.whole <- mydata[,c(34:180)]
alt.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(mydata$progress..yes.no.)
alt.data.whole <- cbind(progress..yes.no., alt.data.whole)

#Seperating Whole dataset into A_, B_, C_
alt.data.whole.A <- alt.data.whole[,c(1:72)]
alt.data.whole.B <- alt.data.whole[,c(73:105)]
alt.data.whole.B <- cbind(progress..yes.no., alt.data.whole.B)
alt.data.whole.C <- alt.data.whole[,c(105:146)]
alt.data.whole.C <- cbind(progress..yes.no., alt.data.whole.C)

#Split df by value in the Chemotherapy column
alt.data <- split(mydata, mydata$Chemotherapy)
alt.data0 <- alt.data$`0`
alt.data1 <- alt.data$`1`
#For treatment without chemotherapy
alt.data00 <- alt.data0[,c(34:180)]
alt.data00$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data00$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(alt.data00$progress..yes.no.)
alt.data00 <- cbind(progress..yes.no., alt.data00)
#For treatment with chemotherapy 
alt.data11 <- alt.data1[,c(34:180)]
alt.data11$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data11$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(alt.data1$progress..yes.no.)
alt.data11 <- cbind(progress..yes.no., alt.data11)

#Seperating Chemo dataset into A_, B_, C_
alt.data11.A <- alt.data11[,c(1:72)]
alt.data11.B <- alt.data11[,c(73:105)]
alt.data11.B <- cbind(progress..yes.no., alt.data11.B)
alt.data11.C <- alt.data11[,c(105:146)]
alt.data11.C <- cbind(progress..yes.no., alt.data11.C)

library(caret)
library(gbm)
library(MASS)

######################################################################################################################################
#Whole dataset
#LVQ A_
test.whole.lvq.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.lvq.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF,
                               data = alt.data.whole, method = "lvq",
                               trControl = test.whole.lvq.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.lvq.pro <- varImp(model.whole.lvq.pro, scale = FALSE)
print(importance.whole.lvq.pro)
plot(importance.whole.lvq.pro)

#SVM A_
test.whole.svm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.svm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                               data = alt.data.whole, method = "svmRadial",
                               trControl = test.whole.svm.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.svm.pro <- varImp(model.whole.svm.pro, scale = FALSE)
print(importance.whole.svm.pro)
plot(importance.whole.svm.pro)

#GBM A_
test.whole.gbm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.gbm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                               data = alt.data.whole, method = "gbm",
                               trControl = test.whole.gbm.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)

importance.whole.gbm.pro <- varImp(model.whole.gbm.pro, scale = FALSE)
print(importance.whole.gbm.pro)
plot(importance.whole.gbm.pro)

#Common predictors for A_ progress..yes.no.
LVQ.whole.Pro <- importance.whole.lvq.pro$importance
SVM.whole.Pro <- importance.whole.svm.pro$importance
GBM.whole.Pro <- importance.whole.gbm.pro$importance
all(LVQ.whole.Pro$X0 == LVQ.whole.Pro$X1)
all(SVM.whole.Pro$X0 == SVM.whole.Pro$X1)
LVQ.whole.Pro$X1 <- NULL
SVM.whole.Pro$X1 <- NULL
colnames(LVQ.whole.Pro)[1] <- "Overall"
colnames(SVM.whole.Pro)[1] <- "Overall"
LVQ.whole.Pro <- setDT(LVQ.whole.Pro, keep.rownames = TRUE)[]
SVM.whole.Pro <- setDT(SVM.whole.Pro, keep.rownames = TRUE)[]
GBM.whole.Pro <- setDT(GBM.whole.Pro, keep.rownames = TRUE)[]
LVQ.whole.Pro <- LVQ.whole.Pro[order(-LVQ.whole.Pro$Overall),]
SVM.whole.Pro <- SVM.whole.Pro[order(-SVM.whole.Pro$Overall),]
GBM.whole.Pro <- GBM.whole.Pro[order(-GBM.whole.Pro$Overall),]
LVQ.whole.Pro.Sig <- head(LVQ.whole.Pro, 20)
SVM.whole.Pro.Sig <- head(SVM.whole.Pro, 20)
#Adjust GBM selection based on printed variables
GBM.whole.Pro.Sig <- head(GBM.whole.Pro, 20)
library(data.table)
LVQ.whole.Pro.Sig$Overall <- NULL
SVM.whole.Pro.Sig$Overall <- NULL
GBM.whole.Pro.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.whole.Pro.Sig <- unlist(LVQ.whole.Pro.Sig)
SVM.whole.Pro.Sig <- unlist(SVM.whole.Pro.Sig)
GBM.whole.Pro.Sig <- unlist(GBM.whole.Pro.Sig)
library(VennDiagram)
venn.data.whole <- list(LVQ.whole.Pro.Sig, SVM.whole.Pro.Sig, GBM.whole.Pro.Sig)
grid.newpage()
venn.plot.whole <- venn.diagram(x = list(LVQ.whole.Pro.Sig=LVQ.whole.Pro.Sig, SVM.whole.Pro.Sig=SVM.whole.Pro.Sig, GBM.whole.Pro.Sig=GBM.whole.Pro.Sig),
                                  filename=NULL, 
                                  fill = c("red", "blue", "green"),
                                  alpha = 0.50,
                                  col = "transparent")
grid.draw(venn.plot.whole)
venn.intersect.whole <- calculate.overlap(venn.data.whole)
print(venn.intersect.whole$a5)


#Seperation for A_, B_ and C_ using progress..yes.no. Whole dataset
#LVQ A_
test.whole.lvq.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.lvq.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF,
                               data = alt.data.whole.A, method = "lvq",
                               trControl = test.whole.lvq.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.lvq.A.pro <- varImp(model.whole.lvq.A.pro, scale = FALSE)
print(importance.whole.lvq.A.pro)
plot(importance.whole.lvq.A.pro)

#SVM A_
test.whole.svm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.svm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                               data = alt.data.whole.A, method = "svmRadial",
                               trControl = test.whole.svm.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.svm.A.pro <- varImp(model.whole.svm.A.pro, scale = FALSE)
print(importance.whole.svm.A.pro)
plot(importance.whole.svm.A.pro)

#GBM A_
test.whole.gbm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.gbm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                               data = alt.data.whole.A, method = "gbm",
                               trControl = test.whole.gbm.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)

importance.whole.gbm.A.pro <- varImp(model.whole.gbm.A.pro, scale = FALSE)
print(importance.whole.gbm.A.pro)
plot(importance.whole.gbm.A.pro)

#Common predictors for A_ progress..yes.no.
LVQ.whole.Pro.A <- importance.whole.lvq.A.pro$importance
SVM.whole.Pro.A <- importance.whole.svm.A.pro$importance
GBM.whole.Pro.A <- importance.whole.gbm.A.pro$importance
all(LVQ.whole.Pro.A$X0 == LVQ.whole.Pro.A$X1)
all(SVM.whole.Pro.A$X0 == SVM.whole.Pro.A$X1)
LVQ.whole.Pro.A$X1 <- NULL
SVM.whole.Pro.A$X1 <- NULL
colnames(LVQ.whole.Pro.A)[1] <- "Overall"
colnames(SVM.whole.Pro.A)[1] <- "Overall"
LVQ.whole.Pro.A <- setDT(LVQ.whole.Pro.A, keep.rownames = TRUE)[]
SVM.whole.Pro.A <- setDT(SVM.whole.Pro.A, keep.rownames = TRUE)[]
GBM.whole.Pro.A <- setDT(GBM.whole.Pro.A, keep.rownames = TRUE)[]
LVQ.whole.Pro.A <- LVQ.whole.Pro.A[order(-LVQ.whole.Pro.A$Overall),]
SVM.whole.Pro.A <- SVM.whole.Pro.A[order(-SVM.whole.Pro.A$Overall),]
GBM.whole.Pro.A <- GBM.whole.Pro.A[order(-GBM.whole.Pro.A$Overall),]
LVQ.whole.Pro.Sig.A <- head(LVQ.whole.Pro.A, 20)
SVM.whole.Pro.Sig.A <- head(SVM.whole.Pro.A, 20)
#Adjust GBM selection based on printed variables
GBM.whole.Pro.Sig.A <- head(GBM.whole.Pro.A, 20)
library(data.table)
LVQ.whole.Pro.Sig.A$Overall <- NULL
SVM.whole.Pro.Sig.A$Overall <- NULL
GBM.whole.Pro.Sig.A$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.whole.Pro.Sig.A <- unlist(LVQ.whole.Pro.Sig.A)
SVM.whole.Pro.Sig.A <- unlist(SVM.whole.Pro.Sig.A)
GBM.whole.Pro.Sig.A <- unlist(GBM.whole.Pro.Sig.A)
library(VennDiagram)
venn.data.whole.A <- list(LVQ.whole.Pro.Sig.A, SVM.whole.Pro.Sig.A, GBM.whole.Pro.Sig.A)
grid.newpage()
venn.plot.whole.A <- venn.diagram(x = list(LVQ.whole.Pro.Sig.A=LVQ.whole.Pro.Sig.A, SVM.whole.Pro.Sig.A=SVM.whole.Pro.Sig.A, GBM.whole.Pro.Sig.A=GBM.whole.Pro.Sig.A),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.whole.A)
venn.intersect.whole.A <- calculate.overlap(venn.data.whole.A)
print(venn.intersect.whole.A$a5)

#LVQ B_
test.whole.lvq.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.lvq.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data.whole.B, method = "lvq",
                               trControl = test.whole.lvq.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.lvq.B.pro <- varImp(model.whole.lvq.B.pro, scale = FALSE)
print(importance.whole.lvq.B.pro)
plot(importance.whole.lvq.B.pro)

#SVM B_
test.whole.svm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.svm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data.whole.B, method = "svmRadial",
                               trControl = test.whole.svm.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.svm.B.pro <- varImp(model.whole.svm.B.pro, scale = FALSE)
print(importance.whole.svm.B.pro)
plot(importance.whole.svm.B.pro)

#GBM B_
test.whole.gbm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.gbm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data.whole.B, method = "gbm",
                               trControl = test.whole.gbm.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.gbm.B.pro <- varImp(model.whole.gbm.B.pro, scale = FALSE)
print(importance.whole.gbm.B.pro)
plot(importance.whole.gbm.B.pro)

#Common predictors for A_ progress..yes.no.
LVQ.whole.Pro.B <- importance.whole.lvq.B.pro$importance
SVM.whole.Pro.B <- importance.whole.svm.B.pro$importance
GBM.whole.Pro.B <- importance.whole.gbm.B.pro$importance
all(LVQ.whole.Pro.B$X0 == LVQ.whole.Pro.B$X1)
all(SVM.whole.Pro.B$X0 == SVM.whole.Pro.B$X1)
LVQ.whole.Pro.B$X1 <- NULL
SVM.whole.Pro.B$X1 <- NULL
colnames(LVQ.whole.Pro.B)[1] <- "Overall"
colnames(SVM.whole.Pro.B)[1] <- "Overall"
LVQ.whole.Pro.Sig.B <- setDT(LVQ.whole.Pro.B, keep.rownames = TRUE)[]
SVM.whole.Pro.Sig.B <- setDT(SVM.whole.Pro.B, keep.rownames = TRUE)[]
GBM.whole.Pro.Sig.B <- setDT(GBM.whole.Pro.B, keep.rownames = TRUE)[]
LVQ.whole.Pro.Sig.B <- LVQ.whole.Pro.Sig.B[order(-LVQ.whole.Pro.B$Overall),]
SVM.whole.Pro.Sig.B <- SVM.whole.Pro.Sig.B[order(-SVM.whole.Pro.B$Overall),]
GBM.whole.Pro.Sig.B <- GBM.whole.Pro.Sig.B[order(-GBM.whole.Pro.B$Overall),]
LVQ.whole.Pro.Sig.B <- head(LVQ.whole.Pro.B, 20)
SVM.whole.Pro.Sig.B <- head(SVM.whole.Pro.B, 20)
#Adjust GBM selection based on printed variables
GBM.whole.Pro.Sig.B <- head(GBM.whole.Pro.Sig.B, 20)
LVQ.whole.Pro.Sig.B$Overall <- NULL
SVM.whole.Pro.Sig.B$Overall <- NULL
GBM.whole.Pro.Sig.B$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.whole.Pro.Sig.B <- unlist(LVQ.whole.Pro.Sig.B)
SVM.whole.Pro.Sig.B <- unlist(SVM.whole.Pro.Sig.B)
GBM.whole.Pro.Sig.B <- unlist(GBM.whole.Pro.Sig.B)
library(VennDiagram)
venn.data.whole.B <- list(LVQ.whole.Pro.Sig.B, SVM.whole.Pro.Sig.B, GBM.whole.Pro.Sig.B)
grid.newpage()
venn.plot.whole.B <- venn.diagram(x = list(LVQ.whole.Pro.Sig.B=LVQ.whole.Pro.Sig.B, SVM.whole.Pro.Sig.B=SVM.whole.Pro.Sig.B, GBM.whole.Pro.Sig.B=GBM.whole.Pro.Sig.B),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.whole.B)
venn.intersect.whole.B <- calculate.overlap(venn.data.whole.B)
print(venn.intersect.whole.B$a5)

#LVQ C_
test.whole.lvq.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.lvq.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data.whole.C, method = "lvq",
                               trControl = test.whole.lvq.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.lvq.C.pro <- varImp(model.whole.lvq.C.pro, scale = FALSE)
print(importance.whole.lvq.C.pro)
plot(importance.whole.lvq.C.pro)

#SVM C_
test.whole.svm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.svm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data.whole.C, method = "svmRadial",
                               trControl = test.whole.svm.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.svm.C.pro <- varImp(model.whole.svm.C.pro, scale = FALSE)
print(importance.whole.svm.C.pro)
plot(importance.whole.svm.C.pro)

#GBM C_
test.whole.gbm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.gbm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data.whole.C, method = "gbm",
                               trControl = test.whole.gbm.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.gbm.C.pro <- varImp(model.whole.gbm.C.pro, scale = FALSE)
print(importance.whole.gbm.C.pro)
plot(importance.whole.gbm.C.pro)

#Common predictors for A_ progress..yes.no.
LVQ.whole.Pro.C <- importance.whole.lvq.C.pro$importance
SVM.whole.Pro.C <- importance.whole.svm.C.pro$importance
GBM.whole.Pro.C <- importance.whole.gbm.C.pro$importance
all(LVQ.whole.Pro.C$X0 == LVQ.whole.Pro.C$X1)
all(SVM.whole.Pro.C$X0 == SVM.whole.Pro.C$X1)
LVQ.whole.Pro.C$X1 <- NULL
SVM.whole.Pro.C$X1 <- NULL
colnames(LVQ.whole.Pro.C)[1] <- "Overall"
colnames(SVM.whole.Pro.C)[1] <- "Overall"
LVQ.whole.Pro.Sig.C <- setDT(LVQ.whole.Pro.C, keep.rownames = TRUE)[]
SVM.whole.Pro.Sig.C <- setDT(SVM.whole.Pro.C, keep.rownames = TRUE)[]
GBM.whole.Pro.Sig.C <- setDT(GBM.whole.Pro.C, keep.rownames = TRUE)[]
LVQ.whole.Pro.Sig.C <- LVQ.whole.Pro.Sig.C[order(-LVQ.whole.Pro.Sig.C$Overall),]
SVM.whole.Pro.Sig.C <- SVM.whole.Pro.Sig.C[order(-SVM.whole.Pro.Sig.C$Overall),]
GBM.whole.Pro.Sig.C <- GBM.whole.Pro.Sig.C[order(-GBM.whole.Pro.Sig.C$Overall),]
LVQ.whole.Pro.Sig.C <- head(LVQ.whole.Pro.Sig.C, 20)
SVM.whole.Pro.Sig.C <- head(SVM.whole.Pro.Sig.C, 20)
#Adjust GBM selection based on printed variables
GBM.whole.Pro.Sig.B <- head(GBM.whole.Pro.Sig.B, 20)
LVQ.whole.Pro.Sig.C$Overall <- NULL
SVM.whole.Pro.Sig.C$Overall <- NULL
GBM.whole.Pro.Sig.C$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.whole.Pro.Sig.C <- unlist(LVQ.whole.Pro.Sig.C)
SVM.whole.Pro.Sig.C <- unlist(SVM.whole.Pro.Sig.C)
GBM.whole.Pro.Sig.C <- unlist(GBM.whole.Pro.Sig.C)
library(VennDiagram)
venn.data.whole.C <- list(LVQ.whole.Pro.Sig.C, SVM.whole.Pro.Sig.C, GBM.whole.Pro.Sig.C)
grid.newpage()
venn.plot.whole.C <- venn.diagram(x = list(LVQ.whole.Pro.Sig.C=LVQ.whole.Pro.Sig.C, SVM.whole.Pro.Sig.C=SVM.whole.Pro.Sig.C, GBM.whole.Pro.Sig.C=GBM.whole.Pro.Sig.C),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.whole.C)
venn.intersec.wholet.C <- calculate.overlap(venn.data.whole.C)
print(venn.intersect.whole.C$a5)

#########################################################################################################################################

#With chemo progress..yes.no.
#1. Step.AIC, Correlation matrix and RandomForest not being used

#3 Learning Vector Quantization
test.chemo.lvq.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.lvq.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                      data = alt.data11, method = "lvq",
                      trControl = test.chemo.lvq.pro,
                      preProcess = c("center", "scale"),
                      tuneLength = 10,
                      na.action = na.pass)
importance.chemo.lvq.pro <- varImp(model.chemo.lvq.pro, scale = FALSE)
print(importance.chemo.lvq.pro)
plot(importance.chemo.lvq.pro)
#4 Support Vector Machine
test.chemo.svm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.svm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                      data = alt.data11, method = "svmRadial",
                      trControl = test.chemo.svm.pro,
                      preProcess = c("center", "scale"),
                      tuneLength = 10,
                      na.action = na.pass)
importance.chemo.svm.pro <- varImp(model.chemo.svm.pro, scale = FALSE)
print(importance.chemo.svm.pro)
plot(importance.chemo.svm.pro)
#5 Gradient Boosted Machine
test.chemo.gbm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.gbm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                      data = alt.data11, method = "gbm",
                      trControl = test.chemo.gbm.pro,
                      preProcess = c("center", "scale"),
                      tuneLength = 10,
                      na.action = na.pass)

importance.chemo.gbm.pro <- varImp(model.chemo.gbm.pro, scale = FALSE)
print(importance.chemo.gbm.pro)
plot(importance.chemo.gbm.pro)

#Common predictors for progress..yes.no.
LVQ.chemo.Pro <- importance.chemo.lvq.pro$importance
SVM.chemo.Pro <- importance.chemo.svm.pro$importance
GBM.chemo.Pro <- importance.chemo.gbm.pro$importance
all(LVQ.chemo.Pro$X0 == LVQ.chemo.Pro$X1)
all(SVM.chemo.Pro$X0 == SVM.chemo.Pro$X1)
LVQ.chemo.Pro$X1 <- NULL
SVM.chemo.Pro$X1 <- NULL
colnames(LVQ.chemo.Pro)[1] <- "Overall"
colnames(SVM.chemo.Pro)[1] <- "Overall"
LVQ.chemo.Pro.Sig <- setDT(LVQ.chemo.Pro, keep.rownames = TRUE)[]
SVM.chemo.Pro.Sig <- setDT(SVM.chemo.Pro, keep.rownames = TRUE)[]
GBM.chemo.Pro.Sig <- setDT(GBM.chemo.Pro, keep.rownames = TRUE)[]
LVQ.chemo.Pro.Sig <- LVQ.chemo.Pro.Sig[order(-LVQ.chemo.Pro.Sig$Overall),]
SVM.chemo.Pro.Sig <- SVM.chemo.Pro.Sig[order(-SVM.chemo.Pro.Sig$Overall),]
GBM.chemo.Pro.Sig <- GBM.chemo.Pro.Sig[order(-GBM.chemo.Pro.Sig$Overall),]
LVQ.chemo.Pro.Sig <- head(LVQ.chemo.Pro.Sig, 20)
SVM.chemo.Pro.Sig <- head(SVM.chemo.Pro.Sig, 20)
#Adjust GBM selection based on printed variables
GBM.chemo.Pro.Sig <- head(GBM.chemo.Pro.Sig, 20)
LVQ.chemo.Pro.Sig$Overall <- NULL
SVM.chemo.Pro.Sig$Overall <- NULL
GBM.chemo.Pro.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.chemo.Pro.Sig <- unlist(LVQ.chemo.Pro.Sig)
SVM.chemo.Pro.Sig <- unlist(SVM.chemo.Pro.Sig)
GBM.chemo.Pro.Sig <- unlist(GBM.chemo.Pro.Sig)

library(VennDiagram)

venn.data.chemo <- list(LVQ.chemo.Pro.chemo.Sig, SVM.chemo.Pro.Sig, GBM.chemo.Pro.Sig)

grid.newpage()

venn.plot.chemo <- venn.diagram(x = list(LVQ.chemo.Pro.Sig=LVQ.chemo.Pro.Sig, SVM.chemo.Pro.Sig=SVM.chemo.Pro.Sig, GBM.chemo.Pro.Sig=GBM.chemo.Pro.Sig),
                          filename=NULL, 
                          fill = c("red", "blue", "green"),
                          alpha = 0.50,
                          col = "transparent")

grid.draw(venn.plot.chemo)
venn.intersect.chemo <- calculate.overlap(venn.data.chemo)
print(venn.intersect.chemo$a5)

##################################################################################################################################
#Seperation for A_, B_ and C_ using progress..yes.no.
#LVQ A_
test.chemo.lvq.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.lvq.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF,
                             data = alt.data11.A, method = "lvq",
                             trControl = test.chemo.lvq.A.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.chemo.lvq.A.pro <- varImp(model.chemo.lvq.A.pro, scale = FALSE)
print(importance.chemo.lvq.A.pro)
plot(importance.chemo.lvq.A.pro)

#SVM A_
test.chemo.svm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.svm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                             data = alt.data11.A, method = "svmRadial",
                             trControl = test.chemo.svm.A.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.chemo.svm.A.pro <- varImp(model.chemo.svm.A.pro, scale = FALSE)
print(importance.chemo.svm.A.pro)
plot(importance.chemo.svm.A.pro)

#GBM A_
test.chemo.gbm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.gbm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                             data = alt.data11.A, method = "gbm",
                             trControl = test.chemo.gbm.A.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)

importance.chemo.gbm.A.pro <- varImp(model.chemo.gbm.A.pro, scale = FALSE)
print(importance.chemo.gbm.A.pro)
plot(importance.chemo.gbm.A.pro)

#Common predictors for A_ progress..yes.no.
LVQ.chemo.Pro.A <- importance.chemo.lvq.A.pro$importance
SVM.chemo.Pro.A <- importance.chemo.svm.A.pro$importance
GBM.chemo.Pro.A <- importance.chemo.gbm.A.pro$importance
all(LVQ.chemo.Pro.A$X0 == LVQ.chemo.Pro.A$X1)
all(SVM.chemo.Pro.A$X0 == SVM.chemo.Pro.A$X1)
LVQ.Pro.A$X1 <- NULL
SVM.Pro.A$X1 <- NULL
colnames(LVQ.Pro.A)[1] <- "Overall"
colnames(SVM.Pro.A)[1] <- "Overall"
LVQ.chemo.Pro.Sig <- setDT(LVQ.chemo.Pro.A, keep.rownames = TRUE)[]
SVM.chemo.Pro.Sig <- setDT(SVM.chemo.Pro.A, keep.rownames = TRUE)[]
GBM.chemo.Pro.Sig <- setDT(GBM.chemo.Pro.A, keep.rownames = TRUE)[]
LVQ.chemo.Pro.Sig <- LVQ.chemo.Pro.Sig[order(-LVQ.chemo.Pro.Sig$Overall),]
SVM.chemo.Pro.Sig <- SVM.chemo.Pro.Sig[order(-SVM.chemo.Pro.Sig$Overall),]
GBM.chemo.Pro.Sig <- GBM.chemo.Pro.Sig[order(-GBM.chemo.Pro.Sig$Overall),]
LVQ.chemo.Pro.Sig <- head(LVQ.chemo.Pro.Sig, 20)
SVM.chemo.Pro.Sig <- head(SVM.chemo.Pro.Sig, 20)
#Adjust GBM selection based on printed variables
GBM.chemo.Pro.Sig <- head(GBM.chemo.Pro.Sig, 20)
LVQ.chemo.Pro.Sig.A$Overall <- NULL
SVM.chemo.Pro.Sig.A$Overall <- NULL
GBM..chemoPro.Sig.A$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.chemo.Pro.Sig.A <- unlist(LVQ.chemo.Pro.Sig.A)
SVM.chemo.Pro.Sig.A <- unlist(SVM.chemo.Pro.Sig.A)
GBM.chemo.Pro.Sig.A <- unlist(GBM.chemo.Pro.Sig.A)
library(VennDiagram)
venn.data.chemo.A <- list(LVQ.chemo.Pro.Sig.A, SVM.chemo.Pro.Sig.A, GBM.chemo.Pro.Sig.A)
grid.newpage()
venn.plot.chemo.A <- venn.diagram(x = list(LVQ.chemo.Pro.Sig.A=LVQ.chemo.Pro.Sig.A, SVM.chemo.Pro.Sig.A=SVM.chemo.Pro.Sig.A, GBM.chemo.Pro.Sig.A=GBM.chemo.Pro.Sig.A),
                          filename=NULL, 
                          fill = c("red", "blue", "green"),
                          alpha = 0.50,
                          col = "transparent")
grid.draw(venn.plot.chemo.A)
venn.intersect.chemo.A <- calculate.overlap(venn.data.chemo.A)
print(venn.intersect.chemo.A$a5)
########### Below here isn't done yet
#LVQ B_
test.chemo.lvq.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.lvq.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data11.B, method = "lvq",
                               trControl = test.chemo.lvq.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.chemo.lvq.B.pro <- varImp(model.chemo.lvq.B.pro, scale = FALSE)
print(importance.chemo.lvq.B.pro)
plot(importance.chemo.lvq.B.pro)

#SVM B_
test.chemo.svm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.svm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data11.B, method = "svmRadial",
                               trControl = test.chemo.svm.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.chemo.svm.B.pro <- varImp(model.chemo.svm.B.pro, scale = FALSE)
print(importance.chemo.svm.B.pro)
plot(importance.chemo.svm.B.pro)

#GBM B_
test.chemo.gbm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.gbm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data11.B, method = "gbm",
                               trControl = test.chemo.gbm.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.chemo.gbm.B.pro <- varImp(model.chemo.gbm.B.pro, scale = FALSE)
print(importance.chemo.gbm.B.pro)
plot(importance.chemo.gbm.B.pro)

#Common predictors for A_ progress..yes.no.
LVQ.Pro.B <- importance.chemo.lvq.B.pro$importance
SVM.Pro.B <- importance.chemo.svm.B.pro$importance
GBM.Pro.B <- importance.chemo.gbm.B.pro$importance
all(LVQ.Pro.B$X0 == LVQ.Pro.B$X1)
all(SVM.Pro.B$X0 == SVM.Pro.B$X1)
LVQ.Pro.B$X1 <- NULL
SVM.Pro.B$X1 <- NULL
colnames(LVQ.Pro.B)[1] <- "Overall"
colnames(SVM.Pro.B)[1] <- "Overall"
LVQ.Pro.Sig.B <- subset(LVQ.Pro.B, Overall >= '0.5083', select = c("Overall"))
SVM.Pro.Sig.B <- subset(SVM.Pro.B, Overall >= '0.5083', select = c("Overall"))
GBM.Pro.Sig.B <- subset(GBM.Pro.B, Overall >= '0.02152', select = c("Overall"))
library(data.table)
LVQ.Pro.Sig.B <- setDT(LVQ.Pro.Sig.B, keep.rownames = TRUE)[]
SVM.Pro.Sig.B <- setDT(SVM.Pro.Sig.B, keep.rownames = TRUE)[]
GBM.Pro.Sig.B <- setDT(GBM.Pro.Sig.B, keep.rownames = TRUE)[]
LVQ.Pro.Sig.B$Overall <- NULL
SVM.Pro.Sig.B$Overall <- NULL
GBM.Pro.Sig.B$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.Pro.Sig.B <- unlist(LVQ.Pro.Sig.B)
SVM.Pro.Sig.B <- unlist(SVM.Pro.Sig.B)
GBM.Pro.Sig.B <- unlist(GBM.Pro.Sig.B)
library(VennDiagram)
venn.data.B <- list(LVQ.Pro.Sig.B, SVM.Pro.Sig.B, GBM.Pro.Sig.B)
grid.newpage()
venn.plot.B <- venn.diagram(x = list(LVQ.Pro.Sig.B=LVQ.Pro.Sig.B, SVM.Pro.Sig.B=SVM.Pro.Sig.B, GBM.Pro.Sig.B=GBM.Pro.Sig.B),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.B)
venn.intersect.B <- calculate.overlap(venn.data.B)
print(venn.intersect.B$a5)

#LVQ C_
test.chemo.lvq.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.lvq.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data11.C, method = "lvq",
                               trControl = test.chemo.lvq.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.chemo.lvq.C.pro <- varImp(model.chemo.lvq.C.pro, scale = FALSE)
print(importance.chemo.lvq.C.pro)
plot(importance.chemo.lvq.C.pro)

#SVM C_
test.chemo.svm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.svm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data11.C, method = "svmRadial",
                               trControl = test.chemo.svm.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.chemo.svm.C.pro <- varImp(model.chemo.svm.C.pro, scale = FALSE)
print(importance.chemo.svm.C.pro)
plot(importance.chemo.svm.C.pro)

#GBM C_
test.chemo.gbm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.chemo.gbm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data11.C, method = "gbm",
                               trControl = test.chemo.gbm.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)

importance.chemo.gbm.C.pro <- varImp(model.chemo.gbm.C.pro, scale = FALSE)
print(importance.chemo.gbm.C.pro)
plot(importance.chemo.gbm.C.pro)

#Common predictors for A_ progress..yes.no.
LVQ.Pro.C <- importance.chemo.lvq.C.pro$importance
SVM.Pro.C <- importance.chemo.svm.C.pro$importance
GBM.Pro.C <- importance.chemo.gbm.C.pro$importance
all(LVQ.Pro.C$X0 == LVQ.Pro.C$X1)
all(SVM.Pro.C$X0 == SVM.Pro.C$X1)
LVQ.Pro.C$X1 <- NULL
SVM.Pro.C$X1 <- NULL
colnames(LVQ.Pro.C)[1] <- "Overall"
colnames(SVM.Pro.C)[1] <- "Overall"
LVQ.Pro.Sig.C <- subset(LVQ.Pro.C, Overall >= '0.5082', select = c("Overall"))
SVM.Pro.Sig.C <- subset(SVM.Pro.C, Overall >= '0.5082', select = c("Overall"))
GBM.Pro.Sig.C <- subset(GBM.Pro.C, Overall >= '0.06445', select = c("Overall"))
library(data.table)
LVQ.Pro.Sig.C <- setDT(LVQ.Pro.Sig.C, keep.rownames = TRUE)[]
SVM.Pro.Sig.C <- setDT(SVM.Pro.Sig.C, keep.rownames = TRUE)[]
GBM.Pro.Sig.C <- setDT(GBM.Pro.Sig.C, keep.rownames = TRUE)[]
LVQ.Pro.Sig.C$Overall <- NULL
SVM.Pro.Sig.C$Overall <- NULL
GBM.Pro.Sig.C$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.Pro.Sig.C <- unlist(LVQ.Pro.Sig.C)
SVM.Pro.Sig.C <- unlist(SVM.Pro.Sig.C)
GBM.Pro.Sig.C <- unlist(GBM.Pro.Sig.C)
library(VennDiagram)
venn.data.C <- list(LVQ.Pro.Sig.C, SVM.Pro.Sig.C, GBM.Pro.Sig.C)
grid.newpage()
venn.plot.C <- venn.diagram(x = list(LVQ.Pro.Sig.C=LVQ.Pro.Sig.C, SVM.Pro.Sig.C=SVM.Pro.Sig.C, GBM.Pro.Sig.C=GBM.Pro.Sig.C),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.C)
venn.intersect.C <- calculate.overlap(venn.data.C)
print(venn.intersect.C$a5)

###################################################################################################

###################################################################################################################################

#Without chemo progress..yes.no.
<<<<<<< HEAD
#1. Step.AIC, Correlation matrix and RandomForest not being used 
=======

#1. Step.AIC 
#glm.nochemo.Pro <- glm(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, family = "binomial", data = alt.data00)
#step.glm.nochemo.Pro <- stepAIC(glm.nochemo.Pro, direction = "both", trace = FALSE)
#summary(step.glm.nochemo.Pro)

#2. Correlation matrix
#set.seed(126)
#data.nochemo.pro.cor <- alt.data00[,-1]
#correlationMatrix.pro <- cor(progress..yes.no., data.nochemo.pro.cor)
#print(correlationMatrix.pro2)
#May change cutoff if results are poor
#Have some missing values
#highlyCorrelated.pro2 <- findCorrelation(correlationMatrix.pro2, cutoff=0.75)
#print(highlyCorrelated.pro2)
>>>>>>> a7b3eee50f31cba08cd553bb07e79956dac5c9bb

#3 Learning Vector Quantization
test.nochemo.lvq.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.lvq.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = alt.data00, method = "lvq",
                             trControl = test.nochemo.lvq.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.nochemo.lvq.pro <- varImp(model.nochemo.lvq.pro, scale = FALSE)
print(importance.nochemo.lvq.pro)
plot(importance.nochemo.lvq.pro)

#4 Support Vector Machine
test.nochemo.svm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.svm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = alt.data00, method = "svmRadial",
                             trControl = test.nochemo.svm.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.nochemo.svm.pro <- varImp(model.nochemo.svm.pro, scale = FALSE)
print(importance.nochemo.svm.pro)
plot(importance.nochemo.svm.pro)

#5 Gradient Boosted Machine
test.nochemo.gbm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.gbm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = alt.data00, method = "gbm",
                             trControl = test.nochemo.gbm.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)

importance.nochemo.gbm.pro <- varImp(model.nochemo.gbm.pro, scale = FALSE)
print(importance.nochemo.gbm.pro)
plot(importance.nochemo.gbm.pro)


<<<<<<< HEAD
=======
#Common predictors for progress..yes.no.

LVQ.nochemo.Pro <- importance.nochemo.lvq.pro$importance
SVM.nochemo.Pro <- importance.nochemo.svm.pro$importance
GBM.nochemo.Pro <- importance.nochemo.gbm.pro$importance

all(LVQ.nochemo.Pro$X0 == LVQ.nochemo.Pro$X1)
all(SVM.nochemo.Pro$X0 == SVM.nochemo.Pro$X1)

LVQ.nochemo.Pro$X1 <- NULL
SVM.nochemo.Pro$X1 <- NULL
colnames(LVQ.nochemo.Pro)[1] <- "Overall"
colnames(SVM.nochemo.Pro)[1] <- "Overall"

LVQ.nochemo.Pro.Sig <- subset(LVQ.nochemo.Pro, Overall >= '0.5583', select = c("Overall"))
SVM.nochemo.Pro.Sig <- subset(SVM.nochemo.Pro, Overall >= '0.5417', select = c("Overall"))
GBM.nochemo.Pro.Sig <- subset(GBM.nochemo.Pro, Overall >= '0.09', select = c("Overall"))

library(data.table)

LVQ.nochemo.Pro.Sig <- setDT(LVQ.nochemo.Pro.Sig, keep.rownames = TRUE)[]
SVM.nochemo.Pro.Sig <- setDT(SVM.nochemo.Pro.Sig, keep.rownames = TRUE)[]
GBM.nochemo.Pro.Sig <- setDT(GBM.nochemo.Pro.Sig, keep.rownames = TRUE)[]

LVQ.nochemo.Pro.Sig$Overall <- NULL
SVM.nochemo.Pro.Sig$Overall <- NULL
GBM.nochemo.Pro.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.nochemo.Pro.Sig <- unlist(LVQ.nochemo.Pro.Sig)
SVM.nochemo.Pro.Sig <- unlist(SVM.nochemo.Pro.Sig)
GBM.nochemo.Pro.Sig <- unlist(GBM.nochemo.Pro.Sig)

library(VennDiagram)

venn.data.nochemo <- list(LVQ.nochemo.Pro.Sig, SVM.nochemo.Pro.Sig, GBM.nochemo.Pro.Sig)

grid.newpage()

venn.plot.nochemo <- venn.diagram(x = list(LVQ.nochemo.Pro.Sig=LVQ.nochemo.Pro.Sig, SVM.nochemo.Pro.Sig=SVM.nochemo.Pro.Sig, GBM.nochemo.Pro.Sig=GBM.nochemo.Pro.Sig),
                          filename=NULL, 
                          fill = c("red", "blue", "green"),
                          alpha = 0.50,
                          col = "transparent")

grid.draw(venn.plot.nochemo)

venn.intersect.nochemo <- calculate.overlap(venn.data.nochemo)

print(venn.intersect.nochemo$a5)


#Without chemo PFS..months.
#1. Step.AIC
glm.nochemo.PFS <- glm(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, family = "Gamma"(link=log), data = alt.data00.PFS)
step.glm.nochemo.PFS <- stepAIC(glm.nochemo.PFS, direction = "both", trace = FALSE)
summary(step.glm.nochemo.PFS)

#2 Correlation matrix
set.seed(126)
data.nochemo.PFS.cor <- alt.data00.PFS[,-1]
correlationMatrix.PFS2 <- cor(progress..yes.no., data.nochemo.PFS.cor)
print(correlationMatrix.PFS22)
#Some missing values
#May change cutoff if results are poor
highlyCorrelated.PFS2 <- findCorrelation(correlationMatrix.PFS2, cutoff=0.5)
print(highlyCorrelated.PFS2)

#3 Learning Vector Quantization
test.nochemo.lvq.PFS <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
#Wrong model type for regression
model.nochemo.lvq.PFS <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                         data = alt.data00.PFS, method = "lvq",
                         trControl = test.chemo.lvq.PFS,
                         preProcess = c("center", "scale"),
                         tuneLength = 10,
                         na.action = na.pass)
importance.nochemo,lvq.PFS <- varImp(model.nochemo.lvq.PFS, scale = FALSE)
print(importance.nochemo,lvq.PFS)
plot(importance.nochemo,lvq.PFS)

#4 Support Vector Machine
test.nochemo.svm.PFS <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.svm.PFS <- train(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = alt.data00.PFS, method = "svmRadial",
                             trControl=test.chemo.svm.PFS,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.nochemo.svm.PFS <- varImp(model.nochemo.svm.PFS, scale = FALSE)
print(importance.nochemo.svm.PFS)
plot(importance.nochemo.svm.PFS)

#5 Gradient Boosted Machine
#Remove all the columns with zero variance
alt.data00.PFS.rm <- alt.data00.PFS[, colSums(alt.data00.PFS != 0) > 0]

test.nochemo.gbm.PFS <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.gbm.PFS <- train(progress..yes.no. ~ ., 
                             data = alt.data00.PFS.rm, method = "gbm",
                             trControl = test.chemo.gbm.PFS,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.nochemo.gbm.PFS <- varImp(model.nochemo.gbm.PFS, scale = FALSE)
print(importance.nochemo.gbm.PFS)
plot(importance.nochemo.gbm.PFS)

#6 RandomForest
set.seed(86)
test.nochemo.rf.PFS <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
results.nochemo.rf.PFS <- rfe(alt.data00.PFS[,2:146], alt.data00.PFS[,1], sizes = c(2:146), rfeControl = test.nochemo.rf.PFS)
print(results.nochemo.rf.PFS)
predictors(results.nochemo.rf.PFS)
plot(results.nochemo.rf.PFS, type = c("g", "o"))
>>>>>>> a7b3eee50f31cba08cd553bb07e79956dac5c9bb
