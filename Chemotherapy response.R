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
progress..yes.no. <- as.factor(alt.data0$progress..yes.no.)
alt.data00 <- cbind(progress..yes.no., alt.data00)
#For treatment with chemotherapy 
alt.data11 <- alt.data1[,c(34:180)]
alt.data11$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data11$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(alt.data1$progress..yes.no.)
alt.data11 <- cbind(progress..yes.no., alt.data11)

#Seperating Chemo dataset into A_, B_, C_
progress..yes.no. <- as.factor(alt.data11$progress..yes.no.)
alt.data11.A <- alt.data11[,c(1:72)]
alt.data11.B <- alt.data11[,c(73:105)]
alt.data11.B <- cbind(progress..yes.no., alt.data11.B)
alt.data11.C <- alt.data11[,c(105:146)]
alt.data11.C <- cbind(progress..yes.no., alt.data11.C)

#Seperating NoChemo dataset into A_, B_, C_
progress..yes.no. <- as.factor(alt.data00$progress..yes.no.)
alt.data00.A <- alt.data00[,c(1:72)]
alt.data00.B <- alt.data00[,c(73:105)]
alt.data00.B <- cbind(progress..yes.no., alt.data00.B)
alt.data00.C <- alt.data00[,c(105:146)]
alt.data00.C <- cbind(progress..yes.no., alt.data00.C)

library(caret)
library(gbm)
library(MASS)

######################################################################################################################################
#Whole dataset
#LVQ 
test.whole.lvq.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.lvq.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                               data = alt.data.whole, method = "lvq",
                               trControl = test.whole.lvq.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.lvq.pro <- varImp(model.whole.lvq.pro, scale = FALSE)
print(importance.whole.lvq.pro)
plot(importance.whole.lvq.pro)

#SVM
test.whole.svm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.svm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data.whole, method = "svmRadial",
                               trControl = test.whole.svm.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.svm.pro <- varImp(model.whole.svm.pro, scale = FALSE)
print(importance.whole.svm.pro)
plot(importance.whole.svm.pro)

#GBM 
test.whole.gbm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.whole.gbm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data.whole, method = "gbm",
                               trControl = test.whole.gbm.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.whole.gbm.pro <- varImp(model.whole.gbm.pro, scale = FALSE)
print(importance.whole.gbm.pro)
plot(importance.whole.gbm.pro)

#Common predictors for  progress..yes.no.
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
GBM.whole.Pro.Sig <- head(GBM.whole.Pro, 13)
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
print(venn.intersect.whole$a2)

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
GBM.whole.Pro.Sig.A <- head(GBM.whole.Pro.A, 3)
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
LVQ.whole.Pro.Sig.B <- head(LVQ.whole.Pro.Sig.B, 20)
SVM.whole.Pro.Sig.B <- head(SVM.whole.Pro.Sig.B, 20)
#Adjust GBM selection based on printed variables
GBM.whole.Pro.Sig.B <- head(GBM.whole.Pro.Sig.B, 5)
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
print(venn.intersect.whole.B$a2)

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
GBM.whole.Pro.Sig.C <- head(GBM.whole.Pro.Sig.C, 6)
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
venn.intersect.whole.C <- calculate.overlap(venn.data.whole.C)
print(venn.intersect.whole.C$a5)
print(venn.intersect.whole.C$a2)

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
GBM.chemo.Pro.Sig <- head(GBM.chemo.Pro.Sig, 7)
LVQ.chemo.Pro.Sig$Overall <- NULL
SVM.chemo.Pro.Sig$Overall <- NULL
GBM.chemo.Pro.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.chemo.Pro.Sig <- unlist(LVQ.chemo.Pro.Sig)
SVM.chemo.Pro.Sig <- unlist(SVM.chemo.Pro.Sig)
GBM.chemo.Pro.Sig <- unlist(GBM.chemo.Pro.Sig)
library(VennDiagram)
venn.data.chemo <- list(LVQ.chemo.Pro.Sig, SVM.chemo.Pro.Sig, GBM.chemo.Pro.Sig)
grid.newpage()
venn.plot.chemo <- venn.diagram(x = list(LVQ.chemo.Pro.Sig=LVQ.chemo.Pro.Sig, SVM.chemo.Pro.Sig=SVM.chemo.Pro.Sig, GBM.chemo.Pro.Sig=GBM.chemo.Pro.Sig),
                          filename=NULL, 
                          fill = c("red", "blue", "green"),
                          alpha = 0.50,
                          col = "transparent")

grid.draw(venn.plot.chemo)
venn.intersect.chemo <- calculate.overlap(venn.data.chemo)
print(venn.intersect.chemo$a5)
print(venn.intersect.chemo$a2)

##################################################################################################################################
#Seperation for chemo by A_, B_ and C_ using progress..yes.no.
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
LVQ.chemo.Pro.A$X1 <- NULL
SVM.chemo.Pro.A$X1 <- NULL
colnames(LVQ.chemo.Pro.A)[1] <- "Overall"
colnames(SVM.chemo.Pro.A)[1] <- "Overall"
LVQ.chemo.Pro.Sig.A <- setDT(LVQ.chemo.Pro.A, keep.rownames = TRUE)[]
SVM.chemo.Pro.Sig.A <- setDT(SVM.chemo.Pro.A, keep.rownames = TRUE)[]
GBM.chemo.Pro.Sig.A <- setDT(GBM.chemo.Pro.A, keep.rownames = TRUE)[]
LVQ.chemo.Pro.Sig.A <- LVQ.chemo.Pro.Sig.A[order(-LVQ.chemo.Pro.Sig.A$Overall),]
SVM.chemo.Pro.Sig.A <- SVM.chemo.Pro.Sig.A[order(-SVM.chemo.Pro.Sig.A$Overall),]
GBM.chemo.Pro.Sig.A <- GBM.chemo.Pro.Sig.A[order(-GBM.chemo.Pro.Sig.A$Overall),]
LVQ.chemo.Pro.Sig.A <- head(LVQ.chemo.Pro.Sig.A, 20)
SVM.chemo.Pro.Sig.A <- head(SVM.chemo.Pro.Sig.A, 20)
#Adjust GBM selection based on printed variables
GBM.chemo.Pro.Sig.A <- head(GBM.chemo.Pro.Sig.A, 2)
LVQ.chemo.Pro.Sig.A$Overall <- NULL
SVM.chemo.Pro.Sig.A$Overall <- NULL
GBM.chemo.Pro.Sig.A$Overall <- NULL
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
print(venn.intersect.chemo.A$a2)

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
LVQ.chemo.Pro.B <- importance.chemo.lvq.B.pro$importance
SVM.chemo.Pro.B <- importance.chemo.svm.B.pro$importance
GBM.chemo.Pro.B <- importance.chemo.gbm.B.pro$importance
all(LVQ.chemo.Pro.B$X0 == LVQ.chemo.Pro.B$X1)
all(SVM.chemo.Pro.B$X0 == SVM.chemo.Pro.B$X1)
LVQ.chemo.Pro.B$X1 <- NULL
SVM.chemo.Pro.B$X1 <- NULL
colnames(LVQ.chemo.Pro.B)[1] <- "Overall"
colnames(SVM.chemo.Pro.B)[1] <- "Overall"
LVQ.chemo.Pro.Sig.B <- setDT(LVQ.chemo.Pro.B, keep.rownames = TRUE)[]
SVM.chemo.Pro.Sig.B <- setDT(SVM.chemo.Pro.B, keep.rownames = TRUE)[]
GBM.chemo.Pro.Sig.B <- setDT(GBM.chemo.Pro.B, keep.rownames = TRUE)[]
LVQ.chemo.Pro.Sig.B <- LVQ.chemo.Pro.Sig.B[order(-LVQ.chemo.Pro.Sig.B$Overall),]
SVM.chemo.Pro.Sig.B <- SVM.chemo.Pro.Sig.B[order(-SVM.chemo.Pro.Sig.B$Overall),]
GBM.chemo.Pro.Sig.B <- GBM.chemo.Pro.Sig.B[order(-GBM.chemo.Pro.Sig.B$Overall),]
LVQ.chemo.Pro.Sig.B <- head(LVQ.chemo.Pro.Sig.B, 20)
SVM.chemo.Pro.Sig.B <- head(SVM.chemo.Pro.Sig.B, 20)
#Adjust GBM selection based on printed variables
GBM.chemo.Pro.Sig.B <- head(GBM.chemo.Pro.Sig.B, 3)
LVQ.chemo.Pro.Sig.B$Overall <- NULL
SVM.chemo.Pro.Sig.B$Overall <- NULL
GBM.chemo.Pro.Sig.B$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.chemo.Pro.Sig.B <- unlist(LVQ.chemo.Pro.Sig.B)
SVM.chemo.Pro.Sig.B <- unlist(SVM.chemo.Pro.Sig.B)
GBM.chemo.Pro.Sig.B <- unlist(GBM.chemo.Pro.Sig.B)
library(VennDiagram)
venn.data.chemo.B <- list(LVQ.chemo.Pro.Sig.B, SVM.chemo.Pro.Sig.B, GBM.chemo.Pro.Sig.B)
grid.newpage()
venn.plot.chemo.B <- venn.diagram(x = list(LVQ.chemo.Pro.Sig.B=LVQ.chemo.Pro.Sig.B, SVM.chemo.Pro.Sig.B=SVM.chemo.Pro.Sig.B, GBM.chemo.Pro.Sig.B=GBM.chemo.Pro.Sig.B),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.chemo.B)
venn.intersect.chemo.B <- calculate.overlap(venn.data.chemo.B)
print(venn.intersect.chemo.B$a5)
print(venn.intersect.chemo.B$a2)

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
LVQ.chemo.Pro.C <- importance.chemo.lvq.C.pro$importance
SVM.chemo.Pro.C <- importance.chemo.svm.C.pro$importance
GBM.chemo.Pro.C <- importance.chemo.gbm.C.pro$importance
all(LVQ.chemo.Pro.C$X0 == LVQ.chemo.Pro.C$X1)
all(SVM.chemo.Pro.C$X0 == SVM.chemo.Pro.C$X1)
LVQ.chemo.Pro.C$X1 <- NULL
SVM.chemo.Pro.C$X1 <- NULL
colnames(LVQ.chemo.Pro.C)[1] <- "Overall"
colnames(SVM.chemo.Pro.C)[1] <- "Overall"
LVQ.chemo.Pro.Sig.C <- setDT(LVQ.chemo.Pro.C, keep.rownames = TRUE)[]
SVM.chemo.Pro.Sig.C <- setDT(SVM.chemo.Pro.C, keep.rownames = TRUE)[]
GBM.chemo.Pro.Sig.C <- setDT(GBM.chemo.Pro.C, keep.rownames = TRUE)[]
LVQ.chemo.Pro.Sig.C <- LVQ.chemo.Pro.Sig.C[order(-LVQ.chemo.Pro.Sig.C$Overall),]
SVM.chemo.Pro.Sig.C <- SVM.chemo.Pro.Sig.C[order(-SVM.chemo.Pro.Sig.C$Overall),]
GBM.chemo.Pro.Sig.C <- GBM.chemo.Pro.Sig.C[order(-GBM.chemo.Pro.Sig.C$Overall),]
LVQ.chemo.Pro.Sig.C <- head(LVQ.chemo.Pro.Sig.C, 20)
SVM.chemo.Pro.Sig.C <- head(SVM.chemo.Pro.Sig.C, 20)
#Adjust GBM selection based on printed variables
GBM.chemo.Pro.Sig.C <- head(GBM.chemo.Pro.Sig.C, 4)
LVQ.chemo.Pro.Sig.C$Overall <- NULL
SVM.chemo.Pro.Sig.C$Overall <- NULL
GBM.chemo.Pro.Sig.C$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.chemo.Pro.Sig.C <- unlist(LVQ.chemo.Pro.Sig.C)
SVM.chemo.Pro.Sig.C <- unlist(SVM.chemo.Pro.Sig.C)
GBM.chemo.Pro.Sig.C <- unlist(GBM.chemo.Pro.Sig.C)
library(VennDiagram)
venn.data.chemo.C <- list(LVQ.chemo.Pro.Sig.C, SVM.chemo.Pro.Sig.C, GBM.chemo.Pro.Sig.C)
grid.newpage()
venn.plot.chemo.C <- venn.diagram(x = list(LVQ.chemo.Pro.Sig.C=LVQ.chemo.Pro.Sig.C, SVM.chemo.Pro.Sig.C=SVM.chemo.Pro.Sig.C, GBM.chemo.Pro.Sig.C=GBM.chemo.Pro.Sig.C),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.chemo.C)
venn.intersect.chemo.C <- calculate.overlap(venn.data.chemo.C)
print(venn.intersect.chemo.C$a5)
print(venn.intersect.chemo.C$a2)

###################################################################################################################################
######Removed all the GBM functions because they wouldn't work (possibly due to small datasets)
#Without chemo progress..yes.no.
#1. Step.AIC, Correlation matrix and RandomForest not being used 

#According to author "default LVQ grid creation code doesn't work on a such a small data set" so specified using his grid 
#3 Learning Vector Quantization
test.nochemo.lvq.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.lvq.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = alt.data00, method = "lvq",
                             trControl = test.nochemo.lvq.pro,
                             preProcess = c("center", "scale"),
                             tuneGrid = data.frame(size = 3, k = 1:2),
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

# GBM does not work 
#5 Gradient Boosted Machine
#test.nochemo.gbm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.nochemo.gbm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
#                             data = alt.data00, method = "gbm",
#                             trControl = test.nochemo.gbm.pro,
#                             preProcess = c("center", "scale"),
#                             tuneLength = 10,
#                             na.action = na.pass)
#importance.nochemo.gbm.pro <- varImp(model.nochemo.gbm.pro, scale = FALSE)
#print(importance.nochemo.gbm.pro)
#plot(importance.nochemo.gbm.pro)

#Common predictors for progress..yes.no.
LVQ.nochemo.Pro <- importance.nochemo.lvq.pro$importance
SVM.nochemo.Pro <- importance.nochemo.svm.pro$importance
#GBM.nochemo.Pro <- importance.nochemo.gbm.pro$importance
all(LVQ.nochemo.Pro$X0 == LVQ.nochemo.Pro$X1)
all(SVM.nochemo.Pro$X0 == SVM.nochemo.Pro$X1)
LVQ.nochemo.Pro$X1 <- NULL
SVM.nochemo.Pro$X1 <- NULL
colnames(LVQ.nochemo.Pro)[1] <- "Overall"
colnames(SVM.nochemo.Pro)[1] <- "Overall"
LVQ.nochemo.Pro.Sig <- setDT(LVQ.nochemo.Pro, keep.rownames = TRUE)[]
SVM.nochemo.Pro.Sig <- setDT(SVM.nochemo.Pro, keep.rownames = TRUE)[]
#GBM.nochemo.Pro.Sig <- setDT(GBM.nochemo.Pro, keep.rownames = TRUE)[]
LVQ.nochemo.Pro.Sig <- LVQ.nochemo.Pro.Sig[order(-LVQ.nochemo.Pro.Sig$Overall),]
SVM.nochemo.Pro.Sig <- SVM.nochemo.Pro.Sig[order(-SVM.nochemo.Pro.Sig$Overall),]
#GBM.nochemo.Pro.Sig <- GBM.nochemo.Pro.Sig[order(-GBM.nochemo.Pro.Sig$Overall),]
LVQ.nochemo.Pro.Sig <- head(LVQ.nochemo.Pro.Sig, 20)
SVM.nochemo.Pro.Sig <- head(SVM.nochemo.Pro.Sig, 20)
#Adjust GBM selection based on printed variables
#GBM.nochemo.Pro.Sig <- head(GBM.nochemo.Pro.Sig, 20)
LVQ.nochemo.Pro.Sig$Overall <- NULL
SVM.nochemo.Pro.Sig$Overall <- NULL
#GBM.nochemo.Pro.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.nochemo.Pro.Sig <- unlist(LVQ.nochemo.Pro.Sig)
SVM.nochemo.Pro.Sig <- unlist(SVM.nochemo.Pro.Sig)
#GBM.nochemo.Pro.Sig <- unlist(GBM.nochemo.Pro.Sig)
library(VennDiagram)
venn.data.nochemo <- list(LVQ.nochemo.Pro.Sig, SVM.nochemo.Pro.Sig)
grid.newpage()
venn.plot.nochemo <- venn.diagram(x = list(LVQ.nochemo.Pro.Sig=LVQ.nochemo.Pro.Sig, SVM.nochemo.Pro.Sig=SVM.nochemo.Pro.Sig),
                          filename=NULL, 
                          fill = c("red", "blue"),
                          alpha = 0.50,
                          col = "transparent")
grid.draw(venn.plot.nochemo)
venn.intersect.nochemo <- calculate.overlap(venn.data.nochemo)
print(venn.intersect.nochemo$a2)

##################################################################################################################################

#Seperation for nochemo by A_, B_ and C_ using progress..yes.no.
#LVQ A_
test.nochemo.lvq.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.lvq.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF,
                               data = alt.data00.A, method = "lvq",
                               trControl = test.nochemo.lvq.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.nochemo.lvq.A.pro <- varImp(model.nochemo.lvq.A.pro, scale = FALSE)
print(importance.nochemo.lvq.A.pro)
plot(importance.nochemo.lvq.A.pro)

#SVM A_
test.nochemo.svm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.svm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                               data = alt.data00.A, method = "svmRadial",
                               trControl = test.nochemo.svm.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.nochemo.svm.A.pro <- varImp(model.nochemo.svm.A.pro, scale = FALSE)
print(importance.nochemo.svm.A.pro)
plot(importance.nochemo.svm.A.pro)

#GBM A_
#test.nochemo.gbm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.nochemo.gbm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
#                               data = alt.data00.A, method = "gbm",
#                               trControl = test.nochemo.gbm.A.pro,
#                               preProcess = c("center", "scale"),
#                               tuneLength = 10,
#                               na.action = na.pass)
#importance.nochemo.gbm.A.pro <- varImp(model.nochemo.gbm.A.pro, scale = FALSE)
#print(importance.nochemo.gbm.A.pro)
#plot(importance.nochemo.gbm.A.pro)

#Common predictors for A_ progress..yes.no.
LVQ.nochemo.Pro.A <- importance.nochemo.lvq.A.pro$importance
SVM.nochemo.Pro.A <- importance.nochemo.svm.A.pro$importance
#GBM.nochemo.Pro.A <- importance.nochemo.gbm.A.pro$importance
all(LVQ.nochemo.Pro.A$X0 == LVQ.nochemo.Pro.A$X1)
all(SVM.nochemo.Pro.A$X0 == SVM.nochemo.Pro.A$X1)
LVQ.nochemo.Pro.A$X1 <- NULL
SVM.nochemo.Pro.A$X1 <- NULL
colnames(LVQ.nochemo.Pro.A)[1] <- "Overall"
colnames(SVM.nochemo.Pro.A)[1] <- "Overall"
LVQ.nochemo.Pro.Sig.A <- setDT(LVQ.nochemo.Pro.A, keep.rownames = TRUE)[]
SVM.nochemo.Pro.Sig.A <- setDT(SVM.nochemo.Pro.A, keep.rownames = TRUE)[]
#GBM.nochemo.Pro.Sig.A <- setDT(GBM.nochemo.Pro.A, keep.rownames = TRUE)[]
LVQ.nochemo.Pro.Sig.A <- LVQ.nochemo.Pro.Sig.A[order(-LVQ.nochemo.Pro.Sig.A$Overall),]
SVM.nochemo.Pro.Sig.A <- SVM.nochemo.Pro.Sig.A[order(-SVM.nochemo.Pro.Sig.A$Overall),]
#GBM.nochemo.Pro.Sig.A <- GBM.nochemo.Pro.Sig.A[order(-GBM.nochemo.Pro.Sig.A$Overall),]
LVQ.nochemo.Pro.Sig.A <- head(LVQ.nochemo.Pro.Sig.A, 20)
SVM.nochemo.Pro.Sig.A <- head(SVM.nochemo.Pro.Sig.A, 20)
#Adjust GBM selection based on printed variables
#GBM.nochemo.Pro.Sig.A <- head(GBM.nochemo.Pro.Sig.A, 20)
LVQ.nochemo.Pro.Sig.A$Overall <- NULL
SVM.nochemo.Pro.Sig.A$Overall <- NULL
#GBM.nochemo.Pro.Sig.A$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.nochemo.Pro.Sig.A <- unlist(LVQ.nochemo.Pro.Sig.A)
SVM.nochemo.Pro.Sig.A <- unlist(SVM.nochemo.Pro.Sig.A)
#GBM.nochemo.Pro.Sig.A <- unlist(GBM.nochemo.Pro.Sig.A)
library(VennDiagram)
venn.data.nochemo.A <- list(LVQ.nochemo.Pro.Sig.A, SVM.nochemo.Pro.Sig.A)
grid.newpage()
venn.plot.nochemo.A <- venn.diagram(x = list(LVQ.nochemo.Pro.Sig.A=LVQ.nochemo.Pro.Sig.A, SVM.nochemo.Pro.Sig.A=SVM.nochemo.Pro.Sig.A),
                                  filename=NULL, 
                                  fill = c("red", "blue"),
                                  alpha = 0.50,
                                  col = "transparent")
grid.draw(venn.plot.nochemo.A)
venn.intersect.nochemo.A <- calculate.overlap(venn.data.nochemo.A)
print(venn.intersect.nochemo.A$a5)

#LVQ B_
test.nochemo.lvq.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.lvq.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data00.B, method = "lvq",
                               trControl = test.nochemo.lvq.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.nochemo.lvq.B.pro <- varImp(model.nochemo.lvq.B.pro, scale = FALSE)
print(importance.nochemo.lvq.B.pro)
plot(importance.nochemo.lvq.B.pro)

#SVM B_
test.nochemo.svm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.svm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = alt.data00.B, method = "svmRadial",
                               trControl = test.nochemo.svm.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.nochemo.svm.B.pro <- varImp(model.nochemo.svm.B.pro, scale = FALSE)
print(importance.nochemo.svm.B.pro)
plot(importance.nochemo.svm.B.pro)

#GBM B_
#test.nochemo.gbm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.nochemo.gbm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
#                               data = alt.data00.B, method = "gbm",
#                               trControl = test.nochemo.gbm.B.pro,
#                               preProcess = c("center", "scale"),
#                               tuneLength = 10,
#                               na.action = na.pass)
#importance.nochemo.gbm.B.pro <- varImp(model.nochemo.gbm.B.pro, scale = FALSE)
#print(importance.nochemo.gbm.B.pro)
#plot(importance.nochemo.gbm.B.pro)

#Common predictors for A_ progress..yes.no.
LVQ.nochemo.Pro.B <- importance.nochemo.lvq.B.pro$importance
SVM.nochemo.Pro.B <- importance.nochemo.svm.B.pro$importance
#GBM.nochemo.Pro.B <- importance.nochemo.gbm.B.pro$importance
all(LVQ.nochemo.Pro.B$X0 == LVQ.nochemo.Pro.B$X1)
all(SVM.nochemo.Pro.B$X0 == SVM.nochemo.Pro.B$X1)
LVQ.nochemo.Pro.B$X1 <- NULL
SVM.nochemo.Pro.B$X1 <- NULL
colnames(LVQ.nochemo.Pro.B)[1] <- "Overall"
colnames(SVM.nochemo.Pro.B)[1] <- "Overall"
LVQ.nochemo.Pro.Sig.B <- setDT(LVQ.nochemo.Pro.B, keep.rownames = TRUE)[]
SVM.nochemo.Pro.Sig.B <- setDT(SVM.nochemo.Pro.B, keep.rownames = TRUE)[]
#GBM.nochemo.Pro.Sig.B <- setDT(GBM.nochemo.Pro.B, keep.rownames = TRUE)[]
LVQ.nochemo.Pro.Sig.B <- LVQ.nochemo.Pro.Sig.B[order(-LVQ.nochemo.Pro.Sig.B$Overall),]
SVM.nochemo.Pro.Sig.B <- SVM.nochemo.Pro.Sig.B[order(-SVM.nochemo.Pro.Sig.B$Overall),]
#GBM.nochemo.Pro.Sig.B <- GBM.nochemo.Pro.Sig.B[order(-GBM.nochemo.Pro.Sig.B$Overall),]
LVQ.nochemo.Pro.Sig.B <- head(LVQ.nochemo.Pro.Sig.B, 20)
SVM.nochemo.Pro.Sig.B <- head(SVM.nochemo.Pro.Sig.B, 20)
#Adjust GBM selection based on printed variables
#GBM.nochemo.Pro.Sig.B <- head(GBM.nochemo.Pro.Sig.B, 20)
LVQ.nochemo.Pro.Sig.B$Overall <- NULL
SVM.nochemo.Pro.Sig.B$Overall <- NULL
#GBM.nochemo.Pro.Sig.B$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.nochemo.Pro.Sig.B <- unlist(LVQ.nochemo.Pro.Sig.B)
SVM.nochemo.Pro.Sig.B <- unlist(SVM.nochemo.Pro.Sig.B)
#GBM.nochemo.Pro.Sig.B <- unlist(GBM.nochemo.Pro.Sig.B)
library(VennDiagram)
venn.data.nochemo.B <- list(LVQ.nochemo.Pro.Sig.B, SVM.nochemo.Pro.Sig.B)
grid.newpage()
venn.plot.nochemo.B <- venn.diagram(x = list(LVQ.nochemo.Pro.Sig.B=LVQ.nochemo.Pro.Sig.B, SVM.nochemo.Pro.Sig.B=SVM.nochemo.Pro.Sig.B),
                            filename=NULL, 
                            fill = c("red", "blue"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.nochemo.B)
venn.intersect.nochemo.B <- calculate.overlap(venn.data.nochemo.B)
print(venn.intersect.nochemo.B$a5)

#LVQ C_
test.nochemo.lvq.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.lvq.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data00.C, method = "lvq",
                               trControl = test.nochemo.lvq.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.nochemo.lvq.C.pro <- varImp(model.nochemo.lvq.C.pro, scale = FALSE)
print(importance.nochemo.lvq.C.pro)
plot(importance.nochemo.lvq.C.pro)

#SVM C_
test.nochemo.svm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.nochemo.svm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = alt.data00.C, method = "svmRadial",
                               trControl = test.nochemo.svm.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.nochemo.svm.C.pro <- varImp(model.nochemo.svm.C.pro, scale = FALSE)
print(importance.nochemo.svm.C.pro)
plot(importance.nochemo.svm.C.pro)

#GBM C_
#test.nochemo.gbm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.nochemo.gbm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
#                               data = alt.data00.C, method = "gbm",
#                               trControl = test.nochemo.gbm.C.pro,
#                               preProcess = c("center", "scale"),
#                               tuneLength = 10,
#                               na.action = na.pass)
#importance.nochemo.gbm.C.pro <- varImp(model.nochemo.gbm.C.pro, scale = FALSE)
#print(importance.nochemo.gbm.C.pro)
#plot(importance.nochemo.gbm.C.pro)

#Common predictors for A_ progress..yes.no.
LVQ.nochemo.Pro.C <- importance.nochemo.lvq.C.pro$importance
SVM.nochemo.Pro.C <- importance.nochemo.svm.C.pro$importance
#GBM.nochemo.Pro.C <- importance.nochemo.gbm.C.pro$importance
all(LVQ.nochemo.Pro.C$X0 == LVQ.nochemo.Pro.C$X1)
all(SVM.nochemo.Pro.C$X0 == SVM.nochemo.Pro.C$X1)
LVQ.nochemo.Pro.C$X1 <- NULL
SVM.nochemo.Pro.C$X1 <- NULL
colnames(LVQ.nochemo.Pro.C)[1] <- "Overall"
colnames(SVM.nochemo.Pro.C)[1] <- "Overall"
LVQ.nochemo.Pro.Sig.C <- setDT(LVQ.nochemo.Pro.C, keep.rownames = TRUE)[]
SVM.nochemo.Pro.Sig.C <- setDT(SVM.nochemo.Pro.C, keep.rownames = TRUE)[]
#GBM.nochemo.Pro.Sig.C <- setDT(GBM.nochemo.Pro.C, keep.rownames = TRUE)[]
LVQ.nochemo.Pro.Sig.C <- LVQ.nochemo.Pro.Sig.C[order(-LVQ.nochemo.Pro.Sig.C$Overall),]
SVM.nochemo.Pro.Sig.C <- SVM.nochemo.Pro.Sig.C[order(-SVM.nochemo.Pro.Sig.C$Overall),]
#GBM.nochemo.Pro.Sig.C <- GBM.nochemo.Pro.Sig.C[order(-GBM.nochemo.Pro.Sig.C$Overall),]
LVQ.nochemo.Pro.Sig.C <- head(LVQ.nochemo.Pro.Sig.C, 20)
SVM.nochemo.Pro.Sig.C <- head(SVM.nochemo.Pro.Sig.C, 20)
#Adjust GBM selection based on printed variables
#GBM.nochemo.Pro.Sig.C <- head(GBM.nochemo.Pro.Sig.C, 20)
LVQ.nochemo.Pro.Sig.C$Overall <- NULL
SVM.nochemo.Pro.Sig.C$Overall <- NULL
#GBM.nochemo.Pro.Sig.C$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.nochemo.Pro.Sig.C <- unlist(LVQ.nochemo.Pro.Sig.C)
SVM.nochemo.Pro.Sig.C <- unlist(SVM.nochemo.Pro.Sig.C)
#GBM.nochemo.Pro.Sig.C <- unlist(GBM.nochemo.Pro.Sig.C)
library(VennDiagram)
venn.data.nochemo.C <- list(LVQ.nochemo.Pro.Sig.C, SVM.nochemo.Pro.Sig.C)
grid.newpage()
venn.plot.nochemo.C <- venn.diagram(x = list(LVQ.nochemo.Pro.Sig.C=LVQ.nochemo.Pro.Sig.C, SVM.nochemo.Pro.Sig.C=SVM.nochemo.Pro.Sig.C),
                                  filename=NULL, 
                                  fill = c("red", "blue"),
                                  alpha = 0.50,
                                  col = "transparent")
grid.draw(venn.plot.nochemo.C)
venn.intersect.nochemo.C <- calculate.overlap(venn.data.nochemo.C)
print(venn.intersect.nochemo.C$a5)

#########################################################################################################################################

#Split df by value in the Mitotane column
mitotane.data <- split(mydata, mydata$Mitotane..pall.or.adj.)

mitotane.pos <- rbind(mitotane.data$adj, mitotane.data$`adj (with chemo)`, mitotane.data$`adj+pall`, mitotane.data$pall)
mitotane.neg <- rbind(mitotane.data$`0`)

mitotane.pos.data.whole <- mitotane.pos[,c(34:180)]
mitotane.pos.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.pos.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(mitotane.pos$progress..yes.no.)
mitotane.pos.data.whole <- cbind(progress..yes.no., mitotane.pos.data.whole)

mitotane.neg.data.whole <- mitotane.neg[,c(34:180)]
mitotane.neg.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.neg.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(mitotane.neg$progress..yes.no.)
mitotane.neg.data.whole <- cbind(progress..yes.no., mitotane.neg.data.whole)

#Seperating Mitotane dataset into A_, B_, C_
progress..yes.no. <- as.factor(mitotane.pos.data.whole$progress..yes.no.)
mitotane.pos.data.A <- mitotane.pos.data.whole[,c(1:72)]
mitotane.pos.data.B <- mitotane.pos.data.whole[,c(73:105)]
mitotane.pos.data.B <- cbind(progress..yes.no., mitotane.pos.data.B)
mitotane.pos.data.C <- mitotane.pos.data.whole[,c(105:146)]
mitotane.pos.data.C <- cbind(progress..yes.no., mitotane.pos.data.C)

#Seperating NoMitotane dataset into A_, B_, C_
progress..yes.no. <- as.factor(mitotane.neg.data.whole$progress..yes.no.)
mitotane.neg.data.A <- mitotane.neg.data.whole[,c(1:72)]
mitotane.neg.data.B <- mitotane.neg.data.whole[,c(73:105)]
mitotane.neg.data.B <- cbind(progress..yes.no., mitotane.neg.data.B)
mitotane.neg.data.C <- mitotane.neg.data.whole[,c(105:146)]
mitotane.neg.data.C <- cbind(progress..yes.no., mitotane.neg.data.C)

#Mitotane positive treatment 
#LVQ 
test.mitotane.lvq.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mitotane.lvq.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                             data = mitotane.pos.data.whole, method = "lvq",
                             trControl = test.mitotane.lvq.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.mitotane.lvq.pro <- varImp(model.mitotane.lvq.pro, scale = FALSE)
print(importance.mitotane.lvq.pro)
plot(importance.mitotane.lvq.pro)

#SVM 
test.mitotane.svm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mitotane.svm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = mitotane.pos.data.whole, method = "svmRadial",
                             trControl = test.mitotane.svm.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.mitotane.svm.pro <- varImp(model.mitotane.svm.pro, scale = FALSE)
print(importance.mitotane.svm.pro)
plot(importance.mitotane.svm.pro)

#GBM 
test.mitotane.gbm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mitotane.gbm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = mitotane.pos.data.whole, method = "gbm",
                             trControl = test.mitotane.gbm.pro,
                             preProcess = c("center", "scale"),
                             tuneLength = 10,
                             na.action = na.pass)
importance.mitotane.gbm.pro <- varImp(model.mitotane.gbm.pro, scale = FALSE)
print(importance.mitotane.gbm.pro)
plot(importance.mitotane.gbm.pro)

#Common predictors for  progress..yes.no.
LVQ.mitotane.Pro <- importance.mitotane.lvq.pro$importance
SVM.mitotane.Pro <- importance.mitotane.svm.pro$importance
GBM.mitotane.Pro <- importance.mitotane.gbm.pro$importance
all(LVQ.mitotane.Pro$X0 == LVQ.mitotane.Pro$X1)
all(SVM.mitotane.Pro$X0 == SVM.mitotane.Pro$X1)
LVQ.mitotane.Pro$X1 <- NULL
SVM.mitotane.Pro$X1 <- NULL
colnames(LVQ.mitotane.Pro)[1] <- "Overall"
colnames(SVM.mitotane.Pro)[1] <- "Overall"
LVQ.mitotane.Pro <- setDT(LVQ.mitotane.Pro, keep.rownames = TRUE)[]
SVM.mitotane.Pro <- setDT(SVM.mitotane.Pro, keep.rownames = TRUE)[]
GBM.mitotane.Pro <- setDT(GBM.mitotane.Pro, keep.rownames = TRUE)[]
LVQ.mitotane.Pro <- LVQ.mitotane.Pro[order(-LVQ.mitotane.Pro$Overall),]
SVM.mitotane.Pro <- SVM.mitotane.Pro[order(-SVM.mitotane.Pro$Overall),]
GBM.mitotane.Pro <- GBM.mitotane.Pro[order(-GBM.mitotane.Pro$Overall),]
LVQ.mitotane.Pro.Sig <- head(LVQ.mitotane.Pro, 20)
SVM.mitotane.Pro.Sig <- head(SVM.mitotane.Pro, 20)
#Adjust GBM selection based on printed variables
GBM.mitotane.Pro.Sig <- head(GBM.mitotane.Pro, 9)
LVQ.mitotane.Pro.Sig$Overall <- NULL
SVM.mitotane.Pro.Sig$Overall <- NULL
GBM.mitotane.Pro.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.mitotane.Pro.Sig <- unlist(LVQ.mitotane.Pro.Sig)
SVM.mitotane.Pro.Sig <- unlist(SVM.mitotane.Pro.Sig)
GBM.mitotane.Pro.Sig <- unlist(GBM.mitotane.Pro.Sig)
library(VennDiagram)
venn.data.mitotane.pos <- list(LVQ.mitotane.Pro.Sig, SVM.mitotane.Pro.Sig, GBM.mitotane.Pro.Sig)
grid.newpage()
venn.plot.mitotane.pos <- venn.diagram(x = list(LVQ.mitotane.Pro.Sig=LVQ.mitotane.Pro.Sig, SVM.mitotane.Pro.Sig=SVM.mitotane.Pro.Sig, GBM.mitotane.Pro.Sig=GBM.mitotane.Pro.Sig),
                                filename=NULL, 
                                fill = c("red", "blue", "green"),
                                alpha = 0.50,
                                col = "transparent")
grid.draw(venn.plot.mitotane.pos)
venn.intersect.mitotane.pos <- calculate.overlap(venn.data.mitotane.pos)
print(venn.intersect.mitotane.pos$a5)
print(venn.intersect.mitotane.pos$a2)



#Seperation for mitotane by A_, B_ and C_ using progress..yes.no.
#LVQ A_
test.mito.lvq.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.lvq.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF,
                               data = mitotane.pos.data.A, method = "lvq",
                               trControl = test.mito.lvq.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.lvq.A.pro <- varImp(model.mito.lvq.A.pro, scale = FALSE)
print(importance.mito.lvq.A.pro)
plot(importance.mito.lvq.A.pro)

#SVM A_
test.mito.svm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.svm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                               data = mitotane.pos.data.A, method = "svmRadial",
                               trControl = test.mito.svm.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.svm.A.pro <- varImp(model.mito.svm.A.pro, scale = FALSE)
print(importance.mito.svm.A.pro)
plot(importance.mito.svm.A.pro)

#GBM A_
test.mito.gbm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.gbm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                               data = mitotane.pos.data.A, method = "gbm",
                               trControl = test.mito.gbm.A.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.gbm.A.pro <- varImp(model.mito.gbm.A.pro, scale = FALSE)
print(importance.mito.gbm.A.pro)
plot(importance.mito.gbm.A.pro)

#Common predictors for A_ progress..yes.no.
LVQ.mito.Pro.A <- importance.mito.lvq.A.pro$importance
SVM.mito.Pro.A <- importance.mito.svm.A.pro$importance
GBM.mito.Pro.A <- importance.mito.gbm.A.pro$importance
all(LVQ.mito.Pro.A$X0 == LVQ.mito.Pro.A$X1)
all(SVM.mito.Pro.A$X0 == SVM.mito.Pro.A$X1)
LVQ.mito.Pro.A$X1 <- NULL
SVM.mito.Pro.A$X1 <- NULL
colnames(LVQ.mito.Pro.A)[1] <- "Overall"
colnames(SVM.mito.Pro.A)[1] <- "Overall"
LVQ.mito.Pro.Sig.A <- setDT(LVQ.mito.Pro.A, keep.rownames = TRUE)[]
SVM.mito.Pro.Sig.A <- setDT(SVM.mito.Pro.A, keep.rownames = TRUE)[]
GBM.mito.Pro.Sig.A <- setDT(GBM.mito.Pro.A, keep.rownames = TRUE)[]
LVQ.mito.Pro.Sig.A <- LVQ.mito.Pro.Sig.A[order(-LVQ.mito.Pro.Sig.A$Overall),]
SVM.mito.Pro.Sig.A <- SVM.mito.Pro.Sig.A[order(-SVM.mito.Pro.Sig.A$Overall),]
GBM.mito.Pro.Sig.A <- GBM.mito.Pro.Sig.A[order(-GBM.mito.Pro.Sig.A$Overall),]
LVQ.mito.Pro.Sig.A <- head(LVQ.mito.Pro.Sig.A, 20)
SVM.mito.Pro.Sig.A <- head(SVM.mito.Pro.Sig.A, 20)
#Adjust GBM selection based on printed variables
GBM.mito.Pro.Sig.A <- head(GBM.mito.Pro.Sig.A, 3)
LVQ.mito.Pro.Sig.A$Overall <- NULL
SVM.mito.Pro.Sig.A$Overall <- NULL
GBM.mito.Pro.Sig.A$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.mito.Pro.Sig.A <- unlist(LVQ.mito.Pro.Sig.A)
SVM.mito.Pro.Sig.A <- unlist(SVM.mito.Pro.Sig.A)
GBM.mito.Pro.Sig.A <- unlist(GBM.mito.Pro.Sig.A)
library(VennDiagram)
venn.data.mito.A <- list(LVQ.mito.Pro.Sig.A, SVM.mito.Pro.Sig.A, GBM.mito.Pro.Sig.A)
grid.newpage()
venn.plot.mito.A <- venn.diagram(x = list(LVQ.mito.Pro.Sig.A=LVQ.mito.Pro.Sig.A, SVM.mito.Pro.Sig.A=SVM.mito.Pro.Sig.A, GBM.mito.Pro.Sig.A=GBM.mito.Pro.Sig.A),
                                  filename=NULL, 
                                  fill = c("red", "blue", "green"),
                                  alpha = 0.50,
                                  col = "transparent")
grid.draw(venn.plot.mito.A)
venn.intersect.mito.A <- calculate.overlap(venn.data.mito.A)
print(venn.intersect.mito.A$a5)
print(venn.intersect.mito.A$a2)


#LVQ B_
test.mito.lvq.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.lvq.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = mitotane.pos.data.B, method = "lvq",
                               trControl = test.mito.lvq.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.lvq.B.pro <- varImp(model.mito.lvq.B.pro, scale = FALSE)
print(importance.mito.lvq.B.pro)
plot(importance.mito.lvq.B.pro)

#SVM B_
test.mito.svm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.svm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = mitotane.pos.data.B, method = "svmRadial",
                               trControl = test.mito.svm.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.svm.B.pro <- varImp(model.mito.svm.B.pro, scale = FALSE)
print(importance.mito.svm.B.pro)
plot(importance.mito.svm.B.pro)

#GBM B_
test.mito.gbm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.gbm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                               data = mitotane.pos.data.B, method = "gbm",
                               trControl = test.mito.gbm.B.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.gbm.B.pro <- varImp(model.mito.gbm.B.pro, scale = FALSE)
print(importance.mito.gbm.B.pro)
plot(importance.mito.gbm.B.pro)

#Common predictors for A_ progress..yes.no.
LVQ.mito.Pro.B <- importance.mito.lvq.B.pro$importance
SVM.mito.Pro.B <- importance.mito.svm.B.pro$importance
GBM.mito.Pro.B <- importance.mito.gbm.B.pro$importance
all(LVQ.mito.Pro.B$X0 == LVQ.mito.Pro.B$X1)
all(SVM.mito.Pro.B$X0 == SVM.mito.Pro.B$X1)
LVQ.mito.Pro.B$X1 <- NULL
SVM.mito.Pro.B$X1 <- NULL
colnames(LVQ.mito.Pro.B)[1] <- "Overall"
colnames(SVM.mito.Pro.B)[1] <- "Overall"
LVQ.mito.Pro.Sig.B <- setDT(LVQ.mito.Pro.B, keep.rownames = TRUE)[]
SVM.mito.Pro.Sig.B <- setDT(SVM.mito.Pro.B, keep.rownames = TRUE)[]
GBM.mito.Pro.Sig.B <- setDT(GBM.mito.Pro.B, keep.rownames = TRUE)[]
LVQ.mito.Pro.Sig.B <- LVQ.mito.Pro.Sig.B[order(-LVQ.mito.Pro.Sig.B$Overall),]
SVM.mito.Pro.Sig.B <- SVM.mito.Pro.Sig.B[order(-SVM.mito.Pro.Sig.B$Overall),]
GBM.mito.Pro.Sig.B <- GBM.mito.Pro.Sig.B[order(-GBM.mito.Pro.Sig.B$Overall),]
LVQ.mito.Pro.Sig.B <- head(LVQ.mito.Pro.Sig.B, 20)
SVM.mito.Pro.Sig.B <- head(SVM.mito.Pro.Sig.B, 20)
#Adjust GBM selection based on printed variables
GBM.mito.Pro.Sig.B <- head(GBM.mito.Pro.Sig.B, 3)
LVQ.mito.Pro.Sig.B$Overall <- NULL
SVM.mito.Pro.Sig.B$Overall <- NULL
GBM.mito.Pro.Sig.B$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.mito.Pro.Sig.B <- unlist(LVQ.mito.Pro.Sig.B)
SVM.mito.Pro.Sig.B <- unlist(SVM.mito.Pro.Sig.B)
GBM.mito.Pro.Sig.B <- unlist(GBM.mito.Pro.Sig.B)
library(VennDiagram)
venn.data.mito.B <- list(LVQ.mito.Pro.Sig.B, SVM.mito.Pro.Sig.B, GBM.mito.Pro.Sig.B)
grid.newpage()
venn.plot.mito.B <- venn.diagram(x = list(LVQ.mito.Pro.Sig.B=LVQ.mito.Pro.Sig.B, SVM.mito.Pro.Sig.B=SVM.mito.Pro.Sig.B, GBM.mito.Pro.Sig.B=GBM.mito.Pro.Sig.B),
                                  filename=NULL, 
                                  fill = c("red", "blue", "green"),
                                  alpha = 0.50,
                                  col = "transparent")
grid.draw(venn.plot.mito.B)
venn.intersect.mito.B <- calculate.overlap(venn.data.mito.B)
print(venn.intersect.mito.B$a5)
print(venn.intersect.mito.B$a2)

#LVQ C_
test.mito.lvq.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.lvq.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = mitotane.pos.data.C, method = "lvq",
                               trControl = test.mito.lvq.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.lvq.C.pro <- varImp(model.mito.lvq.C.pro, scale = FALSE)
print(importance.mito.lvq.C.pro)
plot(importance.mito.lvq.C.pro)

#SVM C_
test.mito.svm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.svm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = mitotane.pos.data.C, method = "svmRadial",
                               trControl = test.mito.svm.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.svm.C.pro <- varImp(model.mito.svm.C.pro, scale = FALSE)
print(importance.mito.svm.C.pro)
plot(importance.mito.svm.C.pro)
#GBM C_
test.mito.gbm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.mito.gbm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                               data = mitotane.pos.data.C, method = "gbm",
                               trControl = test.mito.gbm.C.pro,
                               preProcess = c("center", "scale"),
                               tuneLength = 10,
                               na.action = na.pass)
importance.mito.gbm.C.pro <- varImp(model.mito.gbm.C.pro, scale = FALSE)
print(importance.mito.gbm.C.pro)
plot(importance.mito.gbm.C.pro)

#Common predictors for A_ progress..yes.no.
LVQ.mito.Pro.C <- importance.mito.lvq.C.pro$importance
SVM.mito.Pro.C <- importance.mito.svm.C.pro$importance
GBM.mito.Pro.C <- importance.mito.gbm.C.pro$importance
all(LVQ.mito.Pro.C$X0 == LVQ.mito.Pro.C$X1)
all(SVM.mito.Pro.C$X0 == SVM.mito.Pro.C$X1)
LVQ.mito.Pro.C$X1 <- NULL
SVM.mito.Pro.C$X1 <- NULL
colnames(LVQ.mito.Pro.C)[1] <- "Overall"
colnames(SVM.mito.Pro.C)[1] <- "Overall"
LVQ.mito.Pro.Sig.C <- setDT(LVQ.mito.Pro.C, keep.rownames = TRUE)[]
SVM.mito.Pro.Sig.C <- setDT(SVM.mito.Pro.C, keep.rownames = TRUE)[]
GBM.mito.Pro.Sig.C <- setDT(GBM.mito.Pro.C, keep.rownames = TRUE)[]
LVQ.mito.Pro.Sig.C <- LVQ.mito.Pro.Sig.C[order(-LVQ.mito.Pro.Sig.C$Overall),]
SVM.mito.Pro.Sig.C <- SVM.mito.Pro.Sig.C[order(-SVM.mito.Pro.Sig.C$Overall),]
GBM.mito.Pro.Sig.C <- GBM.mito.Pro.Sig.C[order(-GBM.mito.Pro.Sig.C$Overall),]
LVQ.mito.Pro.Sig.C <- head(LVQ.mito.Pro.Sig.C, 20)
SVM.mito.Pro.Sig.C <- head(SVM.mito.Pro.Sig.C, 20)
#Adjust GBM selection based on printed variables
GBM.mito.Pro.Sig.C <- head(GBM.mito.Pro.Sig.C, 3)
LVQ.mito.Pro.Sig.C$Overall <- NULL
SVM.mito.Pro.Sig.C$Overall <- NULL
GBM.mito.Pro.Sig.C$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.mito.Pro.Sig.C <- unlist(LVQ.mito.Pro.Sig.C)
SVM.mito.Pro.Sig.C <- unlist(SVM.mito.Pro.Sig.C)
GBM.mito.Pro.Sig.C <- unlist(GBM.mito.Pro.Sig.C)
library(VennDiagram)
venn.data.mito.C <- list(LVQ.mito.Pro.Sig.C, SVM.mito.Pro.Sig.C, GBM.mito.Pro.Sig.C)
grid.newpage()
venn.plot.mito.C <- venn.diagram(x = list(LVQ.mito.Pro.Sig.C=LVQ.mito.Pro.Sig.C, SVM.mito.Pro.Sig.C=SVM.mito.Pro.Sig.C, GBM.mito.Pro.Sig.C=GBM.mito.Pro.Sig.C),
                                  filename=NULL, 
                                  fill = c("red", "blue", "green"),
                                  alpha = 0.50,
                                  col = "transparent")
grid.draw(venn.plot.mito.C)
venn.intersect.mito.C <- calculate.overlap(venn.data.mito.C)
print(venn.intersect.mito.C$a5)
print(venn.intersect.mito.C$a2)


#Mitotane negative treatment 
#LVQ 
test.Nomitotane.lvq.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomitotane.lvq.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                                data = mitotane.neg.data.whole, method = "lvq",
                                trControl = test.Nomitotane.lvq.pro,
                                preProcess = c("center", "scale"),
                                tuneLength = 10,
                                tuneGrid = data.frame(size = 3, k = 1:2),
                                na.action = na.pass)
importance.Nomitotane.lvq.pro <- varImp(model.Nomitotane.lvq.pro, scale = FALSE)
print(importance.Nomitotane.lvq.pro)
plot(importance.Nomitotane.lvq.pro)

#SVM 
test.Nomitotane.svm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomitotane.svm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                                data = mitotane.neg.data.whole, method = "svmRadial",
                                trControl = test.Nomitotane.svm.pro,
                                preProcess = c("center", "scale"),
                                tuneLength = 10,
                                na.action = na.pass)
importance.Nomitotane.svm.pro <- varImp(model.Nomitotane.svm.pro, scale = FALSE)
print(importance.Nomitotane.svm.pro)
plot(importance.Nomitotane.svm.pro)


#GBM 
#test.Nomitotane.gbm.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.Nomitotane.gbm.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
#                                data = mitotane.neg.data.whole, method = "gbm",
#                                trControl = test.Nomitotane.gbm.pro,
#                                preProcess = c("center", "scale"),
#                                tuneLength = 10,
#                                na.action = na.pass)
#importance.Nomitotane.gbm.pro <- varImp(model.Nomitotane.gbm.pro, scale = FALSE)
#print(importance.Nomitotane.gbm.pro)
#plot(importance.Nomitotane.gbm.pro)

#Common predictors for  progress..yes.no.
LVQ.Nomitotane.Pro <- importance.Nomitotane.lvq.pro$importance
SVM.Nomitotane.Pro <- importance.Nomitotane.svm.pro$importance
#GBM.Nomitotane.Pro <- importance.Nomitotane.gbm.pro$importance
all(LVQ.Nomitotane.Pro$X0 == LVQ.Nomitotane.Pro$X1)
all(SVM.Nomitotane.Pro$X0 == SVM.Nomitotane.Pro$X1)
LVQ.Nomitotane.Pro$X1 <- NULL
SVM.Nomitotane.Pro$X1 <- NULL
colnames(LVQ.Nomitotane.Pro)[1] <- "Overall"
colnames(SVM.Nomitotane.Pro)[1] <- "Overall"
LVQ.Nomitotane.Pro <- setDT(LVQ.Nomitotane.Pro, keep.rownames = TRUE)[]
SVM.Nomitotane.Pro <- setDT(SVM.Nomitotane.Pro, keep.rownames = TRUE)[]
#GBM.Nomitotane.Pro <- setDT(GBM.Nomitotane.Pro, keep.rownames = TRUE)[]
LVQ.Nomitotane.Pro <- LVQ.Nomitotane.Pro[order(-LVQ.Nomitotane.Pro$Overall),]
SVM.Nomitotane.Pro <- SVM.Nomitotane.Pro[order(-SVM.Nomitotane.Pro$Overall),]
#GBM.Nomitotane.Pro <- GBM.Nomitotane.Pro[order(-GBM.Nomitotane.Pro$Overall),]
LVQ.Nomitotane.Pro.Sig <- head(LVQ.Nomitotane.Pro, 20)
SVM.Nomitotane.Pro.Sig <- head(SVM.Nomitotane.Pro, 20)
#Adjust GBM selection based on printed variables
#GBM.Nomitotane.Pro.Sig <- head(GBM.Nomitotane.Pro, 3)
LVQ.Nomitotane.Pro.Sig$Overall <- NULL
SVM.Nomitotane.Pro.Sig$Overall <- NULL
#GBM.Nomitotane.Pro.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.Nomitotane.Pro.Sig <- unlist(LVQ.Nomitotane.Pro.Sig)
SVM.Nomitotane.Pro.Sig <- unlist(SVM.Nomitotane.Pro.Sig)
#GBM.Nomitotane.Pro.Sig <- unlist(GBM.Nomitotane.Pro.Sig)
library(VennDiagram)
venn.data.Nomitotane.pos <- list(LVQ.Nomitotane.Pro.Sig, SVM.Nomitotane.Pro.Sig)
grid.newpage()
venn.plot.Nomitotane.pos <- venn.diagram(x = list(LVQ.Nomitotane.Pro.Sig=LVQ.Nomitotane.Pro.Sig, SVM.Nomitotane.Pro.Sig=SVM.Nomitotane.Pro.Sig),
                                       filename=NULL, 
                                       fill = c("red", "blue"),
                                       alpha = 0.50,
                                       col = "transparent")
grid.draw(venn.plot.Nomitotane.pos)
venn.intersect.Nomitotane.pos <- calculate.overlap(venn.data.Nomitotane.pos)


#Seperation for Nomitotane by A_, B_ and C_ using progress..yes.no.
#LVQ A_
test.Nomito.lvq.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomito.lvq.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF,
                              data = mitotane.neg.data.A, method = "lvq",
                              trControl = test.Nomito.lvq.A.pro,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.Nomito.lvq.A.pro <- varImp(model.Nomito.lvq.A.pro, scale = FALSE)
print(importance.Nomito.lvq.A.pro)
plot(importance.Nomito.lvq.A.pro)

#SVM A_
test.Nomito.svm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomito.svm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
                              data = mitotane.neg.data.A, method = "svmRadial",
                              trControl = test.Nomito.svm.A.pro,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.Nomito.svm.A.pro <- varImp(model.Nomito.svm.A.pro, scale = FALSE)
print(importance.Nomito.svm.A.pro)
plot(importance.Nomito.svm.A.pro)

#GBM A_
#test.Nomito.gbm.A.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.Nomito.gbm.A.pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF, 
#                              data = mitotane.neg.data.A, method = "gbm",
#                              trControl = test.Nomito.gbm.A.pro,
#                              preProcess = c("center", "scale"),
#                              tuneLength = 10,
#                              na.action = na.pass)
#importance.Nomito.gbm.A.pro <- varImp(model.Nomito.gbm.A.pro, scale = FALSE)
#print(importance.Nomito.gbm.A.pro)
#plot(importance.Nomito.gbm.A.pro)

#Common predictors for A_ progress..yes.no.
LVQ.Nomito.Pro.A <- importance.Nomito.lvq.A.pro$importance
SVM.Nomito.Pro.A <- importance.Nomito.svm.A.pro$importance
#GBM.Nomito.Pro.A <- importance.Nomito.gbm.A.pro$importance
all(LVQ.Nomito.Pro.A$X0 == LVQ.Nomito.Pro.A$X1)
all(SVM.Nomito.Pro.A$X0 == SVM.Nomito.Pro.A$X1)
LVQ.Nomito.Pro.A$X1 <- NULL
SVM.Nomito.Pro.A$X1 <- NULL
colnames(LVQ.Nomito.Pro.A)[1] <- "Overall"
colnames(SVM.Nomito.Pro.A)[1] <- "Overall"
LVQ.Nomito.Pro.Sig.A <- setDT(LVQ.Nomito.Pro.A, keep.rownames = TRUE)[]
SVM.Nomito.Pro.Sig.A <- setDT(SVM.Nomito.Pro.A, keep.rownames = TRUE)[]
#GBM.Nomito.Pro.Sig.A <- setDT(GBM.Nomito.Pro.A, keep.rownames = TRUE)[]
LVQ.Nomito.Pro.Sig.A <- LVQ.Nomito.Pro.Sig.A[order(-LVQ.Nomito.Pro.Sig.A$Overall),]
SVM.Nomito.Pro.Sig.A <- SVM.Nomito.Pro.Sig.A[order(-SVM.Nomito.Pro.Sig.A$Overall),]
#GBM.Nomito.Pro.Sig.A <- GBM.Nomito.Pro.Sig.A[order(-GBM.Nomito.Pro.Sig.A$Overall),]
LVQ.Nomito.Pro.Sig.A <- head(LVQ.Nomito.Pro.Sig.A, 20)
SVM.Nomito.Pro.Sig.A <- head(SVM.Nomito.Pro.Sig.A, 20)
#Adjust GBM selection based on printed variables
#GBM.Nomito.Pro.Sig.A <- head(GBM.Nomito.Pro.Sig.A, 3)
LVQ.Nomito.Pro.Sig.A$Overall <- NULL
SVM.Nomito.Pro.Sig.A$Overall <- NULL
#GBM.Nomito.Pro.Sig.A$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.Nomito.Pro.Sig.A <- unlist(LVQ.Nomito.Pro.Sig.A)
SVM.Nomito.Pro.Sig.A <- unlist(SVM.Nomito.Pro.Sig.A)
#GBM.Nomito.Pro.Sig.A <- unlist(GBM.Nomito.Pro.Sig.A)
library(VennDiagram)
venn.data.Nomito.A <- list(LVQ.Nomito.Pro.Sig.A, SVM.Nomito.Pro.Sig.A)
grid.newpage()
venn.plot.Nomito.A <- venn.diagram(x = list(LVQ.Nomito.Pro.Sig.A=LVQ.Nomito.Pro.Sig.A, SVM.Nomito.Pro.Sig.A=SVM.Nomito.Pro.Sig.A),
                                 filename=NULL, 
                                 fill = c("red", "blue"),
                                 alpha = 0.50,
                                 col = "transparent")
grid.draw(venn.plot.Nomito.A)
venn.intersect.Nomito.A <- calculate.overlap(venn.data.Nomito.A)


#LVQ B_
test.Nomito.lvq.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomito.lvq.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                              data = mitotane.neg.data.B, method = "lvq",
                              trControl = test.Nomito.lvq.B.pro,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.Nomito.lvq.B.pro <- varImp(model.Nomito.lvq.B.pro, scale = FALSE)
print(importance.Nomito.lvq.B.pro)
plot(importance.Nomito.lvq.B.pro)

#SVM B_
test.Nomito.svm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomito.svm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
                              data = mitotane.neg.data.B, method = "svmRadial",
                              trControl = test.Nomito.svm.B.pro,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.Nomito.svm.B.pro <- varImp(model.Nomito.svm.B.pro, scale = FALSE)
print(importance.Nomito.svm.B.pro)
plot(importance.Nomito.svm.B.pro)

#GBM B_
#test.Nomito.gbm.B.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.Nomito.gbm.B.pro <- train(progress..yes.no. ~ B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1, 
#                              data = mitotane.neg.data.B, method = "gbm",
#                              trControl = test.Nomito.gbm.B.pro,
#                              preProcess = c("center", "scale"),
#                              tuneLength = 10,
#                              na.action = na.pass)
#importance.Nomito.gbm.B.pro <- varImp(model.Nomito.gbm.B.pro, scale = FALSE)
#print(importance.Nomito.gbm.B.pro)
#plot(importance.Nomito.gbm.B.pro)

#Common predictors for A_ progress..yes.no.
LVQ.Nomito.Pro.B <- importance.Nomito.lvq.B.pro$importance
SVM.Nomito.Pro.B <- importance.Nomito.svm.B.pro$importance
#GBM.Nomito.Pro.B <- importance.Nomito.gbm.B.pro$importance
all(LVQ.Nomito.Pro.B$X0 == LVQ.Nomito.Pro.B$X1)
all(SVM.Nomito.Pro.B$X0 == SVM.Nomito.Pro.B$X1)
LVQ.Nomito.Pro.B$X1 <- NULL
SVM.Nomito.Pro.B$X1 <- NULL
colnames(LVQ.Nomito.Pro.B)[1] <- "Overall"
colnames(SVM.Nomito.Pro.B)[1] <- "Overall"
LVQ.Nomito.Pro.Sig.B <- setDT(LVQ.Nomito.Pro.B, keep.rownames = TRUE)[]
SVM.Nomito.Pro.Sig.B <- setDT(SVM.Nomito.Pro.B, keep.rownames = TRUE)[]
#GBM.Nomito.Pro.Sig.B <- setDT(GBM.Nomito.Pro.B, keep.rownames = TRUE)[]
LVQ.Nomito.Pro.Sig.B <- LVQ.Nomito.Pro.Sig.B[order(-LVQ.Nomito.Pro.Sig.B$Overall),]
SVM.Nomito.Pro.Sig.B <- SVM.Nomito.Pro.Sig.B[order(-SVM.Nomito.Pro.Sig.B$Overall),]
#GBM.Nomito.Pro.Sig.B <- GBM.Nomito.Pro.Sig.B[order(-GBM.Nomito.Pro.Sig.B$Overall),]
LVQ.Nomito.Pro.Sig.B <- head(LVQ.Nomito.Pro.Sig.B, 20)
SVM.Nomito.Pro.Sig.B <- head(SVM.Nomito.Pro.Sig.B, 20)
#Adjust GBM selection based on printed variables
#GBM.Nomito.Pro.Sig.B <- head(GBM.Nomito.Pro.Sig.B, 3)
LVQ.Nomito.Pro.Sig.B$Overall <- NULL
SVM.Nomito.Pro.Sig.B$Overall <- NULL
#GBM.Nomito.Pro.Sig.B$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.Nomito.Pro.Sig.B <- unlist(LVQ.Nomito.Pro.Sig.B)
SVM.Nomito.Pro.Sig.B <- unlist(SVM.Nomito.Pro.Sig.B)
#GBM.Nomito.Pro.Sig.B <- unlist(GBM.Nomito.Pro.Sig.B)
library(VennDiagram)
venn.data.Nomito.B <- list(LVQ.Nomito.Pro.Sig.B, SVM.Nomito.Pro.Sig.B)
grid.newpage()
venn.plot.Nomito.B <- venn.diagram(x = list(LVQ.Nomito.Pro.Sig.B=LVQ.Nomito.Pro.Sig.B, SVM.Nomito.Pro.Sig.B=SVM.Nomito.Pro.Sig.B),
                                 filename=NULL, 
                                 fill = c("red", "blue"),
                                 alpha = 0.50,
                                 col = "transparent")
grid.draw(venn.plot.Nomito.B)
venn.intersect.Nomito.B <- calculate.overlap(venn.data.Nomito.B)


#LVQ C_
test.Nomito.lvq.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomito.lvq.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                              data = mitotane.neg.data.C, method = "lvq",
                              trControl = test.Nomito.lvq.C.pro,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.Nomito.lvq.C.pro <- varImp(model.Nomito.lvq.C.pro, scale = FALSE)
print(importance.Nomito.lvq.C.pro)
plot(importance.Nomito.lvq.C.pro)

#SVM C_
test.Nomito.svm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.Nomito.svm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                              data = mitotane.neg.data.C, method = "svmRadial",
                              trControl = test.Nomito.svm.C.pro,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.Nomito.svm.C.pro <- varImp(model.Nomito.svm.C.pro, scale = FALSE)
print(importance.Nomito.svm.C.pro)
plot(importance.Nomito.svm.C.pro)

#GBM C_
#test.Nomito.gbm.C.pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.Nomito.gbm.C.pro <- train(progress..yes.no. ~ C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
#                              data = mitotane.neg.data.C, method = "gbm",
#                              trControl = test.Nomito.gbm.C.pro,
#                              preProcess = c("center", "scale"),
#                              tuneLength = 10,
#                              na.action = na.pass)
#importance.Nomito.gbm.C.pro <- varImp(model.Nomito.gbm.C.pro, scale = FALSE)
#print(importance.Nomito.gbm.C.pro)
#plot(importance.Nomito.gbm.C.pro)

#Common predictors for A_ progress..yes.no.
LVQ.Nomito.Pro.C <- importance.Nomito.lvq.C.pro$importance
SVM.Nomito.Pro.C <- importance.Nomito.svm.C.pro$importance
#GBM.Nomito.Pro.C <- importance.Nomito.gbm.C.pro$importance
all(LVQ.Nomito.Pro.C$X0 == LVQ.Nomito.Pro.C$X1)
all(SVM.Nomito.Pro.C$X0 == SVM.Nomito.Pro.C$X1)
LVQ.Nomito.Pro.C$X1 <- NULL
SVM.Nomito.Pro.C$X1 <- NULL
colnames(LVQ.Nomito.Pro.C)[1] <- "Overall"
colnames(SVM.Nomito.Pro.C)[1] <- "Overall"
LVQ.Nomito.Pro.Sig.C <- setDT(LVQ.Nomito.Pro.C, keep.rownames = TRUE)[]
SVM.Nomito.Pro.Sig.C <- setDT(SVM.Nomito.Pro.C, keep.rownames = TRUE)[]
#GBM.Nomito.Pro.Sig.C <- setDT(GBM.Nomito.Pro.C, keep.rownames = TRUE)[]
LVQ.Nomito.Pro.Sig.C <- LVQ.Nomito.Pro.Sig.C[order(-LVQ.Nomito.Pro.Sig.C$Overall),]
SVM.Nomito.Pro.Sig.C <- SVM.Nomito.Pro.Sig.C[order(-SVM.Nomito.Pro.Sig.C$Overall),]
#GBM.Nomito.Pro.Sig.C <- GBM.Nomito.Pro.Sig.C[order(-GBM.Nomito.Pro.Sig.C$Overall),]
LVQ.Nomito.Pro.Sig.C <- head(LVQ.Nomito.Pro.Sig.C, 20)
SVM.Nomito.Pro.Sig.C <- head(SVM.Nomito.Pro.Sig.C, 20)
#Adjust GBM selection based on printed variables
#GBM.Nomito.Pro.Sig.C <- head(GBM.Nomito.Pro.Sig.C, 3)
LVQ.Nomito.Pro.Sig.C$Overall <- NULL
SVM.Nomito.Pro.Sig.C$Overall <- NULL
#GBM.Nomito.Pro.Sig.C$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.Nomito.Pro.Sig.C <- unlist(LVQ.Nomito.Pro.Sig.C)
SVM.Nomito.Pro.Sig.C <- unlist(SVM.Nomito.Pro.Sig.C)
#GBM.Nomito.Pro.Sig.C <- unlist(GBM.Nomito.Pro.Sig.C)
library(VennDiagram)
venn.data.Nomito.C <- list(LVQ.Nomito.Pro.Sig.C, SVM.Nomito.Pro.Sig.C)
grid.newpage()
venn.plot.Nomito.C <- venn.diagram(x = list(LVQ.Nomito.Pro.Sig.C=LVQ.Nomito.Pro.Sig.C, SVM.Nomito.Pro.Sig.C=SVM.Nomito.Pro.Sig.C),
                                 filename=NULL, 
                                 fill = c("red", "blue"),
                                 alpha = 0.50,
                                 col = "transparent")
grid.draw(venn.plot.Nomito.C)
venn.intersect.Nomito.C <- calculate.overlap(venn.data.Nomito.C)
print(venn.intersect.Nomito.C$a5)
print(venn.intersect.Nomito.C$a2)