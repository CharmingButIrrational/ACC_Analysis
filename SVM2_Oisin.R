#Load data
mydata <- read.csv("C:/Users/User/Desktop/Analysis Data/SepMutData_Oisin")

#Splitting data for the SVM
mydata$X <- NULL #column numbers appeared as a seperate column

smp_size <- floor(0.8 * nrow(mydata))
set.seed(120)
train_ind <- sample(seq_len(nrow(mydata)), size = smp_size)

alt.data <- mydata[,c(34:180)]
alt.data$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- mydata$progress..yes.no.
PFS..months. <- mydata$PFS..months.
alt.data2 <- alt.data
alt.data <- cbind(progress..yes.no., alt.data)
alt.data2 <- cbind(PFS..months., alt.data2)

#Convert target variable to factor
alt.data[["progress..yes.no."]] = factor(alt.data[["progress..yes.no."]])

train.data <- alt.data[train_ind,]
test.data <- alt.data[-train_ind,]

library(caret)
library(kernlab)
library(e1071)
library(gbm)

#SVM using progress..yes.no.

trctrl2 <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

set.seed(52)

svm_Linear2 <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, data = train.data, method = "svmLinear",
                     trControl=trctrl2,
                     preProcess = c("center", "scale"),
                     tuneLength = 10,
                     na.action = na.pass)

svm_Linear2
#Predict classes with the test data
test_pred2 <- predict(svm_Linear2, newdata = test.data)
test_pred2
#Accuracy of ~ 41%
confusionMatrix(test_pred2, test.data$progress..yes.no.)



#Try using the stepAIC to simply the predictors
library(MASS)
glm.gamma.dist <- glm(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, family = "Gamma"(link=log), data = mydata)
step.glm.gamma.both <- stepAIC(glm.gamma.dist, direction = "both", trace = FALSE)
summary(step.glm.gamma.both)

glm.gamma.dist2 <- glm(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, family = "binomial", data = mydata)
step.glm.gamma.both2 <- stepAIC(glm.gamma.dist2, direction = "both", trace = FALSE)
summary(step.glm.gamma.both2)

svm_Linear_stepAIC <- train(progress..yes.no. ~ A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + 
                       A_MAP3K1 + A_ATRX + A_SMARCB1 + A_MYC + A_DAXX + A_CSFR1 + 
                       A_CREBBP + A_FBXO11 + A_ZNFR3 + A_ARID1A + A_BCL6 + A_BRCA1 + 
                       A_TP53 + A_MSH2 + A_RB1 + A_CTNNB1 + A_KDR + A_KMT2D + A_NF1 + 
                       A_CIC + A_JAK2 + A_ESR1 + A_APC + A_GNAS + A_MEN1 + A_DNMT3A + 
                       A_ZNRF3 + A_TSC1 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + 
                       A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_CIC + 
                       B_GNA11 + B_IL7R + B_GNAS + B_MAP3K1 + B_DDB2 + B_MDM2 + 
                       B_MEN1 + B_none + B_PIK3R1 + B_NOTCH1 + B_TERT + B_CDK4 + 
                       B_VHL + B_DDR2 + B_STK11 + C_SLC7A8 + C_PHF6 + C_FGFR2 + 
                       C_SPOP + C_none + C_ECT2L + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + 
                       C_CHEK2 + C_CSF1R + C_TNFRSF14 + C_VHL, data = train.data, method = "svmLinear",
                     trControl=trctrl2,
                     preProcess = c("center", "scale"),
                     tuneLength = 10,
                     na.action = na.pass)
svm_Linear_stepAIC

svm_Linear_stepAIC_predict <- predict(svm_Linear_stepAIC, newdata = test.data)
svm_Linear_stepAIC_predict

#Accuracy of ~ 50%
confusionMatrix(svm_Linear_stepAIC_predict, test.data$progress..yes.no.)

#Feature selection 1
set.seed(106)

#Error must change data to numeric, will fix 
alt.data.cor <- alt.data[,-1]
correlationMatrix <- cor(progress..yes.no., alt.data.cor)
print(correlationMatrix)
#May change cutoff if results are poor
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
print(highlyCorrelated)

#Feature selection 2
#Learning Vector Quantization
control.test2A <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.test2A <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                    data = alt.data, method = "lvq",
                    trControl=control.test2A,
                    preProcess = c("center", "scale"),
                    tuneLength = 10,
                    na.action = na.pass)
importanceA <- varImp(model.test2A, scale = FALSE)
print(importanceA)
plot(importanceA)
#Support Vector Machine
control.test2B <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.test2B <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                      data = alt.data, method = "svmRadial",
                      trControl=control.test2B,
                      preProcess = c("center", "scale"),
                      tuneLength = 10,
                      na.action = na.pass)
importanceB <- varImp(model.test2B, scale = FALSE)
print(importanceB)
plot(importanceB)
#Gradient Boosted Machine
control.test2C <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.test2C <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                      data = alt.data, method = "gbm",
                      trControl=control.test2C,
                      preProcess = c("center", "scale"),
                      tuneLength = 10,
                      na.action = na.pass)
importanceC <- varImp(model.test2C, scale = FALSE)
print(importanceC)
plot(importanceC)

resultsABC <- resamples(list(LVQ = model.test2A, GBM = model.test2B, SVM = model.test2C))

summary(resultsABC)

bwplot(resultsABC)

dotplot(resultsABC)

#Learning Vector Quantization
#Wrong model type for regression
#control.test3D <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#set.seed(86)
#model.test3D <- train(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
#                     data = alt.data2, method = "lvq",
#                     trControl=control.test3D,
#                     preProcess = c("center", "scale"),
#                     tuneLength = 10,
#                     na.action = na.pass)
#importance3D <- varImp(model.test3D, scale = FALSE)
#print(importance3D)
#plot(importance3D)

#Support Vector Machine
control.test3E <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.test3E <- train(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                     data = alt.data2, method = "svmRadial",
                     trControl=control.test3E,
                     preProcess = c("center", "scale"),
                     tuneLength = 10,
                     na.action = na.pass)
importance3E <- varImp(model.test3E, scale = FALSE)
print(importance3E)
plot(importance3E)
#Gradient Boosted Machine
control.test3F <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.test3F <- train(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                     data = alt.data2, method = "gbm",
                     trControl=control.test3F,
                     preProcess = c("center", "scale"),
                     tuneLength = 10,
                     na.action = na.pass)
importance3F <- varImp(model.test3F, scale = FALSE)
print(importance3F)
plot(importance3F)

resultsDEF <- resamples(list(GBM = model.test3E, SVM = model.test3F))

summary(resultsDEF)

bwplot(resultsDEF)

dotplot(resultsDEF)

#Comparison of different features selected
LVQ.Pro <- importanceA$importance
SVM.Pro <- importanceB$importance
GBM.Pro <- importanceC$importance
#SVM.Pro and LVQ.Pro have two columns X0 and X1 with the same values
all(LVQ.Pro$X0 == LVQ.Pro$X1)
all(SVM.Pro$X0 == SVM.Pro$X1)
#Remove and rename rows for comparison
LVQ.Pro$X1 <- NULL
SVM.Pro$X1 <- NULL
colnames(LVQ.Pro)[1] <- "Overall"
colnames(SVM.Pro)[1] <- "Overall"

SVM.PFS <- importance3E$importance
GBM.PFS <- importance3F$importance

#Select values from a certain level (currently based on 20 most significant factors)
LVQ.Pro.Sig <- subset(LVQ.Pro, Overall >= '0.5357', select = c("Overall"))
SVM.Pro.Sig <- subset(SVM.Pro, Overall >= '0.5357', select = c("Overall"))
GBM.Pro.Sig <- subset(GBM.Pro, Overall >= '0.0817', select = c("Overall"))
SVM.PFS.Sig <- subset(SVM.PFS, Overall >= '0.02062', select = c("Overall"))
GBM.PFS.Sig <- subset(GBM.PFS, Overall >= '0.02062', select = c("Overall"))

library(data.table)
#Convert row names to a column
LVQ.Pro.Sig <- setDT(LVQ.Pro.Sig, keep.rownames = TRUE)[]
SVM.Pro.Sig <- setDT(SVM.Pro.Sig, keep.rownames = TRUE)[]
GBM.Pro.Sig <- setDT(GBM.Pro.Sig, keep.rownames = TRUE)[]
SVM.PFS.Sig <- setDT(SVM.PFS.Sig, keep.rownames = TRUE)[]
GBM.PFS.Sig <- setDT(GBM.PFS.Sig, keep.rownames = TRUE)[]
#Remove values to make the intersect function work
LVQ.Pro.Sig$Overall <- NULL
SVM.Pro.Sig$Overall <- NULL
GBM.Pro.Sig$Overall <- NULL
SVM.PFS.Sig$Overall <- NULL
GBM.PFS.Sig$Overall <- NULL

LVQ.Pro.Sig <- unlist(LVQ.Pro.Sig)
SVM.Pro.Sig <- unlist(SVM.Pro.Sig)
GBM.Pro.Sig <- unlist(GBM.Pro.Sig)
SVM.PFS.Sig <- unlist(SVM.PFS.Sig)
GBM.PFS.Sig <- unlist(GBM.PFS.Sig)

library(VennDiagram)

venn.data <- list(LVQ.Pro.Sig, SVM.Pro.Sig, GBM.Pro.Sig, SVM.PFS.Sig, GBM.PFS.Sig)

grid.newpage()

venn.plot <- venn.diagram(x = list(LVQ.Pro.Sig=LVQ.Pro.Sig, SVM.Pro.Sig=SVM.Pro.Sig, GBM.Pro.Sig=GBM.Pro.Sig, SVM.PFS.Sig=SVM.PFS.Sig, GBM.PFS.Sig=GBM.PFS.Sig),
                          filename=NULL, 
                          fill = c("red", "blue", "green", "yellow", "darkorchid1"),
                          alpha = 0.50,
                          col = "transparent")

grid.draw(venn.plot)
venn.intersect <- calculate.overlap(venn.data)
print(venn.intersect$a31)

venn.data.pro <- list(LVQ.Pro.Sig, SVM.Pro.Sig, GBM.Pro.Sig)


venn.plot.Pro <- venn.diagram(x = list(LVQ.Pro.Sig=LVQ.Pro.Sig, SVM.Pro.Sig=SVM.Pro.Sig, GBM.Pro.Sig=GBM.Pro.Sig),
                          filename=NULL, 
                          fill = c("red", "blue", "green"),
                          alpha = 0.50,
                          col = "transparent")
grid.draw(venn.plot.Pro)
venn.intersect.pro <- calculate.overlap(venn.data.pro)
print(venn.intersect.pro$a5)

#Feature selection 3
set.seed(57)
control.test4 <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
results.test4 <- rfe(alt.data[,2:146], alt.data[,1], sizes = c(2:146), rfeControl = control.test4)
print(results.test4)
predictors(results.test4)
plot(results.test4, type = c("g", "o"))

set.seed(436)
control.test5 <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
results.test5 <- rfe(alt.data2[,2:146], alt.data2[,1], sizes = c(2:146), rfeControl = control.test5)
print(results.test5)
predictors(results.test5)
plot(results.test5, type = c("g", "o"))




#################################################################################################################################

#Load data
mydata <- read.csv("C:/Users/User/Desktop/Analysis Data/SepMutData_Oisin")

#Splitting data for the SVM
mydata$X <- NULL #column numbers appeared as a seperate column
alt.data11 <- mydata[,c(34:180)]
alt.data11$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data11$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(mydata$progress..yes.no.)
alt.data11 <- cbind(progress..yes.no., alt.data11)

#Seperating data into A_, B_, C_
alt.data11.A <- alt.data11[,c(1:72)]

alt.data11.B <- alt.data11[,c(73:105)]
alt.data11.B <- cbind(progress..yes.no., alt.data11.B)

alt.data11.C <- alt.data11[,c(105:146)]
alt.data11.C <- cbind(progress..yes.no., alt.data11.C)

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
LVQ.Pro.A <- importance.chemo.lvq.A.pro$importance
SVM.Pro.A <- importance.chemo.svm.A.pro$importance
GBM.Pro.A <- importance.chemo.gbm.A.pro$importance
all(LVQ.Pro.A$X0 == LVQ.Pro.A$X1)
all(SVM.Pro.A$X0 == SVM.Pro.A$X1)
LVQ.Pro.A$X1 <- NULL
SVM.Pro.A$X1 <- NULL
colnames(LVQ.Pro.A)[1] <- "Overall"
colnames(SVM.Pro.A)[1] <- "Overall"
LVQ.Pro.Sig.A <- subset(LVQ.Pro.A, Overall >= '0.5179', select = c("Overall"))
SVM.Pro.Sig.A <- subset(SVM.Pro.A, Overall >= '0.5167', select = c("Overall"))
GBM.Pro.Sig.A <- subset(GBM.Pro.A, Overall >= '0.1337', select = c("Overall"))
library(data.table)
LVQ.Pro.Sig.A <- setDT(LVQ.Pro.Sig.A, keep.rownames = TRUE)[]
SVM.Pro.Sig.A <- setDT(SVM.Pro.Sig.A, keep.rownames = TRUE)[]
GBM.Pro.Sig.A <- setDT(GBM.Pro.Sig.A, keep.rownames = TRUE)[]
LVQ.Pro.Sig.A$Overall <- NULL
SVM.Pro.Sig.A$Overall <- NULL
GBM.Pro.Sig.A$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.Pro.Sig.A <- unlist(LVQ.Pro.Sig.A)
SVM.Pro.Sig.A <- unlist(SVM.Pro.Sig.A)
GBM.Pro.Sig.A <- unlist(GBM.Pro.Sig.A)
library(VennDiagram)
venn.data.A <- list(LVQ.Pro.Sig.A, SVM.Pro.Sig.A, GBM.Pro.Sig.A)
grid.newpage()
venn.plot.A <- venn.diagram(x = list(LVQ.Pro.Sig.A=LVQ.Pro.Sig.A, SVM.Pro.Sig.A=SVM.Pro.Sig.A, GBM.Pro.Sig.A=GBM.Pro.Sig.A),
                            filename=NULL, 
                            fill = c("red", "blue", "green"),
                            alpha = 0.50,
                            col = "transparent")
grid.draw(venn.plot.A)
venn.intersect.A <- calculate.overlap(venn.data.A)
print(venn.intersect.A$a5)

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
LVQ.Pro.Sig.B <- subset(LVQ.Pro.B, Overall >= '0.5060', select = c("Overall"))
SVM.Pro.Sig.B <- subset(SVM.Pro.B, Overall >= '0.5060', select = c("Overall"))
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
LVQ.Pro.Sig.C <- subset(LVQ.Pro.C, Overall >= '0.5119', select = c("Overall"))
SVM.Pro.Sig.C <- subset(SVM.Pro.C, Overall >= '0.5119', select = c("Overall"))
GBM.Pro.Sig.C <- subset(GBM.Pro.C, Overall >= '0.03851', select = c("Overall"))
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

