#Load data
mydata <- read.csv("C:/Users/User/Desktop/Analysis Data/SepMutData_Oisin")
mydata$X <- NULL #column numbers appeared as a seperate column

smp_size <- floor(0.8 * nrow(mydata))
set.seed(167)
train_ind <- sample(seq_len(nrow(mydata)), size = smp_size)

alt.data <- mydata[,c(34:180)]
alt.data$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- mydata$progress..yes.no.
PFS..months. <- mydata$PFS..months.
alt.data.PFS <- alt.data
alt.data.Pro <- cbind(progress..yes.no., alt.data)
alt.data.PFS <- cbind(PFS..months., alt.data.PFS)

alt.data.Pro[["progress..yes.no."]] = factor(alt.data.Pro[["progress..yes.no."]])

train.data.Pro <- alt.data.Pro[train_ind,]
test.data.Pro <- alt.data.Pro[-train_ind,]

train.data.PFS <- alt.data.PFS[train_ind,]
test.data.PFS <- alt.data.PFS[-train_ind,]

library(caret)
library(e1071)
library(gbm)
library(MASS)
library(randomForest)
library(ROCR)

#########################################################################################################

#GLM using PFS..months
glm.gamma.PFS <- glm(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + 
                        A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + 
                        A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + 
                        A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + 
                        A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + 
                        A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + 
                        A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + 
                        A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + 
                        A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + 
                        A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + 
                        A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + 
                        A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                        A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + 
                        B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + 
                        B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + 
                        B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + 
                        B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + 
                        B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + 
                        C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                        C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + 
                        C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + 
                        C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + 
                        C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + 
                        C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                        family = "Gamma"(link=log), data = train.data.PFS)
summary(glm.gamma.PFS)
predict.glm.PFS <- predict(glm.gamma.PFS, test.data.PFS)

#Step.AIC using PFS..months. (from MASS package)
step.glm.gamma.PFS.both <- stepAIC(glm.gamma.PFS, direction = "both", trace = FALSE)
summary(step.glm.gamma.PFS.both)
predict.step.glm.PFS <- predict(step.glm.gamma.PFS.both, test.data.PFS)

#SVM using PFS..months. (from caret package)
control.svm.PFS <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.svm.PFS <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + 
                         A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + 
                         A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + 
                         A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + 
                         A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + 
                         A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + 
                         A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + 
                         A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + 
                         A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                         A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + 
                         B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + 
                         B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + 
                         B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + 
                         B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + 
                         C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                         C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + 
                         C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + 
                         C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + 
                         C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                      data = train.data.PFS, method = "svmRadial",
                      trControl = control.svm.PFS,
                      preProcess = c("center", "scale"),
                      tuneLength = 10,
                      na.action = na.pass)
predict.svm.PFS <- predict(model.svm.PFS, test.data.PFS)

#GBM using PFS..months. (from gbm package)
control.gbm.PFS <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.gbm.PFS <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + 
                         A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + 
                         A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + 
                         A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + 
                         A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + 
                         A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + 
                         A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + 
                         A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + 
                         A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                         A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + 
                         B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + 
                         B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + 
                         B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + 
                         B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + 
                         C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                         C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + 
                         C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + 
                         C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + 
                         C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                       data = train.data.PFS, method = "svmRadial",
                       trControl = control.gbm.PFS,
                       preProcess = c("center", "scale"),
                       tuneLength = 10,
                       na.action = na.pass)
predict.gbm.PFS <- predict.gbm(model.gbm.PFS, test.data.PFS)

#RandomForest with PFS..months. (from randomForest package)
model.rf.PFS <- randomForest(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + 
                               A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + 
                               A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + 
                               A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + 
                               A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + 
                               A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + 
                               A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + 
                               A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + 
                               A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                               A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + 
                               B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + 
                               B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + 
                               B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + 
                               B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + 
                               C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                               C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + 
                               C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + 
                               C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + 
                               C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                             data = train.data.PFS, importance = TRUE)
predict.rf.PFS <- predict(model.rf.PFS, test.data.PFS)

#########################################################################################################

#GLM using progress..yes.no
glm.binomial.Pro <- glm(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + 
                        A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + 
                        A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + 
                        A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + 
                        A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + 
                        A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + 
                        A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + 
                        A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + 
                        A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + 
                        A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + 
                        A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + 
                        A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                        A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + 
                        B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + 
                        B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + 
                        B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + 
                        B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + 
                        B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + 
                        C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                        C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + 
                        C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + 
                        C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + 
                        C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + 
                        C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                        family = "binomial"(link=log), data = train.data.Pro)
summary(glm.binomial.Pro)
predict.glm.Pro <- predict(glm.binomial.Pro, test.data.Pro)
confusion.glm.Pro <- confusionMatrix(predict.glm.Pro, test.data.Pro$progress..yes.no.)
#Plotting ROC graph
ROC.glm.Pro <- performance(predict.glm.Pro, "tpr", "fpr")
plot(ROC.glm.Pro)
plot(ROC.glm.Pro, add = TRUE, col = "green")
legend("right", legend = c("rf"), col = c("green"), lty = 1:2, cex = 0.6)
#Calculating AUC
AUC.glm.Pro <- performance(predict.glm.Pro, "AUC")
AUC.glm.test.Pro <- as.numeric(AUC.glm.Pro@y.values)
AUC.glm.test.Pro

#Step.AIC using progress..yes.no.  (from MASS package)
step.glm.binomial.Pro.both <- stepAIC(glm.binomial.Pro, direction = "both", trace = FALSE)
summary(step.glm.binomial.Pro.both)
predict.step.glm.Pro <- predict(step.glm.binomial.Pro.both, test.data.Pro)
confusion.step.glm.Pro <- confusionMatrix(predict.step.glm.Pro, test.data.Pro$progress..yes.no.)
#Plotting ROC graph
ROC.AIC.Pro <- performance(predict.step.glm.Pro, "tpr", "fpr")
plot(ROC.AIC.Pro)
plot(ROC.AIC.Pro, add = TRUE, col = "green")
legend("right", legend = c("rf"), col = c("green"), lty = 1:2, cex = 0.6)
#Calculating AUC
AUC.AIC.Pro <- performance(predict.step.glm.Pro, "AUC")
AUC.AIC.test.Pro <- as.numeric(AUC.AIC.Pro@y.values)
AUC.AIC.test.Pro

#SVM using progress..yes.no.(from caret package)
control.svm.Pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.svm.Pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + 
                         A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + 
                         A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + 
                         A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + 
                         A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + 
                         A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + 
                         A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + 
                         A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + 
                         A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                         A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + 
                         B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + 
                         B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + 
                         B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + 
                         B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + 
                         C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                         C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + 
                         C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + 
                         C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + 
                         C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                       data = train.data.Pro, method = "svmRadial",
                       trControl = control.svm.Pro,
                       preProcess = c("center", "scale"),
                       tuneLength = 10,
                       na.action = na.pass)
predict.svm.Pro <- predict(model.svm.Pro, test.data.Pro)
confusion.svm.Pro <- confusionMatrix(predict.svm.Pro, test.data.Pro$progress..yes.no.)
#Plotting ROC graph
ROC.svm.Pro <- performance(predict.svm.Pro, "tpr", "fpr")
plot(ROC.svm.Pro)
plot(ROC.svm.Pro, add = TRUE, col = "green")
legend("right", legend = c("rf"), col = c("green"), lty = 1:2, cex = 0.6)
#Calculating AUC
AUC.svm.Pro <- performance(predict.svm.Pro, "AUC")
AUC.svm.test.Pro <- as.numeric(AUC.svm.Pro@y.values)
AUC.svm.test.Pro

#GBM using progress..yes.no. (from gbm package)
control.gbm.Pro <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
model.gbm.Pro <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + 
                         A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + 
                         A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + 
                         A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + 
                         A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + 
                         A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + 
                         A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + 
                         A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + 
                         A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                         A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + 
                         B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + 
                         B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + 
                         B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + 
                         B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + 
                         C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                         C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + 
                         C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + 
                         C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + 
                         C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                       data = train.data.Pro, method = "svmRadial",
                       trControl = control.gbm.Pro,
                       preProcess = c("center", "scale"),
                       tuneLength = 10,
                       na.action = na.pass)
predict.gbm.Pro <- predict.gbm(model.gbm.Pro, test.data.Pro)
confusion.gbm.Pro <- confusionMatrix(predict.gbm.Pro, test.data.Pro$progress..yes.no.)
#Plotting ROC graph
ROC.gbm.Pro <- performance(predict.gbm.Pro, "tpr", "fpr")
plot(ROC.gbm.Pro)
plot(ROC.gbm.Pro, add = TRUE, col = "green")
legend("right", legend = c("rf"), col = c("green"), lty = 1:2, cex = 0.6)
#Calculating AUC
AUC.gbm.Pro <- performance(predict.gbm.Pro, "AUC")
AUC.gbm.test.Pro <- as.numeric(AUC.gbm.Pro@y.values)
AUC.gbm.test.Pro

#RandomForest with progress..yes.no (from randomForest package)
model.rf.Pro <- randomForest(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + 
                              A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + 
                              A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + 
                              A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + 
                              A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + 
                              A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + 
                              A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + 
                              A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + 
                              A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                              A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + 
                              B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + 
                              B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + 
                              B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + 
                              B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + 
                              C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + 
                              C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + 
                              C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + 
                              C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + 
                              C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, 
                            data = train.data.Pro, importance = TRUE)
predict.rf.Pro <- predict(model.rf.Pro, test.data.Pro)
confusion.rf.Pro <- confusionMatrix(predict.rf.Pro, test.data.Pro$progress..yes.no.)
#Plotting ROC graph
ROC.rf.Pro <- performance(predict.rf.Pro, "tpr", "fpr")
plot(ROC.rf.Pro)
plot(ROC.rf.Pro, add = TRUE, col = "green")
legend("right", legend = c("rf"), col = c("green"), lty = 1:2, cex = 0.6)
#Calculating AUC
AUC.rf.Pro <- performance(predict.rf.Pro, "AUC")
AUC.rf.test.Pro <- as.numeric(AUC.rf.Pro@y.values)
AUC.rf.test.Pro