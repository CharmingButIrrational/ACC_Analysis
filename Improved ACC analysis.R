#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/SepMutData_Oisin")

mydata$X <- NULL #column numbers appeared as a seperate column
alt.data.whole <- mydata[,c(34:180)]
alt.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(mydata$progress..yes.no.)
alt.data.whole <- cbind(progress..yes.no., alt.data.whole)

#Split df by value in the Chemotherapy column
alt.data <- split(mydata, mydata$Chemotherapy)
alt.data.nochemo <- alt.data$`0`
alt.data.chemo <- alt.data$`1`

#Split nochemo data by mitotane (All forms of therapy)
no.chemo.mito.split <- split(alt.data.nochemo, alt.data.nochemo$Mitotane..pall.or.adj.)
#Data with mitotane only
mitotane.nochemo <- rbind(no.chemo.mito.split$adj, no.chemo.mito.split$`adj+pall`, no.chemo.mito.split$pall)
progress..yes.no. <- as.factor(mitotane.nochemo$progress..yes.no.)
mitotane.nochemo <- mitotane.nochemo[,c(34:180)]
mitotane.nochemo$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.nochemo$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
mitotane.nochemo <- cbind(progress..yes.no., mitotane.nochemo)
#Data which recieved no therapy
no.therapy <- no.chemo.mito.split$`0`
progress..yes.no. <- as.factor(no.therapy$progress..yes.no.)
no.therapy <- no.therapy[,c(34:180)]
no.therapy$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
no.therapy$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
no.therapy <- cbind(progress..yes.no., no.therapy)

#Split chemo data by mitotane (All forms of therapy)
chemo.mito.split <- split(alt.data.chemo, alt.data.chemo$Mitotane..pall.or.adj.)
#Data with both chemotherapy and mitotane therapy
mitotane.chemo <- rbind(chemo.mito.split$adj, chemo.mito.split$`adj+pall`, chemo.mito.split$pall, chemo.mito.split$`adj (with chemo)`)
progress..yes.no. <- as.factor(mitotane.chemo$progress..yes.no.)
mitotane.chemo <- mitotane.chemo[,c(34:180)]
mitotane.chemo$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.chemo$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
mitotane.chemo <- cbind(progress..yes.no., mitotane.chemo)
#Data with chemotherapy only
chemo.nomito <- chemo.mito.split$`0`
progress..yes.no. <- as.factor(chemo.nomito$progress..yes.no.)
chemo.nomito <- chemo.nomito[,c(34:180)]
chemo.nomito$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
chemo.nomito$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
chemo.nomito <- cbind(progress..yes.no., chemo.nomito)


#Chemotherapy only
#LVQ 

chemo.nomito.lvq <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
chemo.nomito.lvq <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                                data = chemo.nomito, method = "lvq",
                                trControl = chemo.nomito.lvq,
                                preProcess = c("center", "scale"),
                                tuneGrid = ,
                                tuneLength = 10,
                                na.action = na.pass)
importance.chemo.nomito.lvq <- varImp(chemo.nomito.lvq, scale = FALSE)
print(importance.chemo.nomito.lvq)
plot(importance.chemo.nomito.lvq)

#SVM 
set.seed(86)
chemo.nomito.svm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                          data = chemo.nomito, method = "svmLinear",
                          trControl = chemo.nomito.svm,
                          preProcess = c("center", "scale"),
                          tuneGrid = expand.grid(C = seq(0, 2, length = 20)),
                          tuneLength = 10,
                          na.action = na.pass)
importance.chemo.nomito.svm <- varImp(chemo.nomito.svm, scale = FALSE)
print(importance.chemo.nomito.svm)
plot(importance.chemo.nomito.svm)

#Trying a nonlinear kernal
set.seed(86)
chemo.nomito.svm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                          data = chemo.nomito, method = "svmRadial",
                          trControl = chemo.nomito.svm,
                          preProcess = c("center", "scale"),
                          tuneLength = 10,
                          na.action = na.pass)
importance.chemo.nomito.svm <- varImp(chemo.nomito.svm, scale = FALSE)
print(importance.chemo.nomito.svm)
plot(importance.chemo.nomito.svm)

#GBM 
chemo.nomito.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
chemo.nomito.gbm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                          data = chemo.nomito, method = "gbm",
                          trControl = chemo.nomito.gbm,
                          preProcess = c("center", "scale"),
                          tuneLength = 10,
                          na.action = na.pass)
importance.chemo.nomito.gbm <- varImp(chemo.nomito.gbm, scale = FALSE)
print(importance.chemo.nomito.gbm)
plot(importance.chemo.nomito.gbm)

#Common predictors for  progress..yes.no.
LVQ.chemo.nomito <- importance.chemo.nomito.lvq$importance
SVM.chemo.nomito <- importance.chemo.nomito.svm$importance
GBM.chemo.nomito <- importance.chemo.nomito.gbm$importance
all(LVQ.chemo.nomito$X0 == LVQ.chemo.nomito$X1)
all(SVM.chemo.nomito$X0 == SVM.chemo.nomito$X1)
LVQ.chemo.nomito$X1 <- NULL
SVM.chemo.nomito$X1 <- NULL
colnames(LVQ.chemo.nomito)[1] <- "Overall"
colnames(SVM.chemo.nomito)[1] <- "Overall"
LVQ.chemo.nomito <- setDT(LVQ.chemo.nomito, keep.rownames = TRUE)[]
SVM.chemo.nomito <- setDT(SVM.chemo.nomito, keep.rownames = TRUE)[]
GBM.chemo.nomito <- setDT(GBM.chemo.nomito, keep.rownames = TRUE)[]
LVQ.chemo.nomito <- LVQ.chemo.nomito[order(-LVQ.chemo.nomito$Overall),]
SVM.chemo.nomito <- SVM.chemo.nomito[order(-SVM.chemo.nomito$Overall),]
GBM.chemo.nomito <- GBM.chemo.nomito[order(-GBM.chemo.nomito$Overall),]
LVQ.chemo.nomito.Sig <- head(LVQ.chemo.nomito, 20)
SVM.chemo.nomito.Sig <- head(SVM.chemo.nomito, 20)
#Adjust GBM selection based on printed variables
GBM.chemo.nomito.Sig <- head(GBM.chemo.nomito, 9)
LVQ.chemo.nomito.Sig$Overall <- NULL
SVM.chemo.nomito.Sig$Overall <- NULL
GBM.chemo.nomito.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.chemo.nomito.Sig <- unlist(LVQ.chemo.nomito.Sig)
SVM.chemo.nomito.Sig <- unlist(SVM.chemo.nomito.Sig)
GBM.chemo.nomito.Sig <- unlist(GBM.chemo.nomito.Sig)
library(VennDiagram)
venn.data.chemo.nomito <- list(LVQ.chemo.nomito.Sig, SVM.chemo.nomito.Sig, GBM.chemo.nomito.Sig)
grid.newpage()
venn.plot.chemo.nomito <- venn.diagram(x = list(LVQ.chemo.nomito.Sig=LVQ.chemo.nomito.Sig, SVM.chemo.nomito.Sig=SVM.chemo.nomito.Sig, GBM.chemo.nomito.Sig=GBM.chemo.nomito.Sig),
                                       filename=NULL, 
                                       fill = c("red", "blue", "green"),
                                       alpha = 0.50,
                                       col = "transparent")
grid.draw(venn.plot.chemo.nomito)
venn.intersect.chemo.nomito <- calculate.overlap(venn.data.chemo.nomito)
print(venn.intersect.chemo.nomito$a5)
print(venn.intersect.chemo.nomito$a2)


#Mitotane only
#LVQ 
mitotane.nochemo.lvq <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
mitotane.nochemo.lvq <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                          data = mitotane.nochemo, method = "lvq",
                          trControl = mitotane.nochemo.lvq,
                          preProcess = c("center", "scale"),
                          tuneGrid = data.frame(size = 20, k = 1:100),
                          tuneLength = 10,
                          na.action = na.pass)
importance.mitotane.nochemo.lvq <- varImp(mitotane.nochemo.lvq, scale = FALSE)
print(importance.mitotane.nochemo.lvq)
plot(importance.mitotane.nochemo.lvq, 20)

#SVM 
mitotane.nochemo.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
mitotane.nochemo.svm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                          data = mitotane.nochemo, method = "svmLinear",
                          trControl = mitotane.nochemo.svm,
                          preProcess = c("center", "scale"),
                          tuneLength = 10,
                          na.action = na.pass)
importance.mitotane.nochemo.svm <- varImp(mitotane.nochemo.svm, scale = FALSE)
print(importance.mitotane.nochemo.svm)
plot(importance.mitotane.nochemo.svm, 20)

#GBM 
myGrid <- expand.grid(interaction.depth = seq(1, 7, by = 2),
                        n.trees = seq(100, 1000, by = 50),
                        n.minobsinnode = 10,
                        shrinkage = c(0.000001, 0.1))
mitotane.nochemo.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
mitotane.nochemo.gbm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                          data = mitotane.nochemo, method = "gbm",
                          trControl = mitotane.nochemo.gbm,
                          preProcess = c("center", "scale"),
                          tuneGrid = myGrid,
                          tuneLength = 10,
                          na.action = na.pass)
importance.mitotane.nochemo.gbm <- varImp(mitotane.nochemo.gbm, scale = FALSE)
print(importance.mitotane.nochemo.gbm)
plot(importance.mitotane.nochemo.gbm)

#Common predictors for  progress..yes.no.
LVQ.mitotane.nochemo <- importance.mitotane.nochemo.lvq$importance
SVM.mitotane.nochemo <- importance.mitotane.nochemo.svm$importance
#GBM.mitotane.nochemo <- importance.mitotane.nochemo.gbm$importance
all(LVQ.mitotane.nochemo$X0 == LVQ.mitotane.nochemo$X1)
all(SVM.mitotane.nochemo$X0 == SVM.mitotane.nochemo$X1)
LVQ.mitotane.nochemo$X1 <- NULL
SVM.mitotane.nochemo$X1 <- NULL
colnames(LVQ.mitotane.nochemo)[1] <- "Overall"
colnames(SVM.mitotane.nochemo)[1] <- "Overall"
LVQ.mitotane.nochemo <- setDT(LVQ.mitotane.nochemo, keep.rownames = TRUE)[]
SVM.mitotane.nochemo <- setDT(SVM.mitotane.nochemo, keep.rownames = TRUE)[]
#GBM.mitotane.nochemo <- setDT(GBM.mitotane.nochemo, keep.rownames = TRUE)[]
LVQ.mitotane.nochemo <- LVQ.mitotane.nochemo[order(-LVQ.mitotane.nochemo$Overall),]
SVM.mitotane.nochemo <- SVM.mitotane.nochemo[order(-SVM.mitotane.nochemo$Overall),]
#GBM.mitotane.nochemo <- GBM.mitotane.nochemo[order(-GBM.mitotane.nochemo$Overall),]
LVQ.mitotane.nochemo.Sig <- head(LVQ.mitotane.nochemo, 20)
SVM.mitotane.nochemo.Sig <- head(SVM.mitotane.nochemo, 20)
#Adjust GBM selection based on printed variables
#GBM.mitotane.nochemo.Sig <- head(GBM.mitotane.nochemo, 9)
LVQ.mitotane.nochemo.Sig$Overall <- NULL
SVM.mitotane.nochemo.Sig$Overall <- NULL
#GBM.mitotane.nochemo.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.mitotane.nochemo.Sig <- unlist(LVQ.mitotane.nochemo.Sig)
SVM.mitotane.nochemo.Sig <- unlist(SVM.mitotane.nochemo.Sig)
#GBM.mitotane.nochemo.Sig <- unlist(GBM.mitotane.nochemo.Sig)
library(VennDiagram)
venn.data.mitotane.nochemo <- list(LVQ.mitotane.nochemo.Sig, SVM.mitotane.nochemo.Sig)
grid.newpage()
venn.plot.mitotane.nochemo <- venn.diagram(x = list(LVQ.mitotane.nochemo.Sig=LVQ.mitotane.nochemo.Sig, SVM.mitotane.nochemo.Sig=SVM.mitotane.nochemo.Sig),
                                       filename=NULL, 
                                       fill = c("red", "blue"),
                                       alpha = 0.50,
                                       col = "transparent")
grid.draw(venn.plot.mitotane.nochemo)
venn.intersect.mitotane.nochemo <- calculate.overlap(venn.data.mitotane.nochemo)
print(venn.intersect.mitotane.nochemo$a5)
print(venn.intersect.mitotane.nochemo$a2)


#No therapy
#LVQ 
no.therapy.lvq <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
no.therapy.lvq <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                              data = no.therapy, method = "lvq",
                              trControl = no.therapy.lvq,
                              preProcess = c("center", "scale"),
                              tuneGrid = data.frame(size = 20, k = 1:100),
                              tuneLength = 10,
                              na.action = na.pass)
importance.no.therapy.lvq <- varImp(no.therapy.lvq, scale = FALSE)
print(importance.no.therapy.lvq)
plot(importance.no.therapy.lvq, 20)

#SVM 
no.therapy.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
no.therapy.svm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                              data = no.therapy, method = "svmLinear",
                              trControl = no.therapy.svm,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.no.therapy.svm <- varImp(no.therapy.svm, scale = FALSE)
print(importance.no.therapy.svm)
plot(importance.no.therapy.svm, 20)

#GBM 
no.therapy.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
no.therapy.gbm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                              data = no.therapy, method = "gbm",
                              trControl = no.therapy.gbm,
                              preProcess = c("center", "scale"),
                              tuneLength = 10,
                              na.action = na.pass)
importance.no.therapy.gbm <- varImp(no.therapy.gbm, scale = FALSE)
print(importance.no.therapy.gbm)
plot(importance.no.therapy.gbm)

#Common predictors for  progress..yes.no.
LVQ.no.therapy <- importance.no.therapy.lvq$importance
SVM.no.therapy <- importance.no.therapy.svm$importance
#GBM.no.therapy <- importance.no.therapy.gbm$importance
all(LVQ.no.therapy$X0 == LVQ.no.therapy$X1)
all(SVM.no.therapy$X0 == SVM.no.therapy$X1)
LVQ.no.therapy$X1 <- NULL
SVM.no.therapy$X1 <- NULL
colnames(LVQ.no.therapy)[1] <- "Overall"
colnames(SVM.no.therapy)[1] <- "Overall"
LVQ.no.therapy <- setDT(LVQ.no.therapy, keep.rownames = TRUE)[]
SVM.no.therapy <- setDT(SVM.no.therapy, keep.rownames = TRUE)[]
#GBM.no.therapy <- setDT(GBM.no.therapy, keep.rownames = TRUE)[]
LVQ.no.therapy <- LVQ.no.therapy[order(-LVQ.no.therapy$Overall),]
SVM.no.therapy <- SVM.no.therapy[order(-SVM.no.therapy$Overall),]
#GBM.no.therapy <- GBM.no.therapy[order(-GBM.no.therapy$Overall),]
LVQ.no.therapy.Sig <- head(LVQ.no.therapy, 20)
SVM.no.therapy.Sig <- head(SVM.no.therapy, 20)
#Adjust GBM selection based on printed variables
#GBM.no.therapy.Sig <- head(GBM.no.therapy, 9)
LVQ.no.therapy.Sig$Overall <- NULL
SVM.no.therapy.Sig$Overall <- NULL
#GBM.no.therapy.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.no.therapy.Sig <- unlist(LVQ.no.therapy.Sig)
SVM.no.therapy.Sig <- unlist(SVM.no.therapy.Sig)
#GBM.no.therapy.Sig <- unlist(GBM.no.therapy.Sig)
library(VennDiagram)
venn.data.no.therapy <- list(LVQ.no.therapy.Sig, SVM.no.therapy.Sig)
grid.newpage()
venn.plot.no.therapy <- venn.diagram(x = list(LVQ.no.therapy.Sig=LVQ.no.therapy.Sig, SVM.no.therapy.Sig=SVM.no.therapy.Sig),
                                           filename=NULL, 
                                           fill = c("red", "blue"),
                                           alpha = 0.50,
                                           col = "transparent")
grid.draw(venn.plot.no.therapy)
venn.intersect.no.therapy <- calculate.overlap(venn.data.no.therapy)
print(venn.intersect.no.therapy$a5)
print(venn.intersect.no.therapy$a2)


#Chemotherapy and Mitotane combined therapies
#LVQ 
mitotane.chemo.lvq <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
mitotane.chemo.lvq <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                        data = mitotane.chemo, method = "lvq",
                        trControl = mitotane.chemo.lvq,
                        preProcess = c("center", "scale"),
                        tuneGrid = data.frame(size = 20, k = 1:100),
                        tuneLength = 10,
                        na.action = na.pass)
importance.mitotane.chemo.lvq <- varImp(mitotane.chemo.lvq, scale = FALSE)
print(importance.mitotane.chemo.lvq)
plot(importance.mitotane.chemo.lvq, 20)

#SVM 
mitotane.chemo.svm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
mitotane.chemo.svm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                        data = mitotane.chemo, method = "svmLinear",
                        trControl = mitotane.chemo.svm,
                        preProcess = c("center", "scale"),
                        tuneLength = 10,
                        na.action = na.pass)
importance.mitotane.chemo.svm <- varImp(mitotane.chemo.svm, scale = FALSE)
print(importance.mitotane.chemo.svm)
plot(importance.mitotane.chemo.svm, 20)

#GBM 
mitotane.chemo.gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(86)
mitotane.chemo.gbm <- train(progress..yes.no. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL,
                        data = mitotane.chemo, method = "gbm",
                        trControl = mitotane.chemo.gbm,
                        preProcess = c("center", "scale"),
                        tuneLength = 10,
                        na.action = na.pass)
importance.mitotane.chemo.gbm <- varImp(mitotane.chemo.gbm, scale = FALSE)
print(importance.mitotane.chemo.gbm)
plot(importance.mitotane.chemo.gbm, 20)

#Common predictors for  progress..yes.no.
LVQ.mitotane.chemo <- importance.mitotane.chemo.lvq$importance
SVM.mitotane.chemo <- importance.mitotane.chemo.svm$importance
GBM.mitotane.chemo <- importance.mitotane.chemo.gbm$importance
all(LVQ.mitotane.chemo$X0 == LVQ.mitotane.chemo$X1)
all(SVM.mitotane.chemo$X0 == SVM.mitotane.chemo$X1)
LVQ.mitotane.chemo$X1 <- NULL
SVM.mitotane.chemo$X1 <- NULL
colnames(LVQ.mitotane.chemo)[1] <- "Overall"
colnames(SVM.mitotane.chemo)[1] <- "Overall"
LVQ.mitotane.chemo <- setDT(LVQ.mitotane.chemo, keep.rownames = TRUE)[]
SVM.mitotane.chemo <- setDT(SVM.mitotane.chemo, keep.rownames = TRUE)[]
GBM.mitotane.chemo <- setDT(GBM.mitotane.chemo, keep.rownames = TRUE)[]
LVQ.mitotane.chemo <- LVQ.mitotane.chemo[order(-LVQ.mitotane.chemo$Overall),]
SVM.mitotane.chemo <- SVM.mitotane.chemo[order(-SVM.mitotane.chemo$Overall),]
GBM.mitotane.chemo <- GBM.mitotane.chemo[order(-GBM.mitotane.chemo$Overall),]
LVQ.mitotane.chemo.Sig <- head(LVQ.mitotane.chemo, 20)
SVM.mitotane.chemo.Sig <- head(SVM.mitotane.chemo, 20)
#Adjust GBM selection based on printed variables
GBM.mitotane.chemo.Sig <- head(GBM.mitotane.chemo, 4)
LVQ.mitotane.chemo.Sig$Overall <- NULL
SVM.mitotane.chemo.Sig$Overall <- NULL
GBM.mitotane.chemo.Sig$Overall <- NULL
#The data need to be in a vector. As.vector doesn't work, but unlist reduces data to simple
LVQ.mitotane.chemo.Sig <- unlist(LVQ.mitotane.chemo.Sig)
SVM.mitotane.chemo.Sig <- unlist(SVM.mitotane.chemo.Sig)
GBM.mitotane.chemo.Sig <- unlist(GBM.mitotane.chemo.Sig)
library(VennDiagram)
venn.data.mitotane.chemo <- list(LVQ.mitotane.chemo.Sig, SVM.mitotane.chemo.Sig, GBM.mitotane.chemo.Sig)
grid.newpage()
venn.plot.mitotane.chemo <- venn.diagram(x = list(LVQ.mitotane.chemo.Sig=LVQ.mitotane.chemo.Sig, SVM.mitotane.chemo.Sig=SVM.mitotane.chemo.Sig, GBM.mitotane.chemo.Sig=GBM.mitotane.chemo.Sig),
                                     filename=NULL, 
                                     fill = c("red", "blue", "green"),
                                     alpha = 0.50,
                                     col = "transparent")
grid.draw(venn.plot.mitotane.chemo)
venn.intersect.mitotane.chemo <- calculate.overlap(venn.data.mitotane.chemo)
print(venn.intersect.mitotane.chemo$a5)
print(venn.intersect.mitotane.chemo$a2)

####################################################################################

#Load data
mydata <- read.csv("C:/Users/oisin/Desktop/Analysis Data/SepMutData_Oisin")

mydata$X <- NULL #column numbers appeared as a seperate column
alt.data.whole <- mydata[,c(34:180)]
alt.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(mydata$progress..yes.no.)
alt.data.whole <- cbind(progress..yes.no., alt.data.whole)

#Split df by value in the Chemotherapy column
alt.data <- split(mydata, mydata$Chemotherapy)
alt.data.nochemo <- alt.data$`0`
alt.data.chemo <- alt.data$`1`

#Split nochemo data by mitotane (All forms of therapy)
no.chemo.mito.split <- split(alt.data.nochemo, alt.data.nochemo$Mitotane..pall.or.adj.)
#Data with mitotane only
mitotane.nochemo <- rbind(no.chemo.mito.split$adj, no.chemo.mito.split$`adj+pall`, no.chemo.mito.split$pall)
progress..yes.no. <- mitotane.nochemo$progress..yes.no.
mitotane.nochemo <- mitotane.nochemo[,c(34:180)]
mitotane.nochemo$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.nochemo$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
mitotane.nochemo <- cbind(progress..yes.no., mitotane.nochemo)
#Data which recieved no therapy
no.therapy <- no.chemo.mito.split$`0`
progress..yes.no. <- no.therapy$progress..yes.no.
no.therapy <- no.therapy[,c(34:180)]
no.therapy$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
no.therapy$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
no.therapy <- cbind(progress..yes.no., no.therapy)

#Split chemo data by mitotane (All forms of therapy)
chemo.mito.split <- split(alt.data.chemo, alt.data.chemo$Mitotane..pall.or.adj.)
#Data with both chemotherapy and mitotane therapy
mitotane.chemo <- rbind(chemo.mito.split$adj, chemo.mito.split$`adj+pall`, chemo.mito.split$pall, chemo.mito.split$`adj (with chemo)`)
progress..yes.no. <- mitotane.chemo$progress..yes.no.
mitotane.chemo <- mitotane.chemo[,c(34:180)]
mitotane.chemo$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.chemo$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
mitotane.chemo <- cbind(progress..yes.no., mitotane.chemo)
#Data with chemotherapy only
chemo.nomito <- chemo.mito.split$`0`
progress..yes.no. <- chemo.nomito$progress..yes.no.
chemo.nomito <- chemo.nomito[,c(34:180)]
chemo.nomito$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
chemo.nomito$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
chemo.nomito <- cbind(progress..yes.no., chemo.nomito)


library(corrplot)
#No therapy dataset
corr.no.therapy.sig <- cor(no.therapy[c(1, 16, 19, 30, 38, 43, 46, 51, 52, 64, 66, 68, 77, 78, 88, 94, 97, 114, 129, 136, 141)], no.therapy[c(1, 16, 19, 30, 38, 43, 46, 51, 52, 64, 66, 68, 77, 78, 88, 94, 97, 114, 129, 136, 141)], method = "spearman")
corr.no.therapy.sig
corrplot(corr.no.therapy.sig, type = "upper", method="number", is.corr=FALSE)

#Mitotane only dataset
corr.mitotane.nochemo.sig <- cor(mitotane.nochemo[c(1, 9, 30, 46, 55, 57, 60, 68, 72, 77, 81, 88, 93, 94, 102, 112, 113, 116, 122, 126, 129)], mitotane.nochemo[c(1, 9, 30, 46, 55, 57, 60, 68, 72, 77, 81, 88, 93, 94, 102, 112, 113, 116, 122, 126, 129)], method = "spearman")
corr.mitotane.nochemo.sig
corrplot(corr.mitotane.nochemo.sig, type = "upper", method="number", is.corr=FALSE)

#Mito and Chemo dataset
corr.mitotane.chemo.sig <- cor(mitotane.chemo[c(1, 18, 30, 33, 38, 46, 47, 68, 78, 79, 88, 92, 93, 94, 102, 115, 128, 129, 136, 141, 142, 145)], mitotane.chemo[c(1, 18, 30, 33, 38, 46, 47, 68, 78, 79, 88, 92, 93, 94, 102, 115, 128, 129, 136, 141, 142, 145)], method = "spearman")
corr.mitotane.chemo.sig
corrplot(corr.mitotane.chemo.sig, type = "upper", method="number", is.corr=FALSE)
