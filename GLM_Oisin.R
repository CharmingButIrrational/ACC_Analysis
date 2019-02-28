#Load data
mydata <- read.csv("C:/Users/User/Desktop/Analysis Data/SepMutData_Oisin")
mydata$X <- NULL

library(randomForest)
library(caret)

smp_size <- floor(0.8 * nrow(mydata))
set.seed(345)
train_ind <- sample(seq_len(nrow(mydata)), size = smp_size)

alt.data <- mydata[,c(34:180)]
alt.data$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- mydata$progress..yes.no. 
alt.data <- cbind(progress..yes.no., alt.data)

alt.data[["progress..yes.no."]] = factor(alt.data[["progress..yes.no."]])

train.data <- alt.data[train_ind,]
test.data <- alt.data[-train_ind,]

glm.gamma.dist <- glm(PFS..months. ~ A_EGFR + A_ATM + A_FLCN + A_CDK4 + A_TERT + 
                        A_PBRM1 + A_MAP3K1 + A_ATRX + A_PPP2R1A + A_SMARCB1 + 
                        A_ABL1 + A_PALB2 + A_MYC + A_MED12 + A_DAXX + A_CSFR1 + 
                        A_U2AF1 + A_CREBBP + A_FBXO11 + A_CDH1 + A_ZNFR3 + 
                        A_ARID1A + A_GRIN2A + A_BCL6 + A_TSC2 + A_BRCA1 + 
                        A_CDK12 + A_PTCH1 + A_TP53 + A_EZH2 + A_MSH2 + A_AR + 
                        A_GNA11 + A_SMAD4 + A_BRCA2 + A_RB1 + A_CTNNB1 + A_KIT + 
                        A_SETD2 + A_KDR + A_JAK1 + A_WT1 + A_SLC7A8 + A_KMT2D + 
                        A_NF1 + A_NOTCH1 + A_KDM6A + A_CIC + A_PRDM1 + A_JAK2 + 
                        A_RET + A_BCOR + A_ESR1 + A_APC + A_FBXW7 + A_GNAS + 
                        A_MEN1 + A_DNMT3A + A_ZNRF3 + A_FANCE + A_TSC1 + 
                        A_SMARC4A + A_MSH6 + A_EP300 + A_PRKAR1A + A_BUB1B + 
                        A_none + A_MLH1 + A_H3F3A + A_SMO + A_BRAF + B_APC + 
                        B_SMARC4 + B_BAP1 + B_NPM1 + B_FGFR3 + B_CIC + B_GNA11 + 
                        B_IL7R + B_GNAS + B_KMT2D + B_MAP3K1 + B_ZRSR2 + B_DDB2 + 
                        B_MDM2 + B_MEN1 + B_none + B_PIK3R1 + B_ARID2 + B_ABL1 + 
                        B_NOTCH1 + B_TERT + B_CDK4 + B_CSF1R + B_VHL + B_MYC + 
                        B_TSC1 + B_PRKAR1A + B_DDR2 + B_SLC7A8 + B_STK11 + 
                        B_IL6ST + B_U2AF1 + B_ASXL1 + C_CARD11 + C_SDHB + 
                        C_TNRFSF14 + C_CYLD + C_SLC7A8 + C_PHF6 + C_FGFR2 + 
                        C_SPOP + C_MLH1 + C_none + C_BRIP1 + C_MSH2 + C_CD79 + 
                        C_SMO + C_GNAQ + C_ECT2L + C_CDC73 + C_PIK3R1 + C_TNFSRF14 + 
                        C_BRAF + C_RB1 + C_SMAD4 + C_MAP2K4 + C_SMARCB1 + C_DNMT3A + 
                        C_BRCA2 + C_JAK1 + C_EZH2 + C_FAS + C_NF1 + C_CHEK2 + 
                        C_PALB2 + C_CSF1R + C_TNFSRSF14 + C_CDK12 + C_FLCN + 
                        C_FLNC + C_MAP4K3 + C_EPCAM + C_TNFRSF14 + C_VHL, family = "Gamma"(link=log), data = train.data)

predict(glm.gamma.dist, )

step.glm.gamma.both <- stepAIC(glm.gamma.dist, direction = "both", trace = FALSE)

summary(test.data, step.glm.gamma.both)
                            