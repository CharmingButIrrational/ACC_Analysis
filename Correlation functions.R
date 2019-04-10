#Load data
mydata <- read.csv("C:/Users/User/Desktop/Analysis Data/SepMutData_Oisin")

mydata$X <- NULL #column numbers appeared as a seperate column
alt.data.whole <- mydata[,c(34:180)]
alt.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- mydata$progress..yes.no.
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
progress..yes.no. <- alt.data0$progress..yes.no.
alt.data.nochemo <- cbind(progress..yes.no., alt.data00)
#For treatment with chemotherapy 
alt.data11 <- alt.data1[,c(34:180)]
alt.data11$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data11$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- alt.data1$progress..yes.no.
alt.data.chemo <- cbind(progress..yes.no., alt.data11)

#Seperating Chemo dataset into A_, B_, C_
progress..yes.no. <- as.factor(alt.data1$progress..yes.no.)
alt.data11.A <- alt.data.chemo[,c(1:72)]
alt.data11.B <- alt.data.chemo[,c(73:105)]
alt.data11.B <- cbind(progress..yes.no., alt.data11.B)
alt.data11.C <- alt.data.chemo[,c(105:146)]
alt.data11.C <- cbind(progress..yes.no., alt.data11.C)

#Seperating NoChemo dataset into A_, B_, C_
progress..yes.no. <- as.factor(alt.data0$progress..yes.no.)
alt.data00.A <- alt.data.nochemo[,c(1:72)]
alt.data00.B <- alt.data.nochemo[,c(73:105)]
alt.data00.B <- cbind(progress..yes.no., alt.data00.B)
alt.data00.C <- alt.data.nochemo[,c(105:146)]
alt.data00.C <- cbind(progress..yes.no., alt.data00.C)

#Split df by value in the Mitotane column
mitotane.data <- split(mydata, mydata$Mitotane..pall.or.adj.)

mitotane.pos <- rbind(mitotane.data$adj, mitotane.data$`adj (with chemo)`, mitotane.data$`adj+pall`, mitotane.data$pall)
mitotane.neg <- rbind(mitotane.data$`0`)

mitotane.pos.data.whole <- mitotane.pos[,c(34:180)]
mitotane.pos.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.pos.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- mitotane.pos$progress..yes.no.
mitotane.pos.data.whole <- cbind(progress..yes.no., mitotane.pos.data.whole)

mitotane.neg.data.whole <- mitotane.neg[,c(34:180)]
mitotane.neg.data.whole$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
mitotane.neg.data.whole$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- mitotane.neg$progress..yes.no.
mitotane.neg.data.whole <- cbind(progress..yes.no., mitotane.neg.data.whole)

#Seperating Mitotane dataset into A_, B_, C_
progress..yes.no. <- mitotane.pos.data.whole$progress..yes.no.
mitotane.pos.data.A <- mitotane.pos.data.whole[,c(1:72)]
mitotane.pos.data.B <- mitotane.pos.data.whole[,c(73:105)]
mitotane.pos.data.B <- cbind(progress..yes.no., mitotane.pos.data.B)
mitotane.pos.data.C <- mitotane.pos.data.whole[,c(105:146)]
mitotane.pos.data.C <- cbind(progress..yes.no., mitotane.pos.data.C)

#Seperating NoMitotane dataset into A_, B_, C_
progress..yes.no. <- mitotane.neg.data.whole$progress..yes.no.
mitotane.neg.data.A <- mitotane.neg.data.whole[,c(1:72)]
mitotane.neg.data.B <- mitotane.neg.data.whole[,c(73:105)]
mitotane.neg.data.B <- cbind(progress..yes.no., mitotane.neg.data.B)
mitotane.neg.data.C <- mitotane.neg.data.whole[,c(105:146)]
mitotane.neg.data.C <- cbind(progress..yes.no., mitotane.neg.data.C)

library(corrplot)

corr.whole <- cor(alt.data.whole[-1], alt.data.whole$progress..yes.no., method = "spearman")
corr.whole
corr.whole.df <- as.data.frame(corr.whole)
corr.whole.order <- corr.whole.df[order(-corr.whole.df$V1)]

corr.whole.A <- cor(alt.data.whole.A[-1], alt.data.whole.A$progress..yes.no., method = "spearman")
corr.whole.A
corr.whole.df.A <- as.data.frame(corr.whole.A)
corr.whole.order.A <- corr.whole.df.A[order(-corr.whole.df.A$V1)]
corr.whole.B <- cor(alt.data.whole.B[-1], alt.data.whole.B$progress..yes.no., method = "spearman")
corr.whole.B
corr.whole.df.B <- as.data.frame(corr.whole.B)
corr.whole.order.B <- corr.whole.df.B[order(-corr.whole.df.B$V1)]
corr.whole.C <- cor(alt.data.whole.C[-1], alt.data.whole.C$progress..yes.no., method = "spearman")
corr.whole.C 
corr.whole.df.C <- as.data.frame(corr.whole.C)
corr.whole.order.C <- corr.whole.df.C[order(-corr.whole.df.C$V1)]


corr.chemo <- cor(alt.data.chemo[-1], alt.data.chemo$progress..yes.no., method = "spearman")
corr.chemo 
corr.chemo.df <- as.data.frame(corr.chemo)
corr.chemo.order <- corr.chemo.df[order(-corr.chemo.df$V1)]

corr.chemo.A <- cor(alt.data.whole.A[-1], alt.data.whole.A$progress..yes.no., method = "spearman")
corr.chemo.A
corr.chemo.df.A <- as.data.frame(corr.chemo.A)
corr.chemo.order.A <- corr.chemo.df.A[order(-corr.chemo.df.A$V1)]
corr.chemo.B <- cor(alt.data.chemo.B[-1], alt.data.chemo.B$progress..yes.no., method = "spearman")
corr.chemo.B
corr.chemo.df.B <- as.data.frame(corr.chemo.B)
corr.chemo.order.B <- corr.chemo.df.B[order(-corr.chemo.df.B$V1)]
corr.chemo.C <- cor(alt.data.chemo.C[-1], alt.data.chemo.C$progress..yes.no., method = "spearman")
corr.chemo.C 
corr.chemo.df.C <- as.data.frame(corr.chemo.C)
corr.chemo.order.C <- corr.chemo.df.C[order(-corr.chemo.df.C$V1)]


corr.nochemo <- cor(alt.data.nochemo[-1], alt.data.nochemo$progress..yes.no., method = "spearman")
corr.nochemo 
corr.nochemo.df <- as.data.frame(corr.nochemo)
corr.nochemo.order <- corr.nochemo.df[order(-corr.nochemo.df$V1)]

corr.nochemo.A <- cor(alt.data.nochemo.A[-1], alt.data.nochemo.A$progress..yes.no., method = "spearman")
corr.nochemo.A
corr.nochemo.df.A <- as.data.frame(corr.nochemo.A)
corr.nochemo.order.A <- corr.nochemo.df.A[order(-corr.nochemo.df.A$V1)]
corr.nochemo.B <- cor(alt.data.nochemo.B[-1], alt.data.nochemo.B$progress..yes.no., method = "spearman")
corr.nochemo.B
corr.nochemo.df.B <- as.data.frame(corr.nochemo.B)
corr.nochemo.order.B <- corr.nochemo.df.B[order(-corr.nochemo.df.B$V1)]
corr.nochemo.C <- cor(alt.data.nochemo.C[-1], alt.data.nochemo.C$progress..yes.no., method = "spearman")
corr.nochemo.C 
corr.nochemo.df.C <- as.data.frame(corr.nochemo.C)
corr.nochemo.order.C <- corr.nochemo.df.C[order(-corr.nochemo.df.C$V1)]


corr.Mito <- cor(mitotane.pos.data.whole[-1], mitotane.pos.data.whole$progress..yes.no., method = "spearman")
corr.Mito 
corr.Mito.df <- as.data.frame(corr.Mito)
corr.Mito.order <- corr.Mito.df[order(-corr.Mito.df$V1)]

corr.Mito.A <- cor(mitotane.pos.data.A[-1], mitotane.pos.data.A$progress..yes.no., method = "spearman")
corr.Mito.A
corr.Mito.df.A <- as.data.frame(corr.Mito.A)
corr.Mito.order.A <- corr.Mito.df.A[order(-corr.Mito.df.A$V1)]
corr.Mito.B <- cor(mitotane.pos.data.B[-1], mitotane.pos.data.B$progress..yes.no., method = "spearman")
corr.Mito.B
corr.Mito.df.B <- as.data.frame(corr.Mito.B)
corr.Mito.order.B <- corr.Mito.df.B[order(-corr.Mito.df.B$V1)]
corr.Mito.C <- cor(mitotane.pos.data.C[-1], mitotane.pos.data.C$progress..yes.no., method = "spearman")
corr.Mito.C
corr.Mito.df.C <- as.data.frame(corr.Mito.C)
corr.Mito.order.C <- corr.Mito.df.C[order(-corr.Mito.df.C$V1)]


corr.noMito <- cor(mitotane.neg.data.whole[-1], mitotane.neg.data.whole$progress..yes.no., method = "spearman")
corr.noMito
corr.noMito.df <- as.data.frame(corr.noMito)
corr.noMito.order <- corr.noMito.df[order(-corr.noMito.df$V1)]

corr.noMito.A <- cor(mitotane.neg.data.A[-1], mitotane.neg.data.A$progress..yes.no., method = "spearman")
corr.noMito.A
corr.noMito.df.A <- as.data.frame(corr.noMito.A)
corr.noMito.order.A <- corr.noMito.df.A[order(-corr.noMito.df.A$V1)]
corr.noMito.B <- cor(mitotane.neg.data.B[-1], mitotane.neg.data.B$progress..yes.no., method = "spearman")
corr.noMito.B
corr.noMito.df.B <- as.data.frame(corr.noMito.B)
corr.noMito.order.B <- corr.noMito.df.B[order(-corr.noMito.df.B$V1)]
corr.noMito.C <- cor(mitotane.neg.data.C[-1], mitotane.neg.data.C$progress..yes.no., method = "spearman")
corr.noMito.C 
corr.noMito.df.C <- as.data.frame(corr.noMito.C)
corr.noMito.order.C <- corr.noMito.df.C[order(-corr.noMito.df.C$V1)]

#Whole dataset
corrplot(corr.whole, method="number", is.corr=FALSE)
#Chemo dataset
corrplot(corr.chemo, method="number", is.corr=FALSE)
#NoChemo dataset
corrplot(corr.nochemo, method="number", is.corr=FALSE)
#Mito dataset
corrplot(mitotane.pos.data.whole, method="number", is.corr=FALSE)
#NoMito dataset
corrplot(mitotane.neg.data.whole, method="number", is.corr=FALSE)


#Separated whole data
corrplot(corr.whole.A, method="number", is.corr=FALSE)
corrplot(corr.whole.B, method="number", is.corr=FALSE)
corrplot(corr.whole.C, method="number", is.corr=FALSE)

#Separated chemo data
corrplot(corr.chemo.A, method="number", is.corr=FALSE)
corrplot(corr.chemo.B, method="number", is.corr=FALSE)
corrplot(corr.chemo.C, method="number", is.corr=FALSE)

#Separated nochemo data
corrplot(corr.nochemo.A, method="number", is.corr=FALSE)
corrplot(corr.nochemo.B, method="number", is.corr=FALSE)
corrplot(corr.nochemo.C, method="number", is.corr=FALSE)

#Separated mitotane data
corrplot(mitotane.pos.data.A, method="number", is.corr=FALSE)
corrplot(mitotane.pos.data.B, method="number", is.corr=FALSE)
corrplot(mitotane.pos.data.C, method="number", is.corr=FALSE)

#Separated no mitodata data
corrplot(mitotane.neg.data.A, method="number", is.corr=FALSE)
corrplot(mitotane.neg.data.B, method="number", is.corr=FALSE)
corrplot(mitotane.neg.data.C, method="number", is.corr=FALSE)