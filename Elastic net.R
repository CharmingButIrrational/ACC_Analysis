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
alt.data.nochemo <- alt.data0[,c(34:180)]
alt.data.nochemo$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data.nochemo$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(alt.data0$progress..yes.no.)
alt.data.nochemo <- cbind(progress..yes.no., alt.data.nochemo)
#For treatment with chemotherapy 
alt.data.chemo <- alt.data1[,c(34:180)]
alt.data.chemo$gene.with.CN.gains..80..log..FC.0.5..1.5x..amplicons.min.5. <- NULL
alt.data.chemo$gene.with.CN.losses..80..log..FC.0.25..amplicons.min.5. <- NULL
progress..yes.no. <- as.factor(alt.data1$progress..yes.no.)
alt.data.chemo <- cbind(progress..yes.no., alt.data11)

#Seperating Chemo dataset into A_, B_, C_
progress..yes.no. <- as.factor(alt.data.chemo$progress..yes.no.)
alt.data.chemo.A <- alt.data.chemo[,c(1:72)]
alt.data.chemo.B <- alt.data.chemo[,c(73:105)]
alt.data.chemo.B <- cbind(progress..yes.no., alt.data.chemo.B)
alt.data.chemo.C <- alt.data.chemo[,c(105:146)]
alt.data.chemo.C <- cbind(progress..yes.no., alt.data.chemo.C)

#Seperating NoChemo dataset into A_, B_, C_
progress..yes.no. <- as.factor(alt.data.nochemo$progress..yes.no.)
alt.data.nochemo.A <- alt.data00[,c(1:72)]
alt.data.nochemo.B <- alt.data.nochemo[,c(73:105)]
alt.data.nochemo.B <- cbind(progress..yes.no., alt.data.nochemo.B)
alt.data.nochemo.C <- alt.data.nochemo[,c(105:146)]
alt.data.nochemo.C <- cbind(progress..yes.no., alt.data.nochemo.C)

######################################################################################################
set.seed(233)
training.set <- alt.data.whole$progress..yes.no. %>%
  createDataPartition(p = 0.75, list = FALSE)
train.data.whole  <- alt.data.whole[training.set, ]
test.data.whole <- alt.data.whole[-training.set, ]
# Whole dataset
# Predictor variables
elasticnet.whole <- model.matrix(progress..yes.no.~., train.data.whole)[,-1]
# Outcome variable
y.whole <- train.data.whole$progress..yes.no.

elasticnet.whole <- cv.glmnet(elasticnet.whole, y.whole, alpha = 0.5)