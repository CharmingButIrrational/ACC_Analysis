library(GEOquery)
library(VennDiagram)
library(dplyr)



GSE106358 <- getGEO("GSE106358") #Txt type file
GSE94340 <- getGEO("GSE94340") #CEL type file
GSE9285 <- getGEO("GSE9285") #GPR type file 

eset1 = GSE106358[[1]]
eset2 = GSE94340[[1]]
eset3 = GSE9285[[1]]

eset1b=eset1
eset2b=eset2
eset3b=eset3

dim(pData(eset1b))
dim(pData(eset2b))
dim(pData(eset3b))

colnames(pData(eset1b))
colnames(pData(eset2b))
colnames(pData(eset3b))

pData_eset1b <- colnames(pData(eset1b))
pData_eset2b <- colnames(pData(eset2b))
pData_eset3b <- colnames(pData(eset3b))

#Attempt to find the same columns in the different dfs automatically

list <- list(pData_eset1b, pData_eset2b, pData_eset3b)

overlap <- calculate.overlap(list)
overlap$a5 #the $a5 will be different if there are different numbers of datasets in the list

paste(overlap$a5, collapse = ", ") #Prints out the common columns

#####How to decide which columns are important
#eset1b
select1b <- c("title", "geo_accession", "status", "submission_date", "last_update_date", "type", 
            "channel_count", "source_name_ch1", "organism_ch1", "characteristics_ch1", 
            "treatment_protocol_ch1", "molecule_ch1", "extract_protocol_ch1", "label_ch1", 
            "label_protocol_ch1", "taxid_ch1", "hyb_protocol", "scan_protocol", "data_processing", 
            "platform_id", "contact_name", "contact_department", "contact_institute", 
            "contact_address", "contact_city", "contact_zip/postal_code", "contact_country", 
            "supplementary_file", "data_row_count")

setdiff(colnames(pData(eset1b)), select1b)

eset1b$growth_protocol_ch1 <- NULL
eset1b$source_name_ch2 <- NULL
eset1b$organism_ch2 <- NULL
eset1b$characteristics_ch2 <- NULL
eset1b$characteristics_ch2.1<- NULL
eset1b$characteristics_ch2.2 <- NULL
eset1b$characteristics_ch2.3 <- NULL
eset1b$characteristics_ch2.4 <- NULL
eset1b$characteristics_ch2.5 <- NULL
eset1b$characteristics_ch2.6 <- NULL
eset1b$characteristics_ch2.7 <- NULL
eset1b$treatment_protocol_ch2 <- NULL
eset1b$growth_protocol_ch2 <- NULL
eset1b$molecule_ch2 <- NULL
eset1b$extract_protocol_ch2 <- NULL
eset1b$label_ch2 <- NULL
eset1b$label_protocol_ch2 <- NULL
eset1b$taxid_ch2 <- NULL
eset1b$description <- NULL
eset1b$description.1 <- NULL
eset1b$contact_email <- NULL
eset1b$age:ch2 <- NULL
eset1b$'disease state:ch2' <- NULL
eset1b$gender:ch2 <- NULL
eset1b$mrss:ch2 <- NULL
eset1b$'sample type:ch1' <- NULL
eset1b$subject:ch2 <- NULL
eset1b$`time point:ch2` <- NULL
eset1b$tissue:ch2 <- NULL
eset1b$treatment:ch2 <- NULL


subjectseset1b <- list("GSM2836412", "GSM2836413", "GSM2836416",
                      "GSM2836418", "GSM2836419", "GSM2836420",
                      "GSM2836421", "GSM2836423", "GSM2836426",
                      "GSM2836428", "GSM2836430", "GSM2836431",
                      "GSM2836432", "GSM2836434", "GSM2836436",
                      "GSM2836438", "GSM2836440", "GSM2836442",
                      "GSM2836444", "GSM2836446", "GSM2836448",
                      "GSM2836451", "GSM2836452", "GSM2836454",
                      "GSM2836456", "GSM2836457", "GSM2836459",
                      "GSM2836460", "GSM2836462", "GSM2836464",
                      "GSM2836466", "GSM2836468", "GSM2836470",
                      "GSM2836471", "GSM2836473", "GSM2836475",
                      "GSM2836476", "GSM2836478", "GSM2836479")

#This command should work to subset for specific samples (e.g GSM2836412)
eset1btest = eset1b[, colnames(eset1b) %in% subjectseset1b]


#eset2b
select2b <- c("title", "geo_accession", "status", "submission_date", "last_update_date", "type", 
              "channel_count", "source_name_ch1", "organism_ch1", "characteristics_ch1", 
              "treatment_protocol_ch1", "molecule_ch1", "extract_protocol_ch1", "label_ch1", 
              "label_protocol_ch1", "taxid_ch1", "hyb_protocol", "scan_protocol", "data_processing", 
              "platform_id", "contact_name", "contact_department", "contact_institute", 
              "contact_address", "contact_city", "contact_zip/postal_code", "contact_country", 
              "supplementary_file", "data_row_count")

setdiff(colnames(pData(eset2b)), select2b)

eset2b$characteristics_ch1.1 <- NULL
eset2b$contact_state <- NULL
eset2b$timepoint:ch1 <- NULL
eset2b$treatment:ch1 <- NULL

#eset3b
select3b <- c("title", "geo_accession", "status", "submission_date", "last_update_date", "type", 
              "channel_count", "source_name_ch1", "organism_ch1", "characteristics_ch1", 
              "treatment_protocol_ch1", "molecule_ch1", "extract_protocol_ch1", "label_ch1", 
              "label_protocol_ch1", "taxid_ch1", "hyb_protocol", "scan_protocol", "data_processing", 
              "platform_id", "contact_name", "contact_department", "contact_institute", 
              "contact_address", "contact_city", "contact_zip/postal_code", "contact_country", 
              "supplementary_file", "data_row_count")

setdiff(colnames(pData(eset3b)), select3b)

eset3b$growth_protocol_ch1 <- NULL   
eset3b$source_name_ch2 <- NULL
eset3b$organism_ch2 <- NULL
eset3b$characteristics_ch2 <- NULL
eset3b$treatment_protocol_ch2 <- NULL
eset3b$growth_protocol_ch2 <- NULL
eset3b$molecule_ch2 <- NULL
eset3b$extract_protocol_ch2 <- NULL
eset3b$label_ch2 <- NULL
eset3b$label_protocol_ch2 <- NULL
eset3b$taxid_ch2 <- NULL
eset3b$description <- NULL
eset3b$contact_email <- NULL
eset3b$contact_laboratory <- NULL
eset3b$contact_state <- NULL
eset3b$relation <- NULL
eset3b$relation.1 <- NULL

##### 3. Checking normalization #####

boxplot(exprs())
boxplot(exprs())
boxplot(exprs())


#overview of how many tumor and cancer samples are in esets
table(pData(eset1b)[,"Tumor"])
table(pData(eset2b)[,"Tumor"])
table(pData(eset3b)[,"Tumor"])

