library(GEOquery)
library(VennDiagram)

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
selected1b <- eset1b[select1b]

#eset2b
select2b <- c("title", "geo_accession", "status", "submission_date", "last_update_date", "type", 
              "channel_count", "source_name_ch1", "organism_ch1", "characteristics_ch1", 
              "treatment_protocol_ch1", "molecule_ch1", "extract_protocol_ch1", "label_ch1", 
              "label_protocol_ch1", "taxid_ch1", "hyb_protocol", "scan_protocol", "data_processing", 
              "platform_id", "contact_name", "contact_department", "contact_institute", 
              "contact_address", "contact_city", "contact_zip/postal_code", "contact_country", 
              "supplementary_file", "data_row_count")
selected2b <- eset2b[select2b]

#eset3b
select3b <- c("title", "geo_accession", "status", "submission_date", "last_update_date", "type", 
              "channel_count", "source_name_ch1", "organism_ch1", "characteristics_ch1", 
              "treatment_protocol_ch1", "molecule_ch1", "extract_protocol_ch1", "label_ch1", 
              "label_protocol_ch1", "taxid_ch1", "hyb_protocol", "scan_protocol", "data_processing", 
              "platform_id", "contact_name", "contact_department", "contact_institute", 
              "contact_address", "contact_city", "contact_zip/postal_code", "contact_country", 
              "supplementary_file", "data_row_count")
selected3b <- eset3b[select3b]
