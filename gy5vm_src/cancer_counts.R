# counting cancers from linkage, based on source, must hide tumor_restriction <- ... and info_source <- ... line in pilot.R
# properly define NADC_vect
# also year_restriction <- FALSE, remember to switch back when done!!!

rm(list=ls())

#tumors <- "lung"
tumors <- c("anus","lung","prostate","liver","head and neck","HL","bladder","breast","colon and rectum","kidney","oesophagus","skin melanoma","stomach","testis","other NADCs","KSA","NHL","PBL","cervix")
#tumors <- c("Anus & Anal Canal","Bladder","Bones, Joints, Cartilage","Brain & Central Nerves","Breast","Cervix Uteri","Colon, Rectum","Corpus Uteri & NOS","Eye",
#            "Gallbladder & Extrahepatic Bile Tract","Hodgkin Lymphoma","Kaposi Sarcoma","Kidney","Larynx","Liver & Intrahepatic Bile Ducts","Lung, Bronchus, Trachea","Non Hodgkin Lymphoma","Oesophagus",
#            "Oral Cavity & Pharynx","Other Urinary Organs","Ovary","Pancreas","Prostate","Skin Melanoma","Small Intestine","Soft Tissues","Stomach","Testis","Thyroid")

# Pleura: none found -> removed from list of tumors; 
# including KSA category

# also: table(data_l1_res$source) for unrestricted counts (equivalent to year_restriction=FALSE)

count_matrix <- matrix(NA,nrow=length(tumors),ncol=4)

for(i in 1:length(tumors))
{
  tumor_restriction <- tumors[i]
  info_source <- "linkage"
  source("C:\\HomeDir\\P3RL\\R_code\\pilot.R")
  id_link <- X1$shcsid
  info_source <- "SHCS"
  source("C:\\HomeDir\\P3RL\\R_code\\pilot.R")
  id_shcs <- X1$shcsid
  rm(X1)
  count_matrix[i,1] <- length(intersect(id_link,id_shcs))
  count_matrix[i,2] <- length(setdiff(id_link,id_shcs))
  count_matrix[i,3] <- length(setdiff(id_shcs,id_link))
  count_matrix[i,4] <- length(union(id_link,id_shcs))
}

colnames(count_matrix) <- c("both","linkage only","SHCS only","either")
row.names(count_matrix) <- tumors

write.table(count_matrix,file=paste(filepath_code,"\\cancer_counts.txt",sep=''),sep=",")