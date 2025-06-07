

data_pat <- read.csv(file=paste(filepath_read_SHCS,"\\csv\\pat.csv",sep=''),sep=",", stringsAsFactors=FALSE)
names(data_pat)[1] <- "shcsid"

data_lab <- read.csv(file=paste(filepath_read_SHCS,"\\csv\\lab.csv",sep=''),sep=",",na.strings="",stringsAsFactors=FALSE)
names(data_lab)[1] <- "shcsid"
data_drug <- read.csv(file=paste(filepath_read_SHCS,"\\csv\\drug.csv",sep=''),sep=",", stringsAsFactors=FALSE)
names(data_drug)[1] <- "shcsid"
data_dis <- read.csv(file=paste(filepath_read_SHCS,"\\csv\\dis.csv",sep=''),sep=",", stringsAsFactors=FALSE)
names(data_dis)[1] <- "shcsid"
data_admi <- read.csv(file=paste(filepath_read_SHCS,"\\csv\\admi.csv",sep=''),sep=",", stringsAsFactors=FALSE)
names(data_admi)[1] <- "shcsid"
data_tail <-  read.csv(file=paste(filepath_read_SHCS,"\\csv\\tail.csv",sep=''),sep=",", stringsAsFactors=FALSE)
names(data_tail)[1] <- "shcsid"
id_registry_table <- read.dta(paste(filepath_read_SHCS,"\\stata\\lookup_hiv_origin.dta",sep=""))
data_gyn <-  read.csv(file=paste(filepath_read_SHCS,"\\csv\\gyn.csv",sep=''),sep=",", stringsAsFactors=FALSE)
names(data_gyn)[1] <- "shcsid"
data_cvri <-  read.csv(file=paste(filepath_read_SHCS,"\\csv\\cvri.csv",sep=''),sep=",", stringsAsFactors=FALSE)
names(data_cvri)[1] <- "shcsid"

id_vect <- data_pat$shcsid

data_pat$registry <- id_registry_table$hiv_registry1[match(data_pat$shcsid,id_registry_table$shcsid)]

data_pat$registry2 <- NA
data_pat$registry2[(id_vect>=10000 & id_vect<20000) | (id_vect>=90000 & id_vect <99900)] <- "ZH"
data_pat$registry2[(id_vect>=20000 & id_vect<30000)] <- "BS"
data_pat$registry2[(id_vect>=30000 & id_vect<40000)] <- "BE"
data_pat$registry2[(id_vect>=40000 & id_vect<46000) | (id_vect>=48000 & id_vect<50000)] <- "GE"
data_pat$registry2[(id_vect>=50000 & id_vect<60000)] <- "VD"
data_pat$registry2[(id_vect>=60000 & id_vect<70000)] <- "TI"
data_pat$registry2[(id_vect>=46000 & id_vect<48000)] <- "SG"

center_id <- data_pat$center
data_pat$registry3 <- NA
data_pat$registry3[center_id==10] <- "ZH"
data_pat$registry3[center_id==20] <- "BS"
data_pat$registry3[center_id==30] <- "BE"
data_pat$registry3[center_id==40] <- "GE"
data_pat$registry3[center_id==50] <- "VD"
data_pat$registry3[center_id==60] <- "TI"
data_pat$registry3[center_id==70] <- "SG"

#formatting dates
Sys.setlocale("LC_TIME", "C")
data_pat$regdate <- with(data_pat,as.Date(paste(substr(regdate,1,7),ifelse(as.numeric(substr(regdate,8,9))<50,"20","19"),substr(regdate,8,9),sep=''),"%d %b %Y"))
data_pat$hiv_posdate <- with(data_pat,as.Date(paste(substr(hiv_posdate,1,7),ifelse(as.numeric(substr(hiv_posdate,8,9))<50,"20","19"),substr(hiv_posdate,8,9),sep=''),"%d %b %Y"))
data_pat$hiv_negdate <- with(data_pat,as.Date(paste(substr(hiv_negdate,1,7),ifelse(as.numeric(substr(hiv_negdate,8,9))<50,"20","19"),substr(hiv_negdate,8,9),sep=''),"%d %b %Y"))
data_lab$labdate <- with(data_lab,as.Date(paste(substr(labdate,1,7),ifelse(as.numeric(substr(labdate,8,9))<50,"20","19"),substr(labdate,8,9),sep=''),"%d %b %Y"))
data_lab$cd4date <- with(data_lab,as.Date(paste(substr(cd4date,1,7),ifelse(as.numeric(substr(cd4date,8,9))<50,"20","19"),substr(cd4date,8,9),sep=''),"%d %b %Y"))
data_lab$heb_date <- with(data_lab,as.Date(paste(substr(heb_date,1,7),ifelse(as.numeric(substr(heb_date,8,9))<50,"20","19"),substr(heb_date,8,9),sep=''),"%d %b %Y"))
data_lab$hbvdna_date <- with(data_lab,as.Date(paste(substr(hbvdna_date,1,7),ifelse(as.numeric(substr(hbvdna_date,8,9))<50,"20","19"),substr(hbvdna_date,8,9),sep=''),"%d %b %Y"))
data_lab$hec_date <- with(data_lab,as.Date(paste(substr(hec_date,1,7),ifelse(as.numeric(substr(hec_date,8,9))<50,"20","19"),substr(hec_date,8,9),sep=''),"%d %b %Y"))
data_lab$hcv_gen_date <- with(data_lab,as.Date(paste(substr(hcv_gen_date,1,7),ifelse(as.numeric(substr(hcv_gen_date,8,9))<50,"20","19"),substr(hcv_gen_date,8,9),sep=''),"%d %b %Y"))
data_lab$hcv_rna_date <- with(data_lab,as.Date(paste(substr(hcv_rna_date,1,7),ifelse(as.numeric(substr(hcv_rna_date,8,9))<50,"20","19"),substr(hcv_rna_date,8,9),sep=''),"%d %b %Y"))
data_drug$drstartd <- with(data_drug,as.Date(paste(substr(drstartd,1,7),ifelse(as.numeric(substr(drstartd,8,9))<50,"20","19"),substr(drstartd,8,9),sep=''),"%d %b %Y"))
data_drug$drstopd <- with(data_drug,as.Date(paste(substr(drstopd,1,7),ifelse(as.numeric(substr(drstopd,8,9))<50,"20","19"),substr(drstopd,8,9),sep=''),"%d %b %Y"))
data_dis$newdate <- with(data_dis,as.Date(paste(substr(newdate,1,7),ifelse(as.numeric(substr(newdate,8,9))<50,"20","19"),substr(newdate,8,9),sep=''),"%d %b %Y"))
data_admi$exitdate <- with(data_admi,as.Date(paste(substr(exitdate,1,7),ifelse(as.numeric(substr(exitdate,8,9))<50,"20","19"),substr(exitdate,8,9),sep=''),"%d %b %Y"))
data_admi$stopdate <- with(data_admi,as.Date(paste(substr(stopdate,1,7),ifelse(as.numeric(substr(stopdate,8,9))<50,"20","19"),substr(stopdate,8,9),sep=''),"%d %b %Y"))
data_admi$last_info_date <- with(data_admi,as.Date(paste(substr(last_info_date,1,7),ifelse(as.numeric(substr(last_info_date,8,9))<50,"20","19"),substr(last_info_date,8,9),sep=''),"%d %b %Y"))
data_tail$art_start_date <- with(data_tail,as.Date(paste(substr(art_start_date,1,7),ifelse(as.numeric(substr(art_start_date,8,9))<50,"20","19"),substr(art_start_date,8,9),sep=''),"%d %b %Y"))
data_tail$haart_start_date <- with(data_tail,as.Date(paste(substr(haart_start_date,1,7),ifelse(as.numeric(substr(haart_start_date,8,9))<50,"20","19"),substr(haart_start_date,8,9),sep=''),"%d %b %Y"))
data_tail$rna_first_date <- with(data_tail,as.Date(paste(substr(rna_first_date,1,7),ifelse(as.numeric(substr(rna_first_date,8,9))<50,"20","19"),substr(rna_first_date,8,9),sep=''),"%d %b %Y"))
data_tail$first_b_date <- with(data_tail,as.Date(paste(substr(first_b_date,1,7),ifelse(as.numeric(substr(first_b_date,8,9))<50,"20","19"),substr(first_b_date,8,9),sep=''),"%d %b %Y"))
data_tail$first_c_date <- with(data_tail,as.Date(paste(substr(first_c_date,1,7),ifelse(as.numeric(substr(first_c_date,8,9))<50,"20","19"),substr(first_c_date,8,9),sep=''),"%d %b %Y"))
data_gyn$gyndate <- with(data_gyn,as.Date(paste(substr(gyndate,1,7),ifelse(as.numeric(substr(gyndate,8,9))<50,"20","19"),substr(gyndate,8,9),sep=''),"%d %b %Y"))
data_gyn <- data_gyn[order(data_gyn$shcsid,data_gyn$gyndate),]
data_cvri$cardiodate <- with(data_cvri,as.Date(paste(substr(cardiodate,1,7),ifelse(as.numeric(substr(cardiodate,8,9))<50,"20","19"),substr(cardiodate,8,9),sep=''),"%d %b %Y"))
data_cvri <- data_cvri[order(data_cvri$shcsid,data_cvri$cardiodate),]

#only keeping earliest reported dates for each patient
data_lab_res <- data_lab[!duplicated(data_lab$shcsid),]
data_drug_res <- data_drug[!duplicated(data_drug$shcsid),]
data_dis_res <- data_dis[order(data_dis$shcsid,data_dis$newdate),]
data_dis_res <- data_dis_res[!duplicated(data_dis_res$shcsid),]

#only keeping latest reported dates for each patient
data_lab_res_dcr <- data_lab[!duplicated(data_lab$shcsid,fromLast=TRUE),]
data_drug_res_dcr <- data_drug[!duplicated(data_drug$shcsid,fromLast=TRUE),]
data_dis_res_dcr <- data_dis[order(data_dis$shcsid,data_dis$newdate),]
data_dis_res_dcr <- data_dis_res_dcr[!duplicated(data_dis_res_dcr$shcsid,fromLast=TRUE),]
