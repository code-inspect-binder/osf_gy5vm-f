# appending patient information from SHCS to linkage data

pat_to_l0 <- match(data_l0_res$shcsid,data_pat$shcsid)
pat_to_l1 <- match(data_l1_res$shcsid,data_pat$shcsid)

# first positive/last negative/SHCS registration date
negdate_l0 <- data_pat$hiv_negdate[pat_to_l0]
posdate_l0 <- data_pat$hiv_posdate[pat_to_l0]
negdate_l1 <- data_pat$hiv_negdate[pat_to_l1]
posdate_l1 <- data_pat$hiv_posdate[pat_to_l1]
regdate_l0 <- data_pat$regdate[pat_to_l0]
regdate_l1 <- data_pat$regdate[pat_to_l1]
#earliest/latest lab dates for each patient
labdate_min_l0 <- data_lab_res$labdate[match(data_l0_res$shcsid,data_lab_res$shcsid)]
labdate_min_l1 <- data_lab_res$labdate[match(data_l1_res$shcsid,data_lab_res$shcsid)]
labdate_max_l0 <- data_lab_res_dcr$labdate[match(data_l0_res$shcsid,data_lab_res_dcr$shcsid)]
labdate_max_l1 <- data_lab_res_dcr$labdate[match(data_l1_res$shcsid,data_lab_res_dcr$shcsid)]
#earliest/latest drug dates for each patient
drugdate_min_l0 <- data_drug_res$drstartd[match(data_l0_res$shcsid,data_drug_res$shcsid)]
drugdate_min_l1 <- data_drug_res$drstartd[match(data_l1_res$shcsid,data_drug_res$shcsid)]
drugdate_max_l0 <- data_drug_res_dcr$drstartd[match(data_l0_res$shcsid,data_drug_res_dcr$shcsid)]
drugdate_max_l1 <- data_drug_res_dcr$drstartd[match(data_l1_res$shcsid,data_drug_res_dcr$shcsid)]
#earliest/latest disease dates for each patient
disdate_min_l0 <- data_dis_res$newdate[match(data_l0_res$shcsid,data_dis_res$shcsid)]
disdate_min_l1 <- data_dis_res$newdate[match(data_l1_res$shcsid,data_dis_res$shcsid)]
disdate_max_l0 <- data_dis_res_dcr$newdate[match(data_l0_res$shcsid,data_dis_res_dcr$shcsid)]
disdate_max_l1 <- data_dis_res_dcr$newdate[match(data_l1_res$shcsid,data_dis_res_dcr$shcsid)]

lower_cutoff_full <- as.Date(paste(lower_cutoff_year,"-01-01",sep=''))
data_l0_res$start_yar_full <- NA
data_l1_res$start_yar_full <- NA
if(start_point=="diag")
{
 # start date : midway point between HIV_POSDATE and HIV_NEGDATE, or just HIV_POSDATE if NEGDATE is missing
 data_l0_res$start_yar_full <- ifelse(!is.na(negdate_l0),round((as.numeric(negdate_l0)+as.numeric(posdate_l0))/2),as.numeric(posdate_l0))
 data_l1_res$start_yar_full <- ifelse(!is.na(negdate_l1),round((as.numeric(negdate_l1)+as.numeric(posdate_l1))/2),as.numeric(posdate_l1))
 # if still missing, start date is earliest from several reported dates, i.e. drug, disease, lab, registration in SHCS
 data_l0_res$start_yar_full <- ifelse(!is.na(data_l0_res$start_yar_full),data_l0_res$start_yar_full,pmin(labdate_min_l0,drugdate_min_l0,disdate_min_l0,data_l0_res$hiv_dreg,na.rm=T))
 data_l1_res$start_yar_full <- ifelse(!is.na(data_l1_res$start_yar_full),data_l1_res$start_yar_full,pmin(labdate_min_l1,drugdate_min_l1,disdate_min_l1,data_l1_res$hiv_dreg,na.rm=T))
} else if(start_point=="reg")
{
 data_l0_res$start_yar_full <- pmin(labdate_min_l0,data_l0_res$hiv_dreg,na.rm=T)
 data_l1_res$start_yar_full <- pmin(labdate_min_l1,data_l1_res$hiv_dreg,na.rm=T)
}

data_l0_res$start_yar_full <- as.Date(data_l0_res$start_yar_full,origin="1970-01-01")
data_l1_res$start_yar_full <- as.Date(data_l1_res$start_yar_full,origin="1970-01-01")

# direct left-censoring
data_l0_res$start_yar_full <- pmax(data_l0_res$start_yar_full,lower_cutoff_full)
data_l1_res$start_yar_full <- pmax(data_l1_res$start_yar_full,lower_cutoff_full)

data_l0_res$start_yar <- as.Year(data_l0_res$start_yar_full)
data_l1_res$start_yar <- as.Year(data_l1_res$start_yar_full)
if(tumor_restriction!="none" & year_restriction)
 data_l1_res <- subset(data_l1_res,start_yar_full<cr_dinc)                                  # removing patients who developed the cancer before start of study

# for the stop date of unlinked patients, taking latest out of all recorded visits (last_info_date variable in admi.dta)
data_l0_res$stop_yar_full <- data_admi$last_info_date[match(data_l0_res$shcsid,data_admi$shcsid)]
data_l0_res$stop_yar <- as.Year(data_l0_res$stop_yar_full)
if(tumor_restriction=="none")
 {data_l1_res$stop_yar_full <- data_admi$last_info_date[match(data_l1_res$shcsid,data_admi$shcsid)]} else
 {data_l1_res$stop_yar_full <- data_l1_res$cr_dinc}                                         # timeline stops at incidence of cancer under study

data_l1_res$stop_yar <- as.Year(data_l1_res$stop_yar_full)

data_l1_res <- data_l1_res[names(data_l0_res)]

data_res_mrg <- rbind(data_l0_res,data_l1_res)
if(tumor_restriction=="none")
 data_res_mrg$linked <- 0

#data_res_mrg <- transform(data_res_mrg,morpho=factor(morpho),topo=factor(topo),cr_nicer=factor(cr_nicer))

data_res_mrg$cr_dinc_y <- as.Year(data_res_mrg$cr_dinc)
 
data_res_mrg$stop_yar <- pmin(data_res_mrg$stop_yar,ifelse(data_res_mrg$registry=="BS",2011,2012))
data_res_mrg <- transform(data_res_mrg,stop_yar_full=pmin(stop_yar_full,ifelse(registry=="BS","2011-12-31","2012-12-31")))
if(year_restriction)
{
 data_res_mrg <- subset(data_res_mrg,!((registry=="BS"&start_yar<1981)|(registry=="GE"&start_yar<1984)|(registry=="ZH"&start_yar<1980)))   # removing patients that appear before start of cancer registries
 data_res_mrg <- subset(data_res_mrg,start_yar_full<=stop_yar_full)                      # removing patients whose start date is later than the stop date
 data_res_mrg$linked[data_res_mrg$cr_dinc>data_res_mrg$stop_yar_full] <- 0                  # unlinking patients whose stop year is before incidence of the tumor
}

data_res_mrg$cr_dinc_y[data_res_mrg$linked==0] <- Inf
data_res_mrg$cr_dinc[data_res_mrg$linked==0] <- NA
data_res_mrg$time_to_event <- pmin(data_res_mrg$cr_dinc_y,data_res_mrg$stop_yar)-data_res_mrg$start_yar

tail_to_mrg <- match(data_res_mrg$shcsid,data_tail$shcsid)

data_res_mrg$sex <- data_tail$sex[tail_to_mrg]
data_res_mrg$sex[data_res_mrg$sex==1] <- "Male"
data_res_mrg$sex[data_res_mrg$sex==2] <- "Female"
data_res_mrg$risk_cat <- data_tail$riskgroup[tail_to_mrg]
data_res_mrg$risk_cat[!data_res_mrg$risk_cat%in%c("IDU","MSM")] <- "other"
data_res_mrg$risk_cat <- factor(data_res_mrg$risk_cat,levels=c("MSM","IDU","other"))

data_res_mrg$hiv_sex[is.na(data_res_mrg$hiv_sex)] <- data_res_mrg$cr_sex[is.na(data_res_mrg$hiv_sex)]       #if gender from SHCS missing, using the one from cancer registry

data_res_mrg$dob_full <- as.Date(paste(data_res_mrg$hiv_dob_y,"-07-01",sep=''))

admi_to_res_mrg <- match(data_res_mrg$shcsid,data_admi$shcsid)
data_res_mrg$region <- data_admi$region[admi_to_res_mrg]
data_res_mrg$region_cat <- NA
data_res_mrg$region_cat[data_res_mrg$region%in%africa_region_id] <- "Africa"
data_res_mrg$region_cat[data_res_mrg$region%in%europe_region_id] <- "Europe"
data_res_mrg$region_cat[data_res_mrg$region%in%other_region_id] <- "Other"
data_res_mrg$region_cat <- factor(data_res_mrg$region_cat)

#start of HAART therapy and B/C clinical status
data_res_mrg$haart_start_date <- data_tail$haart_start_date[tail_to_mrg]
data_res_mrg$first_b_date <- data_tail$first_b_date[tail_to_mrg]
data_res_mrg$first_c_date <- data_tail$first_c_date[tail_to_mrg]  # AIDS incidence

pat_to_mrg <- match(data_res_mrg$shcsid,data_pat$shcsid)

data_res_mrg$ethnicity <- data_pat$ethnicity[pat_to_mrg]
data_res_mrg$ethnicity[data_res_mrg$ethnicity==1] <- "white"
data_res_mrg$ethnicity[data_res_mrg$ethnicity==2] <- "black"
data_res_mrg$ethnicity[data_res_mrg$ethnicity%in%c(0,3,4)] <- "other"
data_res_mrg$ethnicity[data_res_mrg$ethnicity%in%c(9,NA)] <- "unknown"
data_res_mrg$ethnicity <- factor(data_res_mrg$ethnicity,levels=c("white","black","other","unknown"))

data_res_mrg$education <- data_pat$education[pat_to_mrg]
data_res_mrg$education[data_res_mrg$education%in%c(1,2)] <- "class 1"
data_res_mrg$education[data_res_mrg$education%in%c(3,4)] <- "class 2"
data_res_mrg$education[data_res_mrg$education%in%c(5,6,7)] <- "class 3"
data_res_mrg$education[data_res_mrg$education%in%c(0,9,NA)] <- "other/unknown"
data_res_mrg$education <- factor(data_res_mrg$education,levels=c("class 1","class 2","class 3","other/unknown"))

data_res_mrg$ever_smoke <- data_pat$cv_smoked[pat_to_mrg]
data_res_mrg$ever_smoke[data_res_mrg$ever_smoke==9] <- NA

X0 <- subset(data_res_mrg,linked==0)
X1 <- subset(data_res_mrg,linked==1)