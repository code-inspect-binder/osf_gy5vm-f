# creating monthly values for risk factors though interpolation, with 'barriers' at cancer incidence and start of haart
# run pilot.R first (with any date/tumor_restriction) to load in the various SHCS datasets
# this code does NOT make use of data_res_mrg
# good examples: shcsid=31256; 10610; 41088
# 31256 - check around breakpoint at 2001-09
# 10015 - no breakpoint

library(plyr)                  # for fast reassembling of data frames

ptm <- proc.time()             # takes over 6 minutes

start_date <- as.Date("1981-01-01")
stop_date <- as.Date("2014-01-01")
max_dist <- 12                                                     # max allowable distance (in months) for interpolation; "Inf" for no restriction
write_df <- TRUE                                                   # save the data frame for repeated use?

data_pat_ext <- subset(data_pat,shcsid%in%data_l$shcsid)                         #restricting to patients for which linkage was attempted
data_pat_ext$dob_full <- as.Date(paste(data_pat_ext$born,"-07-01",sep=''))
data_pat_ext$haart_start_date <- data_tail$haart_start_date[match(data_pat_ext$shcsid,data_tail$shcsid)]
data_pat_ext$first_b_date <- data_tail$first_b_date[match(data_pat_ext$shcsid,data_tail$shcsid)]
data_pat_ext$first_c_date <- data_tail$first_c_date[match(data_pat_ext$shcsid,data_tail$shcsid)]

data_lab_cd4 <- data_lab
data_lab_cd4 <- transform(data_lab_cd4,cd4date=as.Date(ifelse(!is.na(cd4date),cd4date,labdate),origin="1970-01-01"))         
data_lab_cd4 <- subset(data_lab_cd4,!is.na(cd4),select=c(shcsid,cd4,cd4date))
data_lab_cd4$visit_id <- paste(data_lab_cd4$shcsid,substr(data_lab_cd4$cd4date,start=1,stop=7))
data_lab_cd4 <- data_lab_cd4[order(data_lab_cd4$shcsid,data_lab_cd4$cd4date),]
data_lab_cd4 <- subset(data_lab_cd4,!duplicated(data_lab_cd4$visit_id,fromLast=TRUE))  # when multiple measurements in one months, keeping latest one
data_lab_cd4 <- subset(data_lab_cd4,cd4date>=start_date & cd4date<stop_date)

data_lab_cd8 <- data_lab
data_lab_cd8 <- transform(data_lab_cd8,cd8date=as.Date(ifelse(!is.na(cd4date),cd4date,labdate),origin="1970-01-01"))         
data_lab_cd8 <- subset(data_lab_cd8,!is.na(cd8),select=c(shcsid,cd8,cd8date))
data_lab_cd8$visit_id <- paste(data_lab_cd8$shcsid,substr(data_lab_cd8$cd8date,start=1,stop=7))
data_lab_cd8 <- data_lab_cd8[order(data_lab_cd8$shcsid,data_lab_cd8$cd8date),]
data_lab_cd8 <- subset(data_lab_cd8,!duplicated(data_lab_cd8$visit_id,fromLast=TRUE))
data_lab_cd8 <- subset(data_lab_cd8,cd8date>=start_date & cd8date<=stop_date)

data_lab_rna <- data_lab
data_lab_rna <- transform(data_lab_rna,rnadate=labdate)         
data_lab_rna <- subset(data_lab_rna,!is.na(rna),select=c(shcsid,rna,rnadate))
data_lab_rna$visit_id <- paste(data_lab_rna$shcsid,substr(data_lab_rna$rnadate,start=1,stop=7))
data_lab_rna <- data_lab_rna[order(data_lab_rna$shcsid,data_lab_rna$rnadate),]
data_lab_rna <- subset(data_lab_rna,!duplicated(data_lab_rna$visit_id,fromLast=TRUE))
data_lab_rna <- subset(data_lab_rna,rnadate>=start_date & rnadate<=stop_date)

data_lab_hem <- data_lab
data_lab_hem <- transform(data_lab_hem,hemdate=labdate)
data_lab_hem <- subset(data_lab_hem,!is.na(hem),select=c(shcsid,hem,hemdate))
data_lab_hem$visit_id <- paste(data_lab_hem$shcsid,substr(data_lab_hem$hemdate,start=1,stop=7))
data_lab_hem <- data_lab_hem[order(data_lab_hem$shcsid,data_lab_hem$hemdate),]
data_lab_hem <- subset(data_lab_hem,!duplicated(data_lab_hem$visit_id,fromLast=TRUE))
data_lab_hem <- subset(data_lab_hem,hemdate>=start_date & hemdate<=stop_date)

data_lab_weight <- data_lab
data_lab_weight <- transform(data_lab_weight,weightdate=labdate)         
data_lab_weight <- subset(data_lab_weight,!is.na(weight),select=c(shcsid,weight,weightdate))
data_lab_weight$visit_id <- paste(data_lab_weight$shcsid,substr(data_lab_weight$weightdate,start=1,stop=7))
data_lab_weight <- data_lab_weight[order(data_lab_weight$shcsid,data_lab_weight$weightdate),]
data_lab_weight <- subset(data_lab_weight,!duplicated(data_lab_weight$visit_id,fromLast=TRUE))
data_lab_weight <- subset(data_lab_weight,weightdate>=start_date & weightdate<=stop_date)

# data_lab_hepC1_test <- subset(data_lab,!is.na(antihcv))
# data_lab_hepC1_test <- transform(data_lab_hepC1_test,hepCdate=hec_date)
# data_lab_hepC2_test <- subset(data_lab,!is.na(hcv_rna_qual))
# data_lab_hepC2_test <- transform(data_lab_hepC2_test,hepCdate=hcv_rna_date)
# data_lab_hepC3_test <- subset(data_lab,!is.na(hcv_rna))
# data_lab_hepC3_test <- transform(data_lab_hepC3_test,hepCdate=hcv_rna_date)
# data_lab_hepC_test <- rbind(data_lab_hepC1_test,data_lab_hepC2_test,data_lab_hepC3_test)
# rm(data_lab_hepC1_test,data_lab_hepC2_test,data_lab_hepC3_test)
# data_lab_hepC_test <- transform(data_lab_hepC_test,hepCdate=ifelse(is.na(hepCdate),labdate,hepCdate))
# data_lab_hepC_test <- transform(data_lab_hepC_test,hepCdate=as.Date(hepCdate,origin="1970-01-01"))
# data_lab_hepC_test <- data_lab_hepC_test[order(data_lab_hepC_test$shcsid,data_lab_hepC_test$hepCdate),]
# data_lab_hepC_test <- subset(data_lab_hepC_test,hepCdate>=start_date & hepCdate<=stop_date)
data_lab_hepC_aHCVpos <- subset(data_lab,antihcv=="P")
data_lab_hepC_aHCVpos <- rename(data_lab_hepC_aHCVpos,c("hec_date"="hepC_aHCVpos_date"))
data_lab_hepC_aHCVpos <- transform(data_lab_hepC_aHCVpos,hepC_aHCVpos_date=as.Date(ifelse(is.na(hepC_aHCVpos_date),labdate,hepC_aHCVpos_date),origin="1970-01-01"))
data_lab_hepC_aHCVpos <- data_lab_hepC_aHCVpos[order(data_lab_hepC_aHCVpos$shcsid, data_lab_hepC_aHCVpos$hepC_aHCVpos_date),]
data_lab_hepC_aHCVpos <- subset(data_lab_hepC_aHCVpos,!duplicated(shcsid))

data_lab_hepC_RNApos <- subset(data_lab,hcv_rna_qual=="P")
data_lab_hepC_RNApos <- rename(data_lab_hepC_RNApos,c("hcv_rna_date"="hepC_RNApos_date"))
data_lab_hepC_RNApos <- transform(data_lab_hepC_RNApos,hepC_RNApos_date=as.Date(ifelse(is.na(hepC_RNApos_date),labdate,hepC_RNApos_date),"1970-01-01"))
data_lab_hepC_RNApos <- data_lab_hepC_RNApos[order(data_lab_hepC_RNApos$shcsid, data_lab_hepC_RNApos$hepC_RNApos_date),]
data_lab_hepC_RNApos <- subset(data_lab_hepC_RNApos,!duplicated(shcsid))

data_lab_hepC_RNAdet <- subset(data_lab,hcv_rna>0)
data_lab_hepC_RNAdet <- rename(data_lab_hepC_RNAdet,c("hcv_rna_date"="hepC_RNAdet_date"))
data_lab_hepC_RNAdet <- transform(data_lab_hepC_RNAdet,hepC_RNAdet_date=as.Date(ifelse(is.na(hepC_RNAdet_date),labdate,hepC_RNAdet_date),"1970-01-01"))
data_lab_hepC_RNAdet <- data_lab_hepC_RNAdet[order(data_lab_hepC_RNAdet$shcsid, data_lab_hepC_RNAdet$hepC_RNAdet_date),]
data_lab_hepC_RNAdet <- subset(data_lab_hepC_RNAdet,!duplicated(shcsid))

data_pat_ext <- merge(data_pat_ext,data_lab_hepC_aHCVpos[c("shcsid","hepC_aHCVpos_date")],by="shcsid",all.x=TRUE)
data_pat_ext <- merge(data_pat_ext,data_lab_hepC_RNApos[c("shcsid","hepC_RNApos_date")],by="shcsid",all.x=TRUE)
data_pat_ext <- merge(data_pat_ext,data_lab_hepC_RNAdet[c("shcsid","hepC_RNAdet_date")],by="shcsid",all.x=TRUE)

# data_pat_ext$hepCtest_date <- NA
# data_pat_ext$hepCtest_date <- data_lab_hepC_test$hepCdate[match(data_pat_ext$shcsid,data_lab_hepC_test$shcsid)]
# data_pat_ext$hepCaHCV_date <- NA
# data_pat_ext$hepCpos_date <- data_lab_hepC_pos$hepCdate[match(data_pat_ext$shcsid,data_lab_hepC_pos$shcsid)]
# data_pat_ext$hepCRNApos_date <- NA
# data_pat_ext$hepCRNApos_date <- data_lab_hepC_RNApos$hepCdate[match(data_pat_ext$shcsid,data_lab_hepC_RNApos$shcsid)]

data_lab_hepB1_test <- subset(data_lab,!is.na(aghbs))
data_lab_hepB1_test <- transform(data_lab_hepB1_test,hepBdate=heb_date)
data_lab_hepB2_test <- subset(data_lab,!is.na(hbvdna_qual))
data_lab_hepB2_test <- transform(data_lab_hepB2_test,hepBdate=hbvdna_date)
data_lab_hepB_test <- rbind(data_lab_hepB1_test,data_lab_hepB2_test)
rm(data_lab_hepB1_test,data_lab_hepB2_test)
data_lab_hepB_test <- transform(data_lab_hepB_test,hepBdate=ifelse(is.na(hepBdate),labdate,hepBdate))
data_lab_hepB_test <- transform(data_lab_hepB_test,hepBdate=as.Date(hepBdate,origin="1970-01-01"))
data_lab_hepB_test <- data_lab_hepB_test[order(data_lab_hepB_test$shcsid,data_lab_hepB_test$hepBdate),]
data_lab_hepB_test <- subset(data_lab_hepB_test,hepBdate>=start_date & hepBdate<=stop_date)
data_lab_hepB_pos <- subset(data_lab_hepB_test,aghbs=="P" | hbvdna_qual=="P")
data_lab_hepB_test <- subset(data_lab_hepB_test,!duplicated(shcsid))
data_lab_hepB_pos <- subset(data_lab_hepB_pos,!duplicated(shcsid))

data_pat_ext$hepBtest_date <- NA
data_pat_ext$hepBtest_date <- data_lab_hepB_test$hepBdate[match(data_pat_ext$shcsid,data_lab_hepB_test$shcsid)]
data_pat_ext$hepBpos_date <- NA
data_pat_ext$hepBpos_date <- data_lab_hepB_pos$hepBdate[match(data_pat_ext$shcsid,data_lab_hepB_pos$shcsid)]

data_dis_bpn <- subset(data_dis,disease=="BPN")
names(data_dis_bpn)[3] <- "bpn_date"
data_pat_ext$bpn_date <- NA
data_pat_ext$bpn_date <- data_dis_bpn$bpn_date[match(data_pat_ext$shcsid,data_dis_bpn$shcsid)]

rm(data_lab_hepB_test,data_lab_hepB_pos)

data_lab_cd4_split <- split(data_lab_cd4,data_lab_cd4$shcsid)
data_lab_cd8_split <- split(data_lab_cd8,data_lab_cd8$shcsid)
data_lab_rna_split <- split(data_lab_rna,data_lab_rna$shcsid)
data_lab_hem_split <- split(data_lab_hem,data_lab_hem$shcsid)
data_lab_weight_split <- split(data_lab_weight,data_lab_weight$shcsid)

if(write_df)
  save(data_pat_ext,file=paste(filepath_read_R,"\\data_pat_ext.RData",sep=''))

monthly_rf <- function(i,max_dist=6,start,stop,barriers=TRUE)
{
 #print(i)
 pat_record <- data_pat_ext[i,]
 pat_id <- pat_record$shcsid
 lab_record_cd4 <- data_lab_cd4_split[[as.character(pat_id)]]
 lab_record_cd8 <- data_lab_cd8_split[[as.character(pat_id)]]
 lab_record_rna <- data_lab_rna_split[[as.character(pat_id)]]
 lab_record_hem <- data_lab_hem_split[[as.character(pat_id)]]
 lab_record_weight <- data_lab_weight_split[[as.character(pat_id)]]
 bp <- breakpoints[[as.character(pat_id)]]
 df <- data.frame(shcsid=rep(pat_id,month_delta(stop,start)+1))
 
 if(!is.null(lab_record_cd4))
 {
  lab_record_cd4 <- transform(lab_record_cd4,month_d=month_delta(lab_record_cd4$cd4date,start))
  df$cd4 <- NA
  df$cd4[lab_record_cd4$month_d+1] <- lab_record_cd4$cd4                                        # CD4 value updated at beginning of following month
 } else
 {
  df$cd4 <- NA
 }
 if(!is.null(lab_record_cd8))
 {
  lab_record_cd8 <- transform(lab_record_cd8,month_d=month_delta(lab_record_cd8$cd8date,start))
  df$cd8 <- NA
  df$cd8[lab_record_cd8$month_d+1] <- lab_record_cd8$cd8
 } else
 {
  df$cd8 <- NA
 }
 if(!is.null(lab_record_rna))
 {
  lab_record_rna <- transform(lab_record_rna,month_d=month_delta(lab_record_rna$rnadate,start))                      
  df$rna <- NA
  df$rna[lab_record_rna$month_d+1] <- lab_record_rna$rna
 } else
 {
  df$rna <- NA
 }
 if(!is.null(lab_record_hem))
 {
   lab_record_hem <- transform(lab_record_hem,month_d=month_delta(lab_record_hem$hemdate,start))                      
   df$hem <- NA
   df$hem[lab_record_hem$month_d+1] <- lab_record_hem$hem
 } else
 {
   df$hem <- NA
 }
 if(!is.null(lab_record_weight))
 {
  lab_record_weight <- transform(lab_record_weight,month_d=month_delta(lab_record_weight$weightdate,start))                                   
  df$weight <- NA
  df$weight[lab_record_weight$month_d+1] <- lab_record_weight$weight
 } else
 {
  df$weight <- NA
 }
 
 month_vect <- ((as.Month(start)+1):(as.Month(start)+nrow(df)))%%12
 month_vect[month_vect==0] <- 12
 year_vect <- as.Year(start) + ((as.Month(start)+1):(as.Month(start)+nrow(df)))%/%12
 year_vect <- ifelse(month_vect==12,year_vect-1,year_vect)

 df$date <- as.Date(paste(year_vect,"-",month_vect,"-01",sep=''))
 df$haart <- FALSE
 df$CDC_cat <- "A"
 if(!is.na(pat_record$haart_start_date))
  df$haart[df$date>pat_record$haart_start_date] <- TRUE
 if(!is.na(pat_record$first_b_date))
  df$CDC_cat[df$date>pat_record$first_b_date] <- "B"
 if(!is.na(pat_record$first_c_date))
  df$CDC_cat[df$date>pat_record$first_c_date] <- "C"
 df$CDC_cat <- factor(df$CDC_cat,levels=c("A","B","C"))
 df$hepC_aHCVpos <- FALSE
 if(!is.na(pat_record$hepC_aHCVpos_date))
  df$hepC_aHCVpos[df$date>pat_record$hepC_aHCVpos_date] <- TRUE
 df$hepC_RNApos <- FALSE
 if(!is.na(pat_record$hepC_RNApos_date))
   df$hepC_RNApos[df$date>pat_record$hepC_RNApos_date] <- TRUE
 df$hepC_RNAdet <- FALSE
 if(!is.na(pat_record$hepC_RNAdet_date))
   df$hepC_RNAdet[df$date>pat_record$hepC_RNAdet_date] <- TRUE
 df$hepB_test <- FALSE
 if(!is.na(pat_record$hepBtest_date))
  df$hepB_test[df$date>pat_record$hepBtest_date] <- TRUE
 df$hepB_pos <- FALSE
 if(!is.na(pat_record$hepBpos_date))
  df$hepB_pos[df$date>pat_record$hepBpos_date] <- TRUE
 df$bpn <- FALSE
 if(!is.na(pat_record$bpn_date))
  df$bpn[df$date>pat_record$bpn_date] <- TRUE
 
 df$is_na_cd4 <- is.na(df$cd4)
 df$is_na_cd8 <- is.na(df$cd8)
 df$is_na_rna <- is.na(df$rna)
 df$is_na_hem <- is.na(df$hem)
 df$is_na_weight <- is.na(df$weight)

 # filling in the holes (nearest neighbour - in case of ties taking the past value rather than the future one)
 ind_not_na_cd4 <- which(!is.na(df$cd4))
 ind_not_na_cd8 <- which(!is.na(df$cd8))
 ind_not_na_rna <- which(!is.na(df$rna))
 ind_not_na_hem <- which(!is.na(df$hem))
 ind_not_na_weight <- which(!is.na(df$weight))
 if(length(ind_not_na_cd4)>0)
 {
  if(is.null(bp) || barriers==FALSE)                                                                 # no thresholds
  {
   M <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_cd4))-rep(1,nrow(df))%o%ind_not_na_cd4)     # for each month (row), "distance" to visits with recorded CD4 (columns) 
  } else
  {
   M0 <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_cd4))-rep(1,nrow(df))%o%ind_not_na_cd4)
   M <- matrix(Inf,nrow=nrow(df),ncol=length(ind_not_na_cd4))
   bp_temp <- c(as.Date("1900-01-01"),bp,as.Date("2500-01-01"))
   for(j in 1:(length(bp_temp)-1))                                                           # taking thresholds (cancer incidence, start of haart) into account
   {
    thresh <- c(bp_temp[j],bp_temp[j+1])
    A <- df$date>thresh[1] & df$date<=thresh[2]
    B <- df$date[ind_not_na_cd4]>thresh[1] & df$date[ind_not_na_cd4]<=thresh[2]
    M[A,B] <- M0[A,B]
   }
  }
  min_ind <- ind_not_na_cd4[apply(M,1,which.min)]                                                                    # for each month, index of nearest cd4 count
  df$cd4 <- df$cd4[min_ind]
  na_ind <- apply(M>=(max_dist+1),1,all)                                                                                  # months with no close enough visit
  df$cd4[na_ind] <- NA
 }
 if(length(ind_not_na_cd8)>0)
 {
  if(is.null(bp) || barriers==FALSE)
  {
   M <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_cd8))-rep(1,nrow(df))%o%ind_not_na_cd8)
  } else
  {
   M0 <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_cd8))-rep(1,nrow(df))%o%ind_not_na_cd8)
   M <- matrix(Inf,nrow=nrow(df),ncol=length(ind_not_na_cd8))
   bp_temp <- c(as.Date("1900-01-01"),bp,as.Date("2500-01-01"))
   for(j in 1:(length(bp_temp)-1))
   {
    thresh <- c(bp_temp[j],bp_temp[j+1])
    A <- df$date>thresh[1] & df$date<=thresh[2]
    B <- df$date[ind_not_na_cd8]>thresh[1] & df$date[ind_not_na_cd8]<=thresh[2]
    M[A,B] <- M0[A,B]
   }
  }
  min_ind <- ind_not_na_cd8[apply(M,1,which.min)]
  df$cd8 <- df$cd8[min_ind]
  na_ind <- apply(M>=(max_dist+1),1,all)
  df$cd8[na_ind] <- NA
 }
 if(length(ind_not_na_rna)>0)
 {
  if(is.null(bp)|| barriers==FALSE)
  {
   M <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_rna))-rep(1,nrow(df))%o%ind_not_na_rna)
  } else
  {
   M0 <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_rna))-rep(1,nrow(df))%o%ind_not_na_rna)
   M <- matrix(Inf,nrow=nrow(df),ncol=length(ind_not_na_rna))
   bp_temp <- c(as.Date("1900-01-01"),bp,as.Date("2500-01-01"))
   for(j in 1:(length(bp_temp)-1))
   {
    thresh <- c(bp_temp[j],bp_temp[j+1])
    A <- df$date>thresh[1] & df$date<=thresh[2]
    B <- df$date[ind_not_na_rna]>thresh[1] & df$date[ind_not_na_rna]<=thresh[2]
    M[A,B] <- M0[A,B]
   }
  }
  min_ind <- ind_not_na_rna[apply(M,1,which.min)]
  df$rna <- df$rna[min_ind]
  na_ind <- apply(M>=(max_dist+1),1,all)
  df$rna[na_ind] <- NA
 }
 if(length(ind_not_na_hem)>0)
 {
   if(is.null(bp)|| barriers==FALSE)
   {
     M <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_hem))-rep(1,nrow(df))%o%ind_not_na_hem)
   } else
   {
     M0 <- abs(1:nrow(df)%o%rep(1,length(ind_not_na_hem))-rep(1,nrow(df))%o%ind_not_na_hem)
     M <- matrix(Inf,nrow=nrow(df),ncol=length(ind_not_na_hem))
     bp_temp <- c(as.Date("1900-01-01"),bp,as.Date("2500-01-01"))
     for(j in 1:(length(bp_temp)-1))
     {
       thresh <- c(bp_temp[j],bp_temp[j+1])
       A <- df$date>thresh[1] & df$date<=thresh[2]
       B <- df$date[ind_not_na_hem]>thresh[1] & df$date[ind_not_na_hem]<=thresh[2]
       M[A,B] <- M0[A,B]
     }
   }
   min_ind <- ind_not_na_hem[apply(M,1,which.min)]
   df$hem <- df$hem[min_ind]
   na_ind <- apply(M>=(max_dist+1),1,all)
   df$hem[na_ind] <- NA
 }
 if(length(ind_not_na_weight)>0)
 {
  # weight: no thresholds
  M <- abs(cbind(1:nrow(df))%*%rep(1,length(ind_not_na_weight))-cbind(rep(1,nrow(df)))%*%rbind(ind_not_na_weight))
  min_ind <- ind_not_na_weight[apply(M,1,which.min)]
  df$weight <- df$weight[min_ind]
 }
 df$age <- month_delta(df$date,pat_record$dob_full)%/%12
 df
}

res <- lapply(1:nrow(data_pat_ext),FUN=monthly_rf,max_dist=max_dist,start=start_date,stop=stop_date)
names(res) <- factor(data_pat_ext$shcsid,levels=data_pat_ext$shcsid)

# creating variables for lagged values and moving averages

delay_rf <- function(v,lag)
  c(rep(NA,lag),v[1:(length(v)-lag)])

moving_avg <- function(x,n)
  c(stats::filter(x,rep(1/n,n),sides=1))

res <- lapply(res,transform,cd4Rcd8=cd4/cd8)
res <- lapply(res,transform,cd4_ma12=moving_avg(cd4,12),cd8_ma12=moving_avg(cd8,12),cd4Rcd8_ma12=moving_avg(cd4Rcd8,12),cd4_ma24=moving_avg(cd4,24),cd8_ma24=moving_avg(cd8,24),cd4Rcd8_ma24=moving_avg(cd4Rcd8,24),
                            rna_ma12=moving_avg(rna,12),rna_ma24=moving_avg(rna,24))
res <- lapply(res,transform,cd4_lag12=delay_rf(cd4,12),cd4_lag24=delay_rf(cd4,24),cd4_lag36=delay_rf(cd4,36),
                            cd8_lag12=delay_rf(cd8,12),cd8_lag24=delay_rf(cd8,24),cd8_lag36=delay_rf(cd8,36),
                            cd4Rcd8_lag12=delay_rf(cd4Rcd8,12),cd4Rcd8_lag24=delay_rf(cd4Rcd8,24),cd4Rcd8_lag36=delay_rf(cd4Rcd8,36),
                            rna_lag12=delay_rf(rna,12),rna_lag24=delay_rf(rna,24),rna_lag36=delay_rf(rna,36),
                            cd4_ma12_lag12=delay_rf(cd4_ma12,12),cd4_ma24_lag12=delay_rf(cd4_ma24,12),cd4_ma12_lag24=delay_rf(cd4_ma12,24),
                            cd8_ma12_lag12=delay_rf(cd8_ma12,12),cd8_ma24_lag12=delay_rf(cd8_ma24,12),cd8_ma12_lag24=delay_rf(cd8_ma12,24),
                            cd4Rcd8_ma12_lag12=delay_rf(cd4Rcd8_ma12,12),cd4Rcd8_ma24_lag12=delay_rf(cd4Rcd8_ma24,12),cd4Rcd8_ma12_lag24=delay_rf(cd4Rcd8_ma12,24),
                            rna_ma12_lag12=delay_rf(rna_ma12,12),rna_ma24_lag12=delay_rf(rna_ma24,12),rna_ma12_lag24=delay_rf(rna_ma12,24))
res <- lapply(res,transform,haart_lag6=delay_rf(haart,6),haart_lag12=delay_rf(haart,12),hepC_aHCVpos_lag12=delay_rf(hepC_aHCVpos,12),hepC_RNApos_lag12=delay_rf(hepC_RNApos,12),
                            hepC_RNAdet_lag12=delay_rf(hepC_RNAdet,12),hepB_pos_lag12=delay_rf(hepB_pos,12),CDC_lag12_cat=delay_rf(CDC_cat,12),bpn_lag12=delay_rf(bpn,12))



res_df <- rbind.fill(res)
res_df$CDC_lag12_cat[res_df$CDC_lag12_cat==1] <- "A"
res_df$CDC_lag12_cat[res_df$CDC_lag12_cat==2] <- "B"
res_df$CDC_lag12_cat[res_df$CDC_lag12_cat==3] <- "C"
res_df$CDC_lag12_cat <- factor(res_df$CDC_lag12_cat,levels=c("A","B","C"))
res_df$hepC_aHCVpos_lag12[is.na(res_df$hepC_aHCVpos_lag12)] <- FALSE
res_df$hepC_RNApos_lag12[is.na(res_df$hepC_RNApos_lag12)] <- FALSE
res_df$hepC_RNAdet_lag12[is.na(res_df$hepC_RNAdet_lag12)] <- FALSE
res_df <- transform(res_df,hepC_RNA=(hepC_RNApos | hepC_RNAdet),hepC_RNA_lag12=(hepC_RNApos_lag12 | hepC_RNAdet_lag12))
res_df$bpn_lag12[is.na(res_df$bpn_lag12)] <- FALSE

rownames(res_df) <- NULL

if(write_df)
 save(res_df,file=paste(filepath_read_R,"\\monthly_rf_df_",max_dist,".RData",sep=''))

print(proc.time()-ptm)