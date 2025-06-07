
ptm <- proc.time()
# creating data frame for time-updated Poisson or Cox regression, for a specific cancer
# must run pilot.R first with desired tumor restriction, (IMPORTANT) with lower_cutoff_year at, say, 1980 
# the output can also be used to create the table of charactersitics per cancer
# takes up to 6 minutes with full features
# lab6.R for multiple cancers

library(dplyr)                 # for the 'recode' function; this library overwrites the 'rename' function

start_year <- 1996              # 1988 -> remove RNA viral load from risk factors
stop_year <- 2012
min_visits <- 0                 # minimum required number of follow-up visits to keep a patient in the analysis; 0 for no requirement
max_distance <- 12              # max allowable distance (in months) for risk factor interpolation (6 or 12); "Inf" for no restriction
include_missing <- FALSE        # whether to have separate category for missing values
rf <- c("rna","cd4","cd8")      # which risk factors for which to impose min_visits
write_df <- TRUE                # save the data frame for repeated use?
track_cumul <- FALSE            # whether to determine (time-varying) cumulative durations and nadirs
cat_agg <- FALSE                # whether to aggregate low-pop categories?
Cox_modelling <- TRUE           # whether to add a tstart and tstop variable for each month (takes a couple more minutes)

save_filename <- paste("res_df_12m_nores_",tumor_restriction,"_",start_year,".RData",sep='')
#save_filename <- paste("res_df_Inf_2visits+_",tumor_restriction,"_",start_year,".RData",sep='')
#save_filename <- paste("res_df_12m_nores_",tumor_restriction,"_",start_year,"_linkage_only.RData",sep='')
#save_filename <- paste("res_df_12m_nores_",tumor_restriction,"_",start_year,"_SHCS_only.RData",sep='')
#save_filename <- "res_df_12m_nores_all NADCs_1996.RData"

# name of dataframe is 'res_df'; important: this is created independently of cancer-specific dataset data_res_mrg, e.g. no date restrictions have been made

load(file=paste(filepath_read_R,"\\monthly_rf_df_",max_distance,".RData",sep=''))               # file with monthly interpolations (from 'risk_factors_monthly.R')

res_df <- subset(res_df,age>=16)

res_df <- merge(res_df,data_res_mrg[c("shcsid","start_yar_full","stop_yar_full")],by="shcsid",all.x=TRUE)
res_df$date_y <- as.Year(res_df$date)

# restricting to study period (based on data_res_mrg), note: removes patients who developed cancer on the same month they were registered
#res_df <- subset(res_df,shcsid%in%data_res_mrg$shcsid &  date>=start_yar_full & date<=stop_yar_full)                              
res_df <- subset(res_df,shcsid%in%data_res_mrg$shcsid)

if(track_cumul)
{
 res_df$cd4_nomissing <- res_df$cd4
 res_df$cd4_nomissing[is.na(res_df$cd4_nomissing)] <- Inf
 res_df$rna_nomissing <- res_df$rna
 res_df$rna_nomissing[is.na(res_df$rna_nomissing)] <- 0
 res_df_split <- split(res_df,res_df$shcsid)
 res_df_split <- lapply(res_df_split,transform,haart_cumul=cumsum(haart),cd4_nadir=cummin(cd4_nomissing),rna_cumul=cumsum(rna_nomissing>100000),
                        cd4_cumul200=cumsum(cd4_nomissing<200),cd4_cumul350=cumsum(cd4_nomissing<350),cd4_cumul500=cumsum(cd4_nomissing<500),
                        haart_cumul_ntd=sum(haart),cd4_nadir_ntd=min(cd4_nomissing),rna_cumul_ntd=sum(rna_nomissing>100000),ever_haart=any(haart),
                        cd4_cumul200_ntd=sum(cd4_nomissing<200),cd4_cumul350_ntd=sum(cd4_nomissing<350),cd4_cumul500_ntd=sum(cd4_nomissing<500))
 res_df <- rbind.fill(res_df_split)
 res_df$cd4_nomissing <- NULL
 res_df$rna_nomissing <- NULL
} 
#else
#{
#  res_df <- transform(res_df,haart_cumul=NA,cd4_nadir=NA,rna_cumul=NA,cd4_cumul200=NA,cd4_cumul350=NA,cd4_cumul500=NA,haart_cumul_ntd=NA,cd4_nadir_ntd=NA,
#                      rna_cumul_ntd=NA,ever_haart=NA,cd4_cumul200_ntd=NA,cd4_cumul350_ntd=NA,cd4_cumul500_ntd=NA,ever_hepB=NA,ever_hepC=NA)
#}

res_df <- subset(res_df,as.Year(date)>=start_year & as.Year(date)<=stop_year)                                                      # global time period restrictions

res_df <- subset(res_df,date>=start_yar_full & date<=stop_yar_full)                                           # patient-specific study periods (i.e. stops at cancer incidence)

if(min_visits>0)
{
 res_df_split <- split(res_df,res_df$shcsid)
 if("cd4"%in%rf)
 {
 cd4_counts <- sapply(res_df_split,FUN=function(x) sum(!(x$is_na_cd4)))
 res_df <- subset(res_df,shcsid%in%as.numeric(names(cd4_counts[cd4_counts>=min_visits])))
 }
 if("cd8"%in%rf)
 {
  cd8_counts <- sapply(res_df_split,FUN=function(x) sum(!(x$is_na_cd8)))
  res_df <- subset(res_df,shcsid%in%as.numeric(names(cd8_counts[cd8_counts>=min_visits])))
 }
 if("rna"%in%rf)
 {
  rna_counts <- sapply(res_df_split,FUN=function(x) sum(!(x$is_na_rna)))
  res_df <- subset(res_df,shcsid%in%as.numeric(names(rna_counts[rna_counts>=min_visits])))
 }
 if("hem"%in%rf)
 {
  hem_counts <- sapply(res_df_split,FUN=function(x) sum(!(x$is_na_hem)))
  res_df <- subset(res_df,shcsid%in%as.numeric(names(hem_counts[hem_counts>=min_visits])))
 }
}

res_df <- merge(res_df,data_res_mrg[c("shcsid","registry","linked","risk_cat","sex","education","region_cat","ethnicity","ever_smoke","haart_start_date")],by="shcsid",all.x=TRUE)
res_df <- res_df[order(res_df$shcsid,res_df$date),]
res_df <- plyr::rename(res_df,c("sex"="sex_cat","education"="edu_cat"))

#res_df <- transform(res_df,cd4_cat=factor(cd4_cat,levels=c(levels(cd4_cat),"missing")),cd4_sqrt=sqrt(cd4),cd8_cat=factor(cd8_cat,levels=c(levels(cd8_cat),"missing")),
#                    rna_cat=factor(rna_cat,levels=c(levels(rna_cat),"missing")),cd4Rcd8_cat=factor(cd4Rcd8_cat,levels=c(levels(cd4Rcd8_cat),"missing")),
#                    weight_cat=factor(weight_cat,levels=c(levels(weight_cat),"missing")))

res_df$risk_cat <- as.character(res_df$risk_cat)
res_df$risk_cat <- with(res_df,ifelse(risk_cat!="other",risk_cat,ifelse(sex_cat=="Male","Other men","Other women")))
if(tumor_restriction=="liver" && cat_agg)
{
 res_df$risk_cat[res_df$risk_cat%in%c("Other men","Other women")] <- "Other men/women"
 res_df$risk_cat <- factor(res_df$risk_cat,levels=c("MSM","IDU","Other men/women"))
} else
{res_df$risk_cat <- factor(res_df$risk_cat,levels=c("MSM","IDU","Other men","Other women"))}
if(tumor_restriction=="prostate")
 res_df <- subset(res_df,sex_cat=="Male" & risk_cat!="Other women")
if(tumor_restriction=="cervix")
 res_df <- subset(res_df,sex_cat=="Female" & risk_cat%in%c("IDU","Other women"))

res_df$region_cat <- factor(res_df$region_cat,levels=c("Europe","Africa","Other"))
res_df$ever_smoke <- factor(dplyr::recode(res_df$ever_smoke,'1'="yes",'0'="no"),levels=c("no","yes"))

# proper category for missing values
if(include_missing)
 {
  res_df$cd4_cat[is.na(res_df$cd4_cat)] <- "missing"
  res_df$cd8_cat[is.na(res_df$cd8_cat)] <- "missing"
  res_df$rna_cat[is.na(res_df$rna_cat)] <- "missing"
  res_df$cd4Rcd8_cat[is.na(res_df$cd4Rcd8_cat)] <- "missing"
  res_df$weight_cat[is.na(res_df$weight_cat)] <- "missing"
  #res_df$ever_smoke[is.na(res_df$ever_smoke)] <- "missing"
 }

status_at_reg <- res_df[!duplicated(res_df$shcsid),]                  # first visits
reg_to_df_res <- match(res_df$shcsid,status_at_reg$shcsid)
res_df$haart_at_reg <- status_at_reg$haart[reg_to_df_res]
res_df$age_at_reg <- status_at_reg$age[reg_to_df_res]
res_df$age_at_reg_cat <- status_at_reg$age_cat[reg_to_df_res]
res_df$cd4_at_reg <- status_at_reg$cd4[reg_to_df_res]
res_df$cd4_at_reg_cat <- status_at_reg$cd4_cat[reg_to_df_res]
res_df$cd8_at_reg <- status_at_reg$cd8[reg_to_df_res]
res_df$cd8_at_reg_cat <- status_at_reg$cd8_cat[reg_to_df_res]
res_df$rna_at_reg <- status_at_reg$rna[reg_to_df_res]
res_df$rna_at_reg_cat <- status_at_reg$rna_cat[reg_to_df_res]
res_df$hem_at_reg <- status_at_reg$hem[reg_to_df_res]
res_df$date_at_reg <- status_at_reg$date[reg_to_df_res]
res_df$date_at_reg_cat <- status_at_reg$date_cat[reg_to_df_res]
res_df$hepB_pos_at_reg <- status_at_reg$hepB_pos[reg_to_df_res]
res_df$hepC_aHCVpos_at_reg <- status_at_reg$hepC_aHCVpos[reg_to_df_res]
res_df$hepC_RNApos_at_reg <- status_at_reg$hepC_RNApos[reg_to_df_res]
res_df$hepC_RNAdet_at_reg <- status_at_reg$hepC_RNAdet[reg_to_df_res]
res_df$hepC_RNA_at_reg <- status_at_reg$hepC_RNA[reg_to_df_res]
res_df$CDC_at_reg_cat <- status_at_reg$CDC_cat[reg_to_df_res]

if(Cox_modelling)           # adding start and stop time indicators for each month
{
  res_df_split <- split(res_df,res_df$shcsid)
  X <- res_df_split
  for(i in 1:length(X))
  {
    Y <- X[[i]]
    Y$tstart <- 0:(nrow(Y)-1)
    Y$tstop <- Y$tstart+1
    Y$linked[1:(nrow(Y)-1)] <- 0
    X[[i]] <- Y
  }
  res_df <- rbind.fill(X)
  rm(X,Y)
}

print(paste("Time period: ",start_year,"-",stop_year,sep=''))
print(paste("Max interpolation distance:",max_distance))
print(paste("Minimum required number of visits:",min_visits))

if(write_df)
{
 save(res_df,file=paste(filepath_read_R,"\\",save_filename,sep=''))
 print(paste("saving res_df for",tumor_restriction))
 print(paste("filename: ",save_filename))
 rm(res_df)
} else
{print("res_df NOT saved!")}

print(proc.time()-ptm)