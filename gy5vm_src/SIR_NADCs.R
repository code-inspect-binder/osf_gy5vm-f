# Standardized incidence ratios for NADCs, adjusted for age, sex, and (optional) time period
# must first run pilot.R with desired tumor restriction and lower cutoff year (currently only works with 1988 and 1996)
# 'other NADCs' must be specificed manually below for data_cr
# Swiss cancer rates obtained from NICER: www.nicer.org
# forestplots -> lab12.R

library(ggplot2)
library(car)

which_gender <- "both"                   # "males", "females", or "both"
if(tumor_restriction=="prostate")
  which_gender <- "males"

remove_ties <- FALSE
time_adj <- TRUE               # whether to adjust for time period

print(paste("tumor:",tumor_restriction))
print(paste("gender:",which_gender))

# general cancer rates in Switzerland from 1989 to 2013
data_cr <- read.csv(file=paste(filepath_read_NICER,"\\cr_CH_1year.csv",sep=''),sep=";")

tumor_restriction_nicer <- switch(tumor_restriction,
                            "anus"="Anus & Anal Canal",
                            "lung"="Lung, Bronchus, Trachea",
                            "head and neck"=c("Oral Cavity & Pharynx","Larynx"),
                            "prostate"="Prostate",
                            "liver"="Liver & Intrahepatic Bile Ducts",
                            "HL"="Hodgkin Lymphoma",
                            "NHL"="Non Hodgkin Lymphoma")

data_cr <- subset(data_cr,region=="Switzerland" & calendar_year>=lower_cutoff_year & calendar_year<=2012)
if(tumor_restriction!="other NADCs")
{data_cr <- subset(data_cr,cancer%in%tumor_restriction_nicer)} else
{
  #data_cr <- subset(data_cr,!cancer%in%c("Non Hodgkin Lymphoma","Cervix Uteri","All Cancers (See Methods for inclusion criteria)","Anus & Anal Canal",
  #                                      "Lung, Bronchus, Trachea","Oral Cavity & Pharynx","Larynx","Prostate","Liver & Intrahepatic Bile Ducts","Hodgkin Lymphoma"))
  data_cr <- subset(data_cr,!cancer%in%c("Non Hodgkin Lymphoma","Cervix Uteri","All Cancers (See Methods for inclusion criteria)","Anus & Anal Canal",
                                         "Lung, Bronchus, Trachea","Prostate","Liver & Intrahepatic Bile Ducts","Hodgkin Lymphoma"))
}
data_cr$sex <- factor(data_cr$sex,levels=c("males","females"))

if(which_gender!="both")
  data_cr <- subset(data_cr,sex==which_gender)

# carrying over 1989 rates to 1988
if(lower_cutoff_year==1988)
{
 X <- subset(data_cr,calendar_year==1989)
 X$calendar_year <- 1988
 data_cr <- rbind(X,data_cr)
 rm(X)
}

#if(lower_cutoff_year==1988)
#  period_breaks <- c(1988,1996,2002,2008,2013)
#if(lower_cutoff_year==1996)
#  period_breaks <- c(1996,2002,2008,2013)

data_cr$time_period <- NA
if(lower_cutoff_year==1988)
  data_cr$time_period[data_cr$calendar_year%in%1988:1995] <- "[1988,1996)"
data_cr$time_period[data_cr$calendar_year%in%1996:2001] <- "[1996,2002)"
data_cr$time_period[data_cr$calendar_year%in%2002:2007] <- "[2002,2008)"
data_cr$time_period[data_cr$calendar_year%in%2008:2012] <- "[2008,2013)"

data_cr$age_grp <- NA
data_cr$age_grp[data_cr$age_cat%in%c("0to4","5to9","10to14","15to19","20to24","25to29","30to34")] <- "[0,35)"
data_cr$age_grp[data_cr$age_cat%in%c("35to39","40to44")] <- "[35,45)"
data_cr$age_grp[data_cr$age_cat%in%c("45to49","50to54")] <- "[45,55)"
data_cr$age_grp[data_cr$age_cat%in%c("55to59","60to64","65to69","70to74","75to79","80to84","85+")] <- "[55,Inf)"

data_cr <- transform(data_cr,time_period=factor(time_period),age_grp=factor(data_cr$age_grp))

if(time_adj)
{data_cr_agg <- aggregate(cbind(obs_cases,pop)~ age_grp + time_period + sex + cancer,data=data_cr,FUN=sum)} else
{data_cr_agg <- aggregate(cbind(obs_cases,pop)~ age_grp + sex + cancer,data=data_cr,FUN=sum)}
data_cr_agg <- transform(data_cr_agg,rate=obs_cases/pop)

data_cr_agg <- aggregate(rate ~ age_grp +time_period + sex, data=data_cr_agg,FUN=sum)          # aggregating ocross multiple cancers (e.g. 'head and neck', 'other NADCs')

if(tumor_restriction=="other NADCs")
  data_res_mrg <- subset(data_res_mrg,!shcsid%in%c(40325,26686))    # bricolage to remove patients who develop cancer the same month as enrollment
if(tumor_restriction=="lung")
  data_res_mrg <- subset(data_res_mrg,shcsid!=19521)

X <- subset(data_res_mrg,stop_yar_full-start_yar_full>=31)   # removing patients with less than a month of follow-up
data_inc <- with(X,data.frame(shcsid=shcsid,start_yar=start_yar,stop_yar=stop_yar,start_yar_abs=start_yar,stop_yar_abs=stop_yar,cr_dinc=cr_dinc_y,
                                         dob=hiv_dob_y,registry=registry,linked=as.logical(linked),sex=sex))
if(!remove_ties)
  data_inc <- transform(data_inc,stop_yar=ifelse(start_yar==stop_yar & linked,stop_yar+0.5,stop_yar))
data_inc <- subset(data_inc,stop_yar>start_yar)     

data_inc <- transform(data_inc,sex=as.character(sex))
data_inc$sex[data_inc$sex=="Male"] <- "males"
data_inc$sex[data_inc$sex=="Female"] <- "females"
data_inc <- transform(data_inc,sex=factor(sex,levels=c("males","females")))

if(which_gender!="both")
  data_inc <- subset(data_inc,sex==which_gender)

data_inc <- transform(data_inc,start_yar_abs=start_yar,stop_yar_abs=stop_yar)
if(lower_cutoff_year==1988)
{
 data_timesplit <- survSplit(data_inc,cut=c(1996,2002,2008),start="start_yar",end="stop_yar",event="linked",episode="i")
 data_timesplit$time_period <- NA
 data_timesplit$time_period[data_timesplit$i==1] <- "[1988,1996)"
 data_timesplit$time_period[data_timesplit$i==2] <- "[1996,2002)"
 data_timesplit$time_period[data_timesplit$i==3] <- "[2002,2008)"
 data_timesplit$time_period[data_timesplit$i==4] <- "[2008,2013)"
}
if(lower_cutoff_year==1996)
{
  data_timesplit <- survSplit(data_inc,cut=c(2002,2008),start="start_yar",end="stop_yar",event="linked",episode="i")
  data_timesplit$time_period <- NA
  data_timesplit$time_period[data_timesplit$i==1] <- "[1996,2002)"
  data_timesplit$time_period[data_timesplit$i==2] <- "[2002,2008)"
  data_timesplit$time_period[data_timesplit$i==3] <- "[2008,2013)"
}
data_timesplit$i <- NULL
data_timesplit <- transform(data_timesplit,start_age=start_yar-dob,stop_age=stop_yar-dob,time_period=factor(time_period))

data_timesplit <- survSplit(data_timesplit,cut=c(35,45,55),start="start_age",end="stop_age",event="linked",episode="i")
data_timesplit$age_grp <- NA
data_timesplit$age_grp[data_timesplit$i==1] <- "[0,35)"
data_timesplit$age_grp[data_timesplit$i==2] <- "[35,45)"
data_timesplit$age_grp[data_timesplit$i==3] <- "[45,55)"
data_timesplit$age_grp[data_timesplit$i==4] <- "[55,Inf)"
data_timesplit$i <- NULL
data_timesplit <- transform(data_timesplit,yar=stop_age-start_age,age_grp=factor(age_grp),time_period=factor(time_period))

if(time_adj)
{data_inc_agg <- aggregate(cbind(linked,yar) ~ age_grp + time_period + sex, data=data_timesplit, FUN=sum)} else
{data_inc_agg <- aggregate(cbind(linked,yar) ~ age_grp + sex, data=data_timesplit, FUN=sum)}
data_inc_agg <- transform(data_inc_agg,rate=linked/yar)

if(tumor_restriction=="prostate")
  data_inc_agg <- subset(data_inc_agg,sex=="males")

if(any(data_cr_agg[c("age_grp","time_period","sex")]!=data_inc_agg[c("age_grp","time_period","sex")]))
  stop("the factors don't match")

data_inc_agg$exp_events <- data_cr_agg$rate*data_inc_agg$yar

# single SIR, age- and time-period- and sex- standardized

nb_obs_events <- sum(data_inc_agg$linked)
print(paste("occurrences in SHCS:",nb_obs_events))
nb_exp_events <- sum(data_inc_agg$exp_events)
print(paste("number of expected events",round(nb_exp_events,digits=round_digits)))
preg_single <- glm(nb_obs_events ~ offset(log(nb_exp_events)),family=poisson)

SIR <- nb_obs_events/nb_exp_events
ef <- exp(qnorm(0.975)/sqrt(nb_obs_events))
CI_l <- SIR/ef
CI_u <- SIR*ef
p_value <- summary(preg_single)$coefficients[4]

print(paste("single SIR:",round(SIR,digits=round_digits),"[",round(CI_l,digits=round_digits),",",round(CI_u,digits=round_digits),"]"))
print(paste("p-value:",round(p_value,digits=3)))
print("*-*-*-*-*-*")
