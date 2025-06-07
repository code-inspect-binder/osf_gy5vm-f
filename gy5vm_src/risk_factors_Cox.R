# Cox regression using res_df dataframe (created by write_res_df.R)
# For crude rates, see 'crude_incidence_rates.R'
# chi_sq test -> lab8.R
# model selection -> Cox_model_selection.R

library(survival)
library(survminer)
library(Hmisc)
library(mitools)    # for pooling of imputations
library(dplyr)

ptm <- proc.time()

filepath_code <- "C:\\HomeDir\\P3RL\\R_code"
filepath_read_R <- "E:\\P3RL\\R"
start_year <- 1996                    # start of study period: 1988 or 1996
skip_analysis <- FALSE                # if TRUE, abort immediately after construction of res_df object, e.g. for model selection
cat_agg <- TRUE                       # whether to aggregate low-pop categories?
univariate_reg <- FALSE               # perform the univariate regressions for each factor? -> reg_list
load_res_df <- TRUE                   # whether to reload the dataframe, important to do this when switching between cancers
diagnostics <- FALSE                  # proportional hazards test, Martingale residuals, ...
linkage_only <- FALSE                 # whether to remove tumors reported in SHCS but not linked to cancer registries
info_source <- "both"                 # which source of cancer information to use: "linkage" (cancer registries), "SHCS", or "both"
nb_imp <- 5                           # number of imputations (up to 10)

tumor_restriction <- "liver"

if(load_res_df)
{
  load(file=switch(info_source, "both"=paste(filepath_read_R,"\\res_df_12m_nores_",tumor_restriction,"_",start_year,".RData",sep=''),
                                "linkage"=paste(filepath_read_R,"\\res_df_12m_nores_",tumor_restriction,"_",start_year,"_linkage_only.RData",sep=''),
                                "SHCS"=paste(filepath_read_R,"\\res_df_12m_nores_",tumor_restriction,"_",start_year,"_SHCS_only.RData",sep='')))
} else
{warning("res_df NOT reloaded")}
print(paste("data_source:",info_source))

base_factors <- c("cd4_ma24_lag12_cat","cd4Rcd8_lag24_cat")    # for ordinal variables, replace the '_cat' with '_ord'
#base_factors <- c("cd8_lag36_cat")
#base_factors <- c("rna_lag24_ord")
#base_factors <- c()

if(tumor_restriction=="lung")
{adj_factors <- c("age_at_reg_cat","ever_smoke_imp","risk_cat","hepC_RNA_lag12","edu_cat","date_cat")} else
{
  adj_factors <- c("age_at_reg_cat","ever_smoke_imp","risk_cat","hepC_RNA_lag12","hepB_pos_lag12","edu_cat","date_cat")
  #adj_factors <- c("age_at_reg_cat","ever_smoke_imp","risk_cat","hepB_pos_lag12","edu_cat","date_cat")
  #adj_factors <- c("age_at_reg_cat","risk_cat","edu_cat","date_cat")
}

strat_factors <- c()
if(tumor_restriction=="lung")
 strat_factors <- c("age_at_reg_cat")
if(tumor_restriction=="prostate")
 strat_factors <- c("age_at_reg_cat")

for(i in 1:length(strat_factors))
  adj_factors[adj_factors==strat_factors[i]] <- paste("strata(",strat_factors[i],")",sep='')

risk_factors <- c(base_factors,adj_factors)
res_df$mar <- 1

if(load_res_df)
{
  res_df$cd4Rcd8[res_df$cd4Rcd8==Inf] <- 20
  
 if(start_year==1988)
  period_breaks <- c(1988,1996,2002,2008,2013)
 if(start_year==1996)
  period_breaks <- c(1996,2002,2008,2013)
 age_breaks <- c(30,40,50)
 age_at_reg_breaks <- c(30,40,50)
 cd4_breaks <- c(200,350,500)
 cd8_breaks <- c(1000)
 cd4Rcd8_breaks <- c(0.5,1)
 rna_breaks <- c(50,500)

 
 if(tumor_restriction=="anus")
   res_df$risk_cat <- factor(res_df$risk_cat,levels=c("Other men","IDU","MSM","Other women"))
 if(tumor_restriction=="lung" && cat_agg)
 {
   if(start_year==1988)
     period_breaks <- c(1988,2002,2008,2013)
   if(start_year==1996)
     period_breaks <- c(1996,2002,2008,2013)
   age_breaks <- c(45,55)
   cd4_breaks <- c(200,500)
   rna_breaks <- 500
 }
 if(tumor_restriction=="head and neck" && cat_agg)
 {
   age_breaks <- c(45,55)
   cd4_breaks <- c(200,500)
 }
 if(tumor_restriction=="prostate" && cat_agg)
 {
   age_breaks <- c(60,70)
   age_at_reg_breaks <- c(40,50,60)
   res_df$risk_cat <- factor(res_df$risk_cat,levels=c("MSM","IDU","Other men"))
   #cd4_breaks <- c(200,500)
   #res_df$edu_cat[res_df$edu_cat=="other/unknown"] <- NA
 }
 if(tumor_restriction=="liver" && cat_agg)
 {
   age_at_reg_breaks <- c(45,55)
   cd4_breaks <- c(200,500)
   res_df$edu_cat <- as.character(res_df$edu_cat)
   res_df$edu_cat[res_df$edu_cat%in%c("class 2","class 3")] <- "class 2 or 3"
   res_df$edu_cat <- factor(res_df$edu_cat,levels=c("class 1","class 2 or 3","other/unknown"))
   res_df$risk_cat <- as.character(res_df$risk_cat)
   res_df$risk_cat[res_df$risk_cat%in%c("Other women","Other men")] <- "Other men/women"
   res_df$risk_cat <- factor(res_df$risk_cat,levels=c("MSM","IDU","Other men/women"))
 }
 if(tumor_restriction=="HL")
{
  period_breaks <- c(1996,2002,2013)
  age_at_reg_breaks <- c(40,50)
  cd4_breaks <- c(200,350)
  res_df$risk_cat <- as.character(res_df$risk_cat)
  res_df$risk_cat[res_df$risk_cat%in%c("Other women","Other men")] <- "Other men/women"
  res_df$risk_cat <- factor(res_df$risk_cat,levels=c("MSM","IDU","Other men/women"))
}
 
 res_df <- transform(res_df,age_cat=cut(age,breaks=c(0,age_breaks,Inf),right=FALSE),age_at_reg_cat=cut(age_at_reg,breaks=c(0,age_at_reg_breaks,Inf),right=FALSE),date_cat=cut(date_y,breaks=period_breaks,right=FALSE),
                            cd4_cat=cut(cd4,breaks=c(0,cd4_breaks,Inf),right=FALSE),cd8_cat=cut(cd8,breaks=c(0,cd8_breaks,Inf),right=FALSE),cd4Rcd8_cat=cut(cd4Rcd8,breaks=c(0,cd4Rcd8_breaks,Inf),right=FALSE),
                            rna_cat=cut(rna,breaks=c(0,rna_breaks,Inf),right=FALSE))

 res_df <- transform(res_df,cd4_lag12_cat=cut(cd4_lag12,breaks=c(0,cd4_breaks,Inf),right=FALSE),cd4_lag24_cat=cut(cd4_lag24,breaks=c(0,cd4_breaks,Inf),right=FALSE),cd4_lag36_cat=cut(cd4_lag36,breaks=c(0,cd4_breaks,Inf),right=FALSE),
                            cd8_lag12_cat=cut(cd8_lag12,breaks=c(0,cd8_breaks,Inf),right=FALSE),cd8_lag24_cat=cut(cd8_lag24,breaks=c(0,cd8_breaks,Inf),right=FALSE),cd8_lag36_cat=cut(cd8_lag36,breaks=c(0,cd8_breaks,Inf),right=FALSE),
                            cd4Rcd8_lag12_cat=cut(cd4Rcd8_lag12,breaks=c(0,cd4Rcd8_breaks,Inf),right=FALSE),cd4Rcd8_lag24_cat=cut(cd4Rcd8_lag24,breaks=c(0,cd4Rcd8_breaks,Inf),right=FALSE),cd4Rcd8_lag36_cat=cut(cd4Rcd8_lag36,breaks=c(0,cd4Rcd8_breaks,Inf),right=FALSE),
                            rna_lag12_cat=cut(rna_lag12,breaks=c(0,rna_breaks,Inf),right=FALSE),rna_lag24_cat=cut(rna_lag24,breaks=c(0,rna_breaks,Inf),right=FALSE),rna_lag36_cat=cut(rna_lag36,breaks=c(0,rna_breaks,Inf),right=FALSE),
                            cd4_ma12_lag12_cat=cut(cd4_ma12_lag12,breaks=c(0,cd4_breaks,Inf),right=FALSE),cd4_ma12_lag24_cat=cut(cd4_ma12_lag24,breaks=c(0,cd4_breaks,Inf),right=FALSE),cd4_ma24_lag12_cat=cut(cd4_ma24_lag12,breaks=c(0,cd4_breaks,Inf),right=FALSE),
                            cd8_ma12_lag12_cat=cut(cd8_ma12_lag12,breaks=c(0,cd8_breaks,Inf),right=FALSE),cd8_ma12_lag24_cat=cut(cd8_ma12_lag24,breaks=c(0,cd8_breaks,Inf),right=FALSE),cd8_ma24_lag12_cat=cut(cd8_ma24_lag12,breaks=c(0,cd8_breaks,Inf),right=FALSE),
                            cd4Rcd8_ma12_lag12_cat=cut(cd4Rcd8_ma12_lag12,breaks=c(0,cd4Rcd8_breaks,Inf),right=FALSE),
                            cd4Rcd8_ma12_lag24_cat=cut(cd4Rcd8_ma12_lag24,breaks=c(0,cd4Rcd8_breaks,Inf),right=FALSE),cd4Rcd8_ma24_lag12_cat=cut(cd4Rcd8_ma24_lag12,breaks=c(0,cd4Rcd8_breaks,Inf),right=FALSE),
                            rna_ma12_lag12_cat=cut(rna_ma12_lag12,breaks=c(0,rna_breaks,Inf),right=FALSE),rna_ma12_lag24_cat=cut(rna_ma12_lag24,breaks=c(0,rna_breaks,Inf),right=FALSE),rna_ma24_lag12_cat=cut(rna_ma24_lag12,breaks=c(0,rna_breaks,Inf),right=FALSE))

 res_df <- transform(res_df,cd4_lag12_ord=scale(as.numeric(cd4_lag12_cat),scale=FALSE),cd4_lag24_ord=scale(as.numeric(cd4_lag24_cat),scale=FALSE),cd4_lag36_ord=scale(as.numeric(cd4_lag36_cat),scale=FALSE),
                            cd4_ma12_lag12_ord=scale(as.numeric(cd4_ma12_lag12_cat),scale=FALSE),cd4_ma12_lag24_ord=scale(as.numeric(cd4_ma12_lag24_cat),scale=FALSE),cd4_ma24_lag12_ord=scale(as.numeric(cd4_ma24_lag12_cat),scale=FALSE),
                            cd8_lag12_ord=scale(as.numeric(cd8_lag12_cat),scale=FALSE),cd8_lag24_ord=scale(as.numeric(cd8_lag24_cat),scale=FALSE),cd8_lag36_ord=scale(as.numeric(cd8_lag36_cat),scale=FALSE),
                            cd8_ma12_lag12_ord=scale(as.numeric(cd8_ma12_lag12_cat),scale=FALSE),cd8_ma12_lag24_ord=scale(as.numeric(cd8_ma12_lag24_cat),scale=FALSE),cd8_ma24_lag12_ord=scale(as.numeric(cd8_ma24_lag12_cat),scale=FALSE),
                            cd4Rcd8_lag12_ord=scale(as.numeric(cd4Rcd8_lag12_cat),scale=FALSE),cd4Rcd8_lag24_ord=scale(as.numeric(cd4Rcd8_lag24_cat),scale=FALSE),cd4Rcd8_lag36_ord=scale(as.numeric(cd4Rcd8_lag36_cat),scale=FALSE),
                            cd4Rcd8_ma12_lag12_ord=scale(as.numeric(cd4Rcd8_ma12_lag12_cat)),cd4Rcd8_ma24_lag12_ord=scale(as.numeric(cd4Rcd8_ma24_lag12_cat)),cd4Rcd8_ma12_lag24_ord=scale(as.numeric(cd4Rcd8_ma12_lag24_cat)),
                            rna_lag12_ord=scale(as.numeric(rna_lag12_cat)),rna_lag24_ord=scale(as.numeric(rna_lag24_cat)),rna_ma12_lag12_ord=scale(as.numeric(rna_ma12_lag12_cat)),rna_ma24_lag12_ord=scale(as.numeric(rna_ma24_lag12_cat)))

}

agg_strat <- function(factors)                            # aggregation over stratas (mainly for debugging purposes)
{
  form  <- paste("cbind(linked,mar) ~",paste(factors,collapse="+"))
  x <- aggregate(as.formula(form),data=res_df,FUN=sum)
  x$rate <- x$linked/x$mar
  x
}

if(skip_analysis)
{
  print(proc.time()-ptm)
  stop("program intentionally aborted")
}
 
if("risk_cat_MSM"%in%risk_factors)
  res_df$risk_cat_MSM <- res_df$risk_cat=="MSM"
if("risk_cat_IDU"%in%risk_factors)
  res_df$risk_cat_IDU <- res_df$risk_cat=="IDU"
if("cd4_ma24_lag12_below500"%in%risk_factors)
  res_df$cd4_ma24_lag12_below500 <- res_df$cd4_ma24_lag12<500
if("cd8_lag12_below600"%in%risk_factors)
  res_df$cd8_lag12_below600 <- res_df$cd8_lag12<600
if("cd4Rcd8_ma24_lag12_below0.4"%in%risk_factors)
  res_df$cd4Rcd8_ma24_lag12_below0.4 <- res_df$cd4Rcd8_ma24_lag12<0.4
if("cd4Rcd8_ma12_lag24_below0.4"%in%risk_factors)
  res_df$cd4Rcd8_ma12_lag24_below0.4 <- res_df$cd4Rcd8_ma12_lag24<0.4
if("cd4Rcd8_ma12_lag12_below0.4"%in%risk_factors)
  res_df$cd4Rcd8_ma12_lag12_below0.4 <- res_df$cd4Rcd8_ma12_lag12<0.4
if("cd4Rcd8_lag24_below0.4"%in%risk_factors)
  res_df$cd4Rcd8_lag24_below0.4 <- res_df$cd4Rcd8_lag24<0.4

reg_formula <- paste("Surv(tstart,tstop,linked) ~",paste(risk_factors,collapse="+"))
creg_list <- list()
#creg_dummy <- coxph(Surv(tstart,tstop,linked) ~ 1, data=res_df)

if(univariate_reg)
 for(i in 1:length(risk_factors))
{
   creg_i <- coxph(as.formula(paste("Surv(tstart,tstop,linked) ~ ",risk_factors[i])),data=res_df)
   creg_list[[i]] <- creg_i
 }

if(!"ever_smoke_imp"%in%risk_factors && !"strata(ever_smoke_imp)"%in%risk_factors)
{
 if("ever_smoke"%in%risk_factors || "strata(ever_smoke)"%in%risk_factors)
   warning("smoking NOT imputed")
 creg <- coxph(formula=as.formula(reg_formula), data=res_df)
 
 tbl <- cbind(round(summary(creg)$conf.int[,c(1,3,4)],digits=2),round(summary(creg)$coefficients[,5],digits=3))
 colnames(tbl)[4] <- "Wald p-value"
} else                                                                                
{
 res_df$ever_smoke_imp <- NULL
 load(file=paste(filepath_read_R,"\\smoking_imputations.RData",sep=''))                   # using imputed smoking values -> pooling Cox regressions
 mm <- match(res_df$shcsid,smoke_imp_matrix[,1])
 res_df_list <- list()
 for(i in 1:nb_imp)
 {
  res_df_temp <- res_df
  res_df_temp$ever_smoke_imp <- smoke_imp_matrix[mm,i+1]
  res_df_temp <- transform(res_df_temp,ever_smoke_imp=factor(ifelse(ever_smoke_imp==1,"yes","no"),levels=c("no","yes")))
  res_df_list[[i]] <- res_df_temp
 }
 imp_obj <- imputationList(res_df_list)
 creg_models <- with(imp_obj,coxph(formula=as.formula(reg_formula)))                # contains each individual fit
 creg_pooled <- summary(MIcombine(creg_models))
 creg <- creg_models[[5]]
 z <- creg_pooled$results/creg_pooled$se                                            # reverse-engineering p-values (not ideal)
 tbl <- data.frame(HR=round(exp(creg_pooled$results),digits=2),CI_l=round(exp(creg_pooled$'(lower'),digits=2),CI_u=round(exp(creg_pooled$'upper)'),digits=2),
                   p_value=round(exp(-0.717*abs(z)-0.416*(z^2)),digits=3))
 rm(z)
 rownames(tbl) <- rownames(creg_pooled)
}

print(paste("number of ",tumor_restriction," cancer cases: ",creg$nevent," (",sum(res_df$linked),") in ",creg$n," (",nrow(res_df),") person-months",sep=''))
print(tbl)
print(paste("stratified with respect to:",paste(strat_factors,collapse=", ")))
if(tumor_restriction!="lung" && !"hepB_pos_lag12"%in%risk_factors)
  warning("not adjusting for hepB")
if(diagnostics)
{
  creg_test <- cox.zph(creg)
  print(paste("Schoenfeld proportional hazards test (global):",round(creg_test$table[nrow(creg_test$table),3],digits=3)))
  dfbeta <- residuals(creg,type="dfbeta")
  martingres <- residuals(creg,type="martingale")
}


print(proc.time()-ptm)
