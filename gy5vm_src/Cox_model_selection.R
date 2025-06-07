# model comparisons, first run 'risk_factors_Cox.R' with skip_analysis=TRUE to create 'res_df' dataset

ptm <- proc.time()

print(paste("tumor:",tumor_restriction))

immuno_factors <- c("cd4","cd8","cd4Rcd8","rna")
#immuno_factors <- "cd4"
AIC_vect <- c()
creg_list <- list()
for(i in 1:length(immuno_factors))
{
 cand_factors <- paste(immuno_factors[i],c("_lag12","_lag24","_lag36","_ma12_lag12","_ma12_lag24","_ma24_lag12"),"_cat",sep='')
 adj_factors <- c("age_at_reg_cat","ever_smoke","risk_cat","hepC_RNA_lag12","hepB_pos_lag12","edu_cat","date_cat")
 #adj_factors <- c("age_at_reg_cat","ever_smoke","risk_cat","hepC_RNA_lag12","edu_cat","date_cat")
 res_df_i <- na.omit(res_df[c("tstart","tstop","linked",cand_factors,adj_factors)])
 #adj_factors <- c("strata(age_at_reg_cat)","ever_smoke","risk_cat","hepC_pos_at_reg","edu_cat")
 AIC_vect_i <- rep(NA,length(cand_factors))
 for(j in 1:length(cand_factors))
 {
  full_factors <- c(cand_factors[j],adj_factors)
  creg <- coxph(as.formula(paste("Surv(tstart,tstop,linked) ~",paste(full_factors,collapse='+'))),data=res_df_i)
  AIC_vect_i[j] <- extractAIC(creg)[2]
  creg_list[[(i-1)*length(cand_factors)+j]] <- creg
 }
 names(AIC_vect_i) <- cand_factors
 AIC_vect <- c(AIC_vect,AIC_vect_i)
}

print(round(AIC_vect,digits=1))
print(paste(c("adjusting for:",adj_factors),collapse="  "))
if(tumor_restriction!="lung" && !"hepB_pos_lag12"%in%adj_factors)
  warning("not adjusting for hepB")

write.table(round(AIC_vect),file=paste(filepath_code,"\\AIC_vect.txt",sep=''),row.names=FALSE,col.names=FALSE)

print(proc.time()-ptm)