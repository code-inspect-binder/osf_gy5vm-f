# patient characteristics at cancer incidence (or end of study for cancer-free patients), requires pre-made 'res_df' file for cancer in question

library(tableone)
library(dplyr)

tumor_restriction <- "no cancer"       # "no cancer" for looking at patients who did not develop cancer, also sets complement <- TRUE
complement <- FALSE                    # TRUE for looking at characteristics of patients that DIDN'T develop the tumor(s)
write_table <- TRUE

filepath_code <- "C:\\HomeDir\\P3RL\\R_code"
filepath_read_R <- "E:\\P3RL\\R"

if(tumor_restriction=="no cancer")
{
  complement <- TRUE      
  load(file=paste(filepath_read_R,"\\res_df_12m_nores_any_1996.RData",sep=''))
 } else
{load(file=paste(filepath_read_R,"\\res_df_12m_nores_",tumor_restriction,"_1996.RData",sep=''))}

status_at_cr <- subset(res_df,!duplicated(res_df$shcsid,fromLast=TRUE))          # only keeping last record
status_at_cr <- subset(status_at_cr,linked==!complement)                         # only keep those linked (or those not linked...)

status_at_cr <- transform(status_at_cr,cd4_cat=cut(cd4,breaks=c(0,200,350,500,Inf),right=FALSE),cd8_cat=cut(cd8,breaks=c(0,1000,Inf),right=FALSE),
                          cd4Rcd8_cat=cut(cd4/cd8,breaks=c(0,0.5,1,Inf),right=FALSE),rna_cat=cut(rna,breaks=c(0,50,1000,Inf),right=FALSE))

#status_at_cr <- transform(status_at_cr,age_cat=cut(age,breaks=c(0,35,45,55,Inf),right=FALSE),cd4_cat=cut(cd4,breaks=c(0,100,200,500,Inf),right=FALSE),
#                          cd8_cat=cut(cd8,breaks=c(0,600,1000,Inf),right=FALSE),rna_cat=cut(rna,breaks=c(0,500,Inf),right=FALSE))
status_at_cr <- transform(status_at_cr,cd4_cat=factor(cd4_cat,levels=c(levels(cd4_cat),"missing")),cd8_cat=factor(cd8_cat,levels=c(levels(cd8_cat),"missing")),
                                       cd4Rcd8_cat=factor(cd4Rcd8_cat,levels=c(levels(cd4Rcd8_cat),"missing")),
                                       rna_cat=factor(rna_cat,levels=c(levels(rna_cat),"missing")),ever_smoke=factor(ever_smoke,levels=c("no","yes","missing")))
status_at_cr$cd4_cat[is.na(status_at_cr$cd4_cat)] <- "missing"
status_at_cr$cd8_cat[is.na(status_at_cr$cd8_cat)] <- "missing"
status_at_cr$cd4Rcd8_cat[is.na(status_at_cr$cd4Rcd8_cat)] <- "missing"
status_at_cr$rna_cat[is.na(status_at_cr$rna_cat)] <- "missing"
status_at_cr$ever_smoke[is.na(status_at_cr$ever_smoke)] <- "missing"
status_at_cr$region_cat[is.na(status_at_cr$region_cat)] <- "Other"
status_at_cr$sex_cat <- factor(status_at_cr$sex_cat,levels=c("Male","Female"))
#status_at_cr$risk_cat <- dplyr::recode(status_at_cr$risk_cat,'Other men'="Other",'Other women'="Other")

risk_factors <- c("age","sex_cat","risk_cat","haart_lag6","cd4_cat","cd8_cat","cd4Rcd8_cat","rna_cat","ever_smoke","edu_cat","hepC_RNA","hepB_pos","CDC_cat")

print(paste("cancer:",tumor_restriction))
tab <- print(CreateTableOne(vars=risk_factors,data=status_at_cr))
if(write_table)
 write.table(tab,file=paste(filepath_code,"\\table.txt",sep=''),sep=";",col.names=FALSE,quote=FALSE)

# to copy .txt table to Word, select whole text in .txt file, copy directly into Word, then select the text, Insert -> Table -> Convert text to table

