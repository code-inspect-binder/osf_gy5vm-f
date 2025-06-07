# multiple imputation for smoking, run pilot.R first to load data_pat
# about 6 minutes for 5 imputations

library(mice)

ptm <- proc.time()

nb_actions <- 5                  # number of imputations for pooling
data_pat_imp <- data_pat

data_pat_imp$cv_smoked[data_pat_imp$cv_smoked==9] <- NA
data_pat_imp$cv_smoked <- data_pat_imp$cv_smoked==1
data_pat_imp$ethnicity[data_pat_imp$ethnicity==9] <- NA
data_pat_imp$education[data_pat_imp$education==9] <- NA
data_pat_imp$profession[data_pat_imp$profession==9] <- NA
data_pat_imp$risk[data_pat_imp$risk==9] <- NA

data_pat_imp <- transform(data_pat_imp,age=as.Year(regdate)-born,sex=factor(sex),ethnicity=factor(ethnicity),education=factor(education),profession=factor(profession),risk=factor(risk))

data_pat_imp <- data_pat_imp[c("ethnicity","education","profession","risk","sex","age","cv_smoked")]

#p_matrix <- matrix(0,nrow=ncol(data_pat_imp),ncol=ncol(data_pat_imp))
#p_matrix[2:ncol(data_pat_imp),1] <- 1

smoke_imp <- mice(data_pat_imp,m=nb_actions)

smoke_imp_matrix <- cbind(data_pat$shcsid)
for(i in 1:smoke_imp$m)
  smoke_imp_matrix <- cbind(smoke_imp_matrix,complete(smoke_imp,action=i)$cv_smoked)

save(smoke_imp_matrix,file=paste(filepath_read_R,"\\smoking_imputations.RData",sep=''))         # patient ids in first column, imputations in remaining columns

print(proc.time()-ptm)

