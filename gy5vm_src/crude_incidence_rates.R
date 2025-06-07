# crude incidence rates, without removing patients that are missing covariates

tumor_restriction <- "liver"
multiplier <- 10000
round_digits <- 2
filepath_read_R <- "E:\\P3RL\\R"

load(file=paste(filepath_read_R,"\\res_df_6m_nores_",tumor_restriction,"_1996.RData",sep=''))

res_df_linked <- subset(res_df,linked==1)
status_at_cr <- subset(res_df_linked,!duplicated(res_df_linked$shcsid,fromLast=TRUE)) 
cases <- nrow(status_at_cr)
py <- round(nrow(res_df)/12)

print(data.frame(py_ar=py,cases=cases,
                 inc=round(multiplier*cases/py,digits=round_digits),
                 CI_l=round(multiplier*qchisq(0.025,2*cases)/(2*py),digits=round_digits),
                 CI_u=round(multiplier*qchisq(0.975,2*(cases+1))/(2*py),digits=round_digits)))