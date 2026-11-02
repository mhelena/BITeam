#Power calculator for burden test

#N = 20000 #sample size
# maf = 0.05 
# beta = log(2.1) #beta is log of odds ratio
# #approximately GRR for a low prevalence disease
# cc_rat = 0.5 #case control ratio; optimal power at 0.5; simulate multiple scenarios

sig_thres = qchisq(1e-4, df = 1, ncp = 0, lower = F) # set a threshold based on test in the gene panel

##Functions
ncp_est = function(N_inp, beta_inp, maf_inp, cc_inp){
  return(2*maf_inp*(1-maf_inp)*N_inp*cc_inp*(1-cc_inp)*beta_inp^2)
}

pwr_est = function(thresh, N, OR, maf, cc_rat){
  beta = log(OR)
  pwr = pchisq(sig_thres, df = 1, ncp = ncp_est(N, beta, maf, cc_rat), lower = FALSE)
  
  return(cbind.data.frame("Sample_Size" = N, "Power" = pwr, "OR" = OR, "CC_ratio" = cc_rat))
}

ns = seq(100, 3000, 50)
ors = seq(1, 4, 0.25)
ccs = seq(0.1, 0.5, 0.1)
maf = 0.01
df_sim <- expand.grid(Var1 = ns, Var2 = ors, Var3 = ccs)

library(purrr)
res_cmb = pmap(df_sim, function(Var1, Var2, Var3) {
  # You can explicitly name the columns as arguments
  pwr_est(sig_thres, N = Var1, OR = Var2 , maf, cc_rat = Var3)
})
                
                
#plt_df = lapply(ns, function(x) pwr_est(sig_thres, x, beta, maf, cc_rat))

plt_df_cmb = do.call("rbind.data.frame", res_cmb)


ggplot(plt_df_cmb[plt_df_cmb$OR == 4,], aes(x = Sample_Size, y = Power, colour = factor(CC_ratio))) + 
         geom_line() + theme_minimal()
       
       