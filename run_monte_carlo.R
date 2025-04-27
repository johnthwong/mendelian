source("code.R")

m_max = 500
mc_results = matrix(rep(NA, m_max*4), ncol = 4)
set.seed(1000)
for(m in 1:m_max){
  df = gen_solow_panel()
  tsls = ivreg::ivreg(
    formula = y ~ d + d_l1 + y_l1 | z + d_l1 + y_l1,
    data = df
  )
  mc_results[m,] = tsls$coefficients
}
colnames(mc_results) = c("intercept", "d", "d_l1", "y_l1")
write.csv(mc_results, "mc_results.csv")
