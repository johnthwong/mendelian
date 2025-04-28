source("code.R")
# parameters
m_max = 500
t_max = 100

# one instance
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

# biased example
mc_results = matrix(rep(NA, m_max*2), ncol = 2)
for(m in 1:m_max){
  df = gen_solow_panel()
  tsls = ivreg::ivreg(
    formula = y ~ d | z,
    data = df
  )
  mc_results[m,] = tsls$coefficients
}
colnames(mc_results) = c("intercept", "d")
write.csv(mc_results, "mc_results_biased.csv")

# lagged treatment only
mc_results = matrix(rep(NA, m_max*3), ncol = 3)
set.seed(1000)
for(m in 1:m_max){
  df = gen_labrecque_panel()
  tsls = ivreg::ivreg(
    formula = y ~ d + d_l1 | z + d_l1,
    data = df
  )
  mc_results[m,] = tsls$coefficients
}
colnames(mc_results) = c("intercept", "d", "d_l1")
write.csv(mc_results, "mc_results_labrecque.csv")
