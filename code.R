library(tidyverse)
set.seed(1000)

# Parameters -------------------------------------------------------------------
i_max = 50
t_max = 1000
alpha = 0.5
delta = 0.2
beta_0 = 0.3
beta_1 = -0.4
gamma_1 = 0.6
theta_1 = 0.7


# Generate panel ---------------------------------------------------------------

## Generate correlated errors first --------------------------------------------

# Generate z sequence
df_z = data.frame(
  i = 1:i_max,
  z = rnorm(i_max)
)

# Define variance-covariance matrix Sigma
gen_correlated_errors <- function(){
  Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)
  epsilons <- MASS::mvrnorm(t_max, mu = c(0,0), Sigma = Sigma)
  return(epsilons)
}

## Set up loop to iteratively generate each variable ---------------------------
gen_labrecque_panel <- function(){
  df_panel = data.frame()
  for (i in 1:i_max){
    by_i = data.frame()
    z = df_z[i,"z"]
    epsilons = gen_correlated_errors()
    for (t in 1:t_max){
      if (t == 1){
        d = epsilons[t, 1] + delta*z # implicitly assume d = 0 at t = 0
        y = beta_0*d + epsilons[t, 2]
      } else {
        d_l1 = by_i[t-1, "d"]
        d = alpha*d_l1 + delta*z + epsilons[t, 1]
        y = beta_0*d + beta_1*d_l1 + epsilons[t, 2]
      }
      df = data.frame(
        t = t,
        d = d,
        y = y
      )
      by_i = rbind(by_i, df)
    }
    by_i = by_i %>% mutate(
      i = i,
      z = z
    )
    df_panel = rbind(df_panel, by_i)
  }
  df_panel = df_panel %>%
    group_by(i) %>%
    mutate(
      d_l1 = lag(d)
    ) %>%
    na.omit()
  return(df_panel)
}

df_panel <- gen_labrecque_panel()
# Estimation -------------------------------------------------------------------

lm_stage2 = df_panel %>% lm(
  formula = y ~ d + d_l1,
  data = .
)
summary(lm_stage2)
# these coefficients are incorrect

lm_stage1 = df_panel %>% lm(
  formula = d ~ d_l1 + z,
  data = .
)
summary(lm_stage1)
# these coefficients are correct

# df_panel_predicted <- data.frame(
#   d = predict(lm_stage1)
# ) %>%
#   mutate(d_l1 = lag(d)) %>%
#   cbind(df_panel %>% select(i, t, y, z)) %>%
#   na.omit()
# 
# tsls_twostep = df_panel_predicted %>% lm(
#   formula = y ~ d + d_l1,
#   data = .
# )
# summary(tsls_twostep)
# # doesn't work

tsls_d <- ivreg::ivreg(
  formula = y ~ d | z,
  data = df_panel
)
# beta_0 should be 0.3
summary(tsls_d)

# instead, it should be biased to be about this amount
# gamma_2 = delta*(1 + alpha)
beta_0 + beta_1*delta/(delta*alpha + delta)

# if we use d_l1 as the treatment instead
tsls_d_l1 <- ivreg::ivreg(
  formula = y ~ d_l1 | z,
  data = df_panel
)

# beta_1 should be -0.1
summary(tsls_d_l1)

# theoretically it should be biased to be about this amount
beta_1 + beta_0*(alpha + delta/delta)

ivreg::ivreg(
  formula = y ~ d | z + d_l1,
  data = df_panel
) %>% summary()
# this is getting closer.

tsls_dl1_as_control <- ivreg::ivreg(
  formula = y ~ d + d_l1 | z + d_l1,
  data = df_panel
)
summary(tsls_dl1_as_control)
# oh this is correctly estimating.

ivreg::ivreg(
  formula = y ~ d + d_l1 | z + lag(d_l1),
  data = df_panel
) %>% summary()
# using d_l1 on d_l2 also works.

# Tables -----------------------------------------------------------------------
true_model <- texreg::createTexreg(
  coef.names = c("(Intercept)", "d", "d_l1"),
  coef = c(0, beta_0, beta_1)
)

table_omitted <- texreg::texreg(
  l = list(true_model, tsls_d),
  custom.model.names = c("True", "TSLS"),
  custom.coef.names = c(
    "Intercept",
    "$d_t$",
    "$Ld_t$"
  ),
  include.rsquared = FALSE,
  include.adjrs = FALSE
)

table_lagd <- texreg::texreg(
  l = list(true_model, tsls_d, tsls_dl1_as_control),
  custom.model.names = c("True", "TSLS", "TSLS With Lag"),
  custom.coef.names = c(
    "Intercept",
    "$d_t$",
    "$Ld_t$"
  ),
  include.rsquared = FALSE,
  include.adjrs = FALSE
)

# Test with IRFs ---------------------------------------------------------------
delta_hat = lm_stage1$coefficients["log(z)"]
alpha_hat = lm_stage1$coefficients["d_l1"]

df_panel = df_panel %>%
  mutate(
    dd_dz = 1/(z + 1) * alpha_hat^(z + 1) - alpha_hat
  )
  

df_panel %>% lm(
  formula = y ~ d + d_l1 + dd_dz,
  data = .
) %>%
  summary()

df_panel %>% ivreg::ivreg(
  data = .,
  formula = y ~ d | dd_dz
) %>%
  summary()

# Generating a bi-directional relationship (no simultaneity tho) ---------------
gen_solow_panel <- function(){
  df_panel = data.frame()
  for (i in 1:i_max){
    by_i = data.frame()
    z = df_z[i,"z"]
    epsilons = gen_correlated_errors()
    for (t in 1:t_max){
      if (t == 1){
        d = epsilons[t, 1] + delta*z # implicitly assume d = 0 and y = 0 at t = 0
        y = beta_0*d + epsilons[t, 2]
      } else {
        d_l1 = by_i[t-1, "d"]
        y_l1 = by_i[t-1, "y"]
        d = alpha*d_l1 + delta*z + theta_1*y_l1 + epsilons[t, 1]
        y = beta_0*d + beta_1*d_l1 + gamma_1*y_l1 + epsilons[t, 2]
      }
      df = data.frame(
        t = t,
        d = d,
        y = y
      )
      by_i = rbind(by_i, df)
    }
    by_i = by_i %>% mutate(
      i = i,
      z = z
    )
    df_panel = rbind(df_panel, by_i)
  }
  df_panel = df_panel %>%
    group_by(i) %>%
    mutate(
      d_l1 = lag(d),
      y_l1 = lag(y)
    ) %>%
    na.omit()
  return(df_panel)
}

df_panel_lagy = gen_solow_panel()

tsls_lagy <- ivreg::ivreg(
  formula = y ~ d + d_l1 + y_l1 | z + d_l1 + y_l1,
  data = df_panel_lagy
)

summary(tsls_lagy)

lm_stage1_lagy <- lm(
  formula = d ~ z + d_l1 + y_l1,
  data = df_panel_lagy 
)
summary(lm_stage1_lagy)

lm_reduced_form_lagy <- lm(
  formula = y ~ z + d_l1 + y_l1,
  data = df_panel_lagy 
)
summary(lm_reduced_form_lagy)

# the coefficient on z is close to
delta*beta_0

# intentionally omit the lagged terms
tsls_solow_omitted_stage2 <- ivreg::ivreg(
  formula = y ~ d | z,
  data = df_panel_lagy
) 

tsls_solow_omitted_stage1 <- lm(
  formula = d ~ z,
  data = df_panel_lagy
) 
# the coefficient becomes negative

# table ------------------------------------------------------------------------
true_solow_model_stage1 <- texreg::createTexreg(
  coef.names = c("(Intercept)", "z", "d_l1", "y_l1"),
  coef = c(0, delta, alpha, theta_1)
)

true_solow_model_stage2 <- texreg::createTexreg(
  coef.names = c("(Intercept)", "d", "d_l1", "y_l1"),
  coef = c(0, beta_0, beta_1, gamma_1)
)

table_solow <- texreg::texreg(
  l = list(
    true_solow_model_stage1, 
    tsls_solow_omitted_stage1,
    lm_stage1_lagy,
    true_solow_model_stage2,
    tsls_solow_omitted_stage2,
    tsls_lagy
    ),
  custom.header = list("1st-Stage ($d_t$)" = 1:3, "2nd-Stage ($y_t$)" = 4:6),
  custom.model.names = c("True", "(a)", "(b)", "True", "(c)", "(d)"),
  custom.coef.names = c(
    "Intercept",
    "$z_t$",
    "$Ld_t$",
    "$Ly_t$",
    "$d_t$"
  ),
  include.rsquared = FALSE,
  include.adjrs = FALSE,
  no.margin = TRUE
)

# monte-carlo analysis ---------------------------------------------------------
mc_coef = read.csv("mc_results.csv")

plot_density <- function(x){
  mu <- mean(x)
  sigma <- sd(x)
  
  # Compute density
  dens <- density(x)
  
  # Convert to data frame
  dens_df <- data.frame(x = dens$x, y = dens$y)
  
  # Plot
  ggplot(dens_df, aes(x = x, y = y)) +
    # Shade the ±1 SD region
    geom_area(data = subset(dens_df, x >= mu - 2*sigma & x <= mu + 2*sigma),
              fill = "skyblue", alpha = 0.5) +
    # Full density curve
    geom_line(color = "black", linewidth = 0.5)+
    # Add vertical lines for mean and ±2 SD
    geom_vline(xintercept = beta_0, color = "black", linetype = "dashed")  +
    labs(title = "",
         x = "Estimated Coefficient", y = "Density") +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18)
    )
}

plot_unbiased = plot_density(mc_coef$d)

mc_coef_biased = read.csv("mc_results_biased.csv")
plot_biased = plot_density(mc_coef_biased$d)

mc_coef_labreque = read.csv("mc_results_labrecque.csv")
plot_labrecque = plot_density(mc_coef_labreque$d)



# generating y and d simultaneously --------------------------------------------
# i should know the first y_{t-1} = 0. the subsequent y_{t-1} are known. 
# generate d_t (z_t) first.