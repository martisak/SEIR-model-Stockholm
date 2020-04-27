


# Remove workspace
rm(list = ls(all=TRUE))

heading <- "
#------------------------------------------------------------------------------------
#
# FILE:         Estimate_SEIR
# CONTACT:      Lisa Brouwers
# EMAIL:        analysenheten@folkhalsomyndigheten.se
# AFFILIATION:  FD-AN
# PROJECT:      Modelling COVID-19
#
#
# Created:      2020-03-15 
# Updated:      2020-04-22
# R version:    3.5.2
#
# What the script does: This script estimates the parameters of the infectivity described in the report 
#                       'Skattning av peakdagen och antal infekterade i covid-19-utbrottet i Stockholms län'.
#                       With these we are able to estimate the number of infectious individuals at different time points, 
#                       the cumulative number of infected, and the estimated effective basic reproduction number. 
#                       If you want to reproduce the results obtained, or change anything, first write the  
#                       directory where the project folder is saved on your own computer in poject.path below.
#                       
#                       
#                       
#------------------------------------------------------------------------------------
\n"

#---------------------------------------------------------------------------------------------------
# LIBRARIES
#---------------------------------------------------------------------------------------------------
library(reshape2)
library(openxlsx)     # to write tables in excel
library(RColorBrewer)
library(rootSolve)    # to load function multiroot that finds roots to system of equations
library(deSolve)
library(tidyverse)
library(glue)
library(rebus)
library(patchwork)
library(furrr)

plan(multiprocess)

Sys.setenv(LANGUAGE = "en")
Sys.setlocale("LC_TIME", "English")

#-----------------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------------

# input your own path where folder (end with "/") is located in project.path.
# e.g. project.path 	     <- "C:/Users/Modelling/"

project.path 	     <- ""


data.path 		     <- paste0(project.path, "Data")
script.path 	     <- paste0(project.path, "Script")
table.path         <- paste0(project.path, "Results/Tables")
figure.path        <- paste0(project.path, "Results/Figures")


#---------------------------------------------------------------------------------------------------
# Basic functions and set.up
#---------------------------------------------------------------------------------------------------

Blues <- brewer.pal(n = 9, name = "Blues")

roundUp <- function(x, level = c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * level[[which(x <= 10^floor(log10(x)) * level)[[1]]]]
}



## From a sample x, calculate a credibility interval with wanted level. 
CRI <- function(x, level = 0.95){
  n <- length(x)
  L = (1 - level) / 2
  U = 1 - (1 - level) / 2
  x <- sort(x)
  resL <- x[n * L] 
  resU <- x[n * U]
  
  return(c(resL,resU))
}

CRI_90_low <- function(x){
  return(CRI(x, level = 0.9)[1])
}
CRI_90_up <- function(x){
  return(CRI(x, level = 0.9)[2])
}

CRI_95_low <- function(x){
  return(CRI(x, level = 0.95)[1])
}
CRI_95_up <- function(x){
  return(CRI(x, level = 0.95)[2])
}



#---------------------------------------------------------------------------------------------------
# Read data. 
#---------------------------------------------------------------------------------------------------

# This is the data used in the analysis. It differs from some of the reported case data since we 
# removed imported cases. 

Stockholm_Data_10_april <- read_delim(
  paste0(data.path, "/Data_2020-04-10Ny.txt"),
  delim = " ", col_names = TRUE)

# Take out population
pop_dfs <- mget(load(file.path(data.path, "Sverige_population_2019.Rdata")))
dat_pop <- pop_dfs$dat_pop
dat_pop_region <- pop_dfs$dat_pop_region
dat_pop_region_totalt <- pop_dfs$dat_pop_region_totalt

Region_population <- dat_pop_region_totalt %>% 
  mutate(ARegion = str_remove(ARegion, 
                              pattern = or("s ", " ") %R% "län" %R% END))

Region_population <- Region_population %>%
  bind_rows(Region_population %>% 
              summarize(Pop = sum(Pop)) %>% 
              mutate(ARegion = "Riket"))


#---------------------------------------------------------------------------------------------------
# Functions for analysis
#---------------------------------------------------------------------------------------------------

## Some important notes on notations

## p_symp := fraction reported cases. Denoted p_r in report. 
## p_asymp := fraction non-reported cases. Denoted p_o in report.
## p_lower_inf := how much lower infectivity non-reported cases have. Denoted q_o in report.
## t_b := 16 march 2020
## iter := In optimisation, how many iterations/guesses that are run.
## eta := rate of leaving incubation period. Denoted rho in report.



eta_value    <- 1 / 5.1
gammaD_value <- 1 / 5

## Tolerance for ode and optimisation. 
## This tolerance not always needed but better to be safe even though it takes a bit longer!
## For faster analysis use different tolerance or use default values.
Atol <- 1e-8
Rtol <- 1e-10

logit <- function(p) {
  log(p) / log(1 - p)
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

## Daily incidence reported cases and their dates
Incidence <- Stockholm_Data_10_april$Incidens
Datum     <- as.Date(Stockholm_Data_10_april$Datum)
Day       <- as.integer(Datum - as.Date("2019-12-31"))

Namedate <- seq.Date(as.Date("2020-01-01"), 
                     as.Date("2021-01-01"), 
                     by = "month")
dayatyear <- as.integer(Namedate - as.Date("2019-12-31"))

wfh_date <- as.Date("2020-03-16")
days_from_wfh_date <- as.integer(wfh_date - as.Date("2019-12-31"))

#' The time-dependent infectivity rate.
#' 
#' @param t
#' @param t_b
#' @param delta
#' @param epsilon
#' @param theta
beta <- function(t, t_b, delta, epsilon, theta){
  ((1 - delta) / (1 + exp(epsilon * (-(t - t_b)))) + delta) * theta
}

#' The time-dependent basic reproductive number.
#' 
#' @param t
#' @param t_b
#' @param delta
#' @param epsilon
#' @param theta
#' @param gamma
Basic_repr <- function(t, t_b, delta, epsilon, theta, gamma, p_symp, 
                       p_lower_inf) {
  a <- p_symp * beta(t, t_b, delta, epsilon, theta)
  b <- (1 - p_symp) * p_lower_inf * beta(t, t_b, delta, epsilon, theta)
  res <- (a + b) / gamma
  return(res)
}

#' Estimate free SEIR parameters for Stockholm data.
#' 
#' @param p_symp
#' @param p_lower_inf
#' @param gammaD
#' @param eta
#' @param iter
#' @param wfh_date
#' @param non_reported
Stockholm_SEIR <- function(
    p_symp = 0.5, 
    p_lower_inf = 0.5, 
    gammaD = gammaD_value, 
    eta = eta_value, 
    iter = 20,
    wfh_date = as.Date("2020-03-16"),
    non_reported = FALSE,
    verbose = TRUE) {
  
  ## Population size Stockholm
  N <- Region_population %>% filter(ARegion == "Stockholm") %>% pull(Pop)
  
  Opt_par_names <- c("logit_delta", "epsilon", "log_theta")
    
  ## Function to create guesses for the optimisation
  ## The range of the guesses can be changed, 
  ## these are good for the specific dates and parameter combinations of p_symp and p_lower_inf
  Guesses <- function() { 
    
    u_d <- logit(runif(1, 0.05, 0.6)) # guess for logit_delta 
    u_e <- runif(1, -0.6, 0)    # guess for epsilon
    u_t <- log(runif(1, 0, 15))      # guess for log_theta

    return(c(u_d, u_e, u_t))
  }
  
  ## The time-dependent infectivity rate 
  beta_peak_free <- function(t, delta, epsilon, theta){
    t_b <- as.integer(wfh_date - as.Date("2019-12-31"))
    res <- ((1 - delta) / (1 + exp(epsilon * (-(t - t_b)))) + delta) * theta 
    return(res)
  }
  
  ## The time-dependent basic reproductive number
  Basic_repr <- function(t, delta, epsilon, theta, gamma) {
    a <- p_symp * beta_peak_free(t, delta, epsilon, theta)
    b <- (1 - p_symp) * p_lower_inf * beta_peak_free(t, delta, epsilon, theta)
    res <- (a + b) / gamma
    return(res)
  }
 
  # The SEIR model. 
  # Note that the rate of going from latency to infectious is denoted eta here, 
  # rho in the report
  seir_model_asymptomatics <- function(time, state, parameters) {
    # S       <- state[1] # susceptibles
    # E       <- state[2] # latent/exposed but not infectious
    # I_symp  <- state[3] # infected who get reported
    # I_asymp <- state[4] # infected who remain non-reported
    # R       <- state[5] # recovered/immune
    par <- as.list(c(state, parameters))
    with(par, {
      dS <- -beta_peak_free(time, delta, epsilon, theta) * S * I_symp / N -
            p_lower_inf * beta_peak_free(time, delta, epsilon, theta) * S * I_asymp/N
      dE <- beta_peak_free(time, delta, epsilon, theta) * S * I_symp / N + 
            p_lower_inf * beta_peak_free(time, delta, epsilon, theta) * S * I_asymp/N - 
            eta * E
      dI_symp <- p_symp * eta * E - gammaD * I_symp
      dI_asymp <- (1 - p_symp) * eta * E - gammaD * I_asymp
      dR <- gammaD * (I_symp + I_asymp)
      dx <- c(dS, dE, dI_symp, dI_asymp, dR)
      list(dx)
    }
    )
  }
  
  # Assumption on initial number of infected:
  # In our main analysis we start with 1 infectious individual at 
  # t_0 = 2020-02-17. You can instead choose the commented one if you want to 
  # try with non-reported cases as well. 
  
  if (non_reported) {
    init <- c(S = N - Incidence[1] * (1 + (1 - p_symp) / p_symp), 
              E = 0, 
              I_symp = Incidence[1], 
              I_asymp = Incidence[1] * (1 - p_symp) / p_symp, 
              R = 0)
  } else {
    init <- c(S = N - Incidence[1],
              E = 0, 
              I_symp = Incidence[1], 
              I_asymp = 0, 
              R = 0)
  }
  
  model <- seir_model_asymptomatics
  
  RSS <- function(parameters) {
        
    names(parameters) <- Opt_par_names
    pars <- list(delta = expit(parameters["logit_delta"]),
                 epsilon = parameters["epsilon"],
                 theta = exp(parameters["log_theta"]))
    
    Dummy_infectivity <- beta_peak_free(t = c(0:700), 
                                        delta = pars$delta, 
                                        epsilon = pars$epsilon,  
                                        theta = pars$theta)
    # if the infectivity is negative, throw away guess
    if (min(Dummy_infectivity) < 0) {
      res <- 10^12
      #print("negative infectivity")
      return(res)
    } else {
      # choose tolerance atol and rtol
      #out <- ode(y = init, times = Day, func = model, parms = parameters)
      out <- ode(y = init, 
                 times = Day, 
                 func = model, 
                 parms = pars, 
                 atol = Atol, 
                 rtol = Rtol)
      
      fit_S <- out[ , "S"]
      fit_E <- out[ , "E"]
      fit_I_symp <- out[ , "I_symp"]
      fit_I_asymp <- out[ , "I_asymp"]
      
      fitted_incidence  <- p_symp * fit_E * eta
      
      # For transparency, the old incorrect incidence was expressed as:
      
      # fitted_incidence  <- beta_peak_free(
      #   out[,1], 
      #   delta = parameters[1], 
      #   epsilon = parameters[2], 
      #   theta = parameters[3]) * fit_S * fit_I_symp/N
      
      return(sum((Incidence - fitted_incidence)^2))
    }
  }
  print("Optimisation initialised")
  
  fitter <- function(x) {
    Guess <- Guesses()
    conl <- list(maxit = 1000, abstol = Atol, reltol = Rtol)
    Opt <- optim(Guess, RSS, control = list(conl), hessian = TRUE)

    while (Opt$convergence > 0) {
      Guess <- Guesses()
      Opt <- optim(Guess, RSS, control = list(conl), hessian = TRUE)
    }
    return(Opt)
  }
  
  fits <- furrr::future_map(1:iter, fitter, .progress = TRUE)
  best_index <- purrr::map(fits, ~ .x$value) %>% unlist() %>% which.min()
  Opt <- fits[[best_index]]
  
  Opt_par_transformed <- Opt$par
  names(Opt_par_transformed) <- Opt_par_names
  Opt_par <- c(delta = unname(expit(Opt_par_transformed["logit_delta"])),
               epsilon = unname(Opt_par_transformed["epsilon"]),
               theta = unname(exp(Opt_par_transformed["log_theta"])))
  return(list(args = list(p_symp = p_symp, 
                          p_lower_inf = p_lower_inf, 
                          gammaD = gammaD, 
                          eta = eta, 
                          iter = iter,
                          .wfh_date = .wfh_date,
                          non_reported = non_reported),
              Observed_incidence = Incidence, 
              Population_size = N, 
              Day = Day, 
              dayatyear = dayatyear, 
              Namedate = Namedate, 
              Optimisation = Opt, 
              Opt_par = Opt_par,
              all_fitted = fits,
              failed_solves = fails,
              Initial_values = init, 
              SEIR_model = model))
}




#---------------------------------------------------------------------------------------------------
# Analysis
#---------------------------------------------------------------------------------------------------


gammaD <- gammaD_value
eta <- eta_value


# Analysis  p_symp_use <- 0.0127
# Analysis  p_lower_inf_use <- 1, 0.55, 0.11
p_symp_use      <- 0.0127
p_asymp_use     <- 1 - p_symp_use
p_lower_inf_use <- 1


Est_par_model <- Estimate_function_Stockholm_only_local(p_symp = p_symp_use, 
                                                        p_lower_inf = p_lower_inf_use)


# Days of incidence
Day <- Est_par_model$Day

N <- Est_par_model$Population_size

dayatyear <- Est_par_model$dayatyear
Namedate  <- Est_par_model$Namedate 

#Observed incidence based on region
Observed_incidence <-  Est_par_model$Observed_incidence 
Est <- Est_par_model$Optimisation
Est
Opt_par <- Est_par_model$Opt_par

# functions based on model scenario
Basic_repr <- Est_par_model$Basic_reproduction
beta       <- Est_par_model$Infectivity
SEIR_model <- Est_par_model$SEIR_model

# initial values based on model scenario
init       <- Est_par_model$Initial_values


RSS_value <- Est$value

# Note: parameters have now been fitted on transformed scale
H         <- Est$hessian
sigest    <- sqrt(RSS_value/(length(Observed_incidence)-3))
NeginvH2  <- solve(1/(2*sigest^2)*H)
sdParams  <- sqrt(diag(NeginvH2))
names(sdParams) <- names(Opt_par)

options("scipen"=100, "digits"=4)
#default options("scipen"=0, "digits"=7)
RSS_value

#-------------------------------------------
# Look at the results, compare with prevalence 
# 27th March to 3rd April
#-------------------------------------------




t <- (Day[1]):(Day[length(Day)]+14+11) # time in days
t_date <- as.Date("2019-12-31") + t
fit <- ode(y = init, times = t, func = SEIR_model , parms = Opt_par)[ , ] %>%
  as_tibble() %>%
  rename(Day = time) %>%
  mutate(Date = as.Date("2019-12-31") + Day)
  
fit_S <- fit$S
fit_E <- fit$E
fit_I_symp <- fit$I_symp
fit_I_asymp <- fit$I_asymp
fit_R <- fit$R
fit_I <- fit_I_symp + fit_I_asymp
fit_I_E <- fit_E + fit_I
fit_cum_inf <- N - fit_S



## The mean prevalence same days as the Hälsorapport Stockholmsstudien (27th March to 3rd April)
Smittsamma <- fit_I_symp + fit_I_asymp #+ fit_E
SmittsammaF <-  Smittsamma[40:47]
mean(SmittsammaF/N)



## Look at the estimated reported cases and fitted

fitted_incidence <- p_symp_use * fit_E * eta
fitted_incidence_non_report  <- (1 - p_symp_use) * fit_E * eta

df_incidence <- Stockholm_Data_10_april %>% 
  rename(Date = Datum,
         Incidence = Incidens) %>%
  mutate(Type = "Observed") %>%
  bind_rows(tibble(Date = t_date,
                   Incidence = fitted_incidence,
                   Type = "Fitted (Reported)")) %>%
  bind_rows(tibble(Date = t_date,
                   Incidence = fitted_incidence_non_report,
                   Type = "Fitted (Unreported)")) %>%
  mutate(Day = as.integer(Date - min(Date) + 1)) %>%
  mutate(Type = factor(Type, levels = c("Observed", 
                                        "Fitted (Reported)",
                                        "Fitted (Unreported)")))

df_incidence %>%
  filter(Type != "Fitted (Unreported)") %>%
  ggplot() + 
  geom_line(aes(x = Date, y = Incidence, color = Type, size = Type)) + 
  geom_point(aes(x = Date, y = Incidence), 
             data = df_incidence %>% filter(Type == "Observed"),
             size = 2) +
  scale_color_manual(values = c("black", "red")) +
  scale_size_manual(values = c(0.5, 1.05)) +
  theme_minimal()


# Look at the estimated infectivity and basic reproductive number

df_R0 <- fit %>% select(Day, Date) %>%
  mutate(R0 = Basic_repr(Day, 
                         delta = Opt_par["delta"],
                         epsilon = Opt_par["epsilon"],  
                         theta = Opt_par["theta"],
                         gamma = gammaD))

df_infectivity <- fit %>% select(Day, Date) %>%
  mutate(Infectivity = beta(Day, 
                            delta = Opt_par["delta"], 
                            epsilon = Opt_par["epsilon"],  
                            theta = Opt_par["theta"])) 

plot_R0 <- df_R0 %>%
  ggplot() + 
  geom_line(aes(x = Date, y = R0), size = 1.05) +
  ylab("R0(t)") + 
  ggtitle("Estimated reproductive number") +
  theme_minimal()

plot_infectivity <- df_infectivity %>%
  ggplot() + 
  geom_line(aes(x = Date, y = Infectivity), size = 1.05) +
  ylab("Infectivity") + 
  ggtitle("Estimated infectivity") +
  theme_minimal()

plot_R0 + plot_infectivity + plot_layout(nrow = 2)


#-------------------------------------------------
# Calculate results to save in tables
#-------------------------------------------------


## To create bootstrap CI
CI_level_05 <- 0.05 / 2

delta_ci <- expit(qnorm(c(CI_level_05, 1 - CI_level_05), 
                        mean = Est$par[1], sd = sdParams[1], 
                        lower.tail = TRUE, log.p = FALSE))
epsilon_ci <- qnorm(c(CI_level_05, 1 - CI_level_05), 
                     mean = Est$par[2], sd = sdParams[2], 
                     lower.tail = TRUE, log.p = FALSE)
theta_ci <- exp(qnorm(c(CI_level_05, 1 - CI_level_05),
                      mean = Est$par[3], sd = sdParams[3], 
                      lower.tail = TRUE, log.p = FALSE))

n_sims <- 1000
par_sims <- MASS::mvrnorm(n = n_sims,
                          mu = Est$par,
                          Sigma = NeginvH2)
par_sims[, 1] <- expit(par_sims[, 1])
par_sims[, 3] <- exp(par_sims[, 3])
colnames(par_sims) <- names(Opt_par)

R0.v.Dag1 <- numeric(n_sims)
R0.v.DagSista <- numeric(n_sims)


sims_df_raw <- furrr::future_map_dfr(1:n_sims, function(n) {
  df1 <- tibble(
    iteration = n,
    Day = fit$Day,
    Date = fit$Date,
    R0 = map(fit$Day, 
             function(x) Basic_repr(x, 
                                    delta = par_sims[n, "delta"], 
                                    epsilon = par_sims[n, "epsilon"], 
                                    theta = par_sims[n, "theta"], 
                                    gamma = gammaD)) %>% unlist())
  df2 <- ode(y = init, 
             times = fit$Day, 
             func = SEIR_model, 
             parms = par_sims[n, ])[ , ] %>%
    as_tibble() %>%
    rename(Day = time) %>%
    mutate(Date = as.Date("2019-12-31") + Day)
  
  return(df1 %>% full_join(df2, by = c("Day", "Date")))
}, .progress = TRUE)

# Creds: https://tbradley1013.github.io/2018/10/01/calculating-quantiles-for-groups-with-dplyr-summarize-and-purrr-partial/
p <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
p_names <- map_chr(p, ~paste0(.x*100, "%"))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)
sims_df <- sims_df_raw %>% 
  pivot_longer(cols = c(R0, S, E, I_symp, I_asymp, R)) %>%
  select(-iteration) %>% 
  group_by(Day, Date, name) %>% 
  summarize_at(vars(everything()), p_funs)


sims_df %>%
  filter(name == "R0") %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = `50%`), size = 1.05) +
  geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.5) +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.25) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.125) +
  xlab("") + ylab("R0(t)") +
  ggtitle("Estimated reproductive number",
          glue("Black line shows median estimated R0(t),\n",
               "shaded regions show 50%, 90% and 95% uncertainty intervals ",
               "in order from darkest to lightest.")) +
  theme_minimal()

#27th March to 3rd April)
infectious_report <- sims_df_raw %>%
  filter(Date >= "2020-03-27", Date <= "2020-04-03") %>%
  group_by(iteration) %>%
  summarize(infectious = mean(I_symp + I_asymp) / N) %>%
  summarize(mean_infectious = mean(infectious),
            median_infectious = median(infectious),
            sd_infectious = sd(infectious),
            `2.5%` = quantile(infectious, 0.025),
            `97.5%` = quantile(infectious, 0.975))
  


#############################################
## save estimated parameters and their SE  ##
#############################################


res_param <- c(p_0                 = p_asymp_use,
               q_0                 = p_asymp_use,
               `27 mars - 3 april` = mean(SmittsammaF / N),
               "s.e."              = infectious_report$sd_infectious,
               RSS                 = RSS_value,
               delta               = unname(Opt_par["delta"]),
               "s.e."              = unname(sdParams["delta"]),
               epsilon             = unname(Opt_par["epsilon"]),
               "s.e."              = unname(sdParams["epsilon"]),
               theta               = unname(Opt_par["theta"]),
               "s.e."              = unname(Opt_par["theta"]))
# Round afterwards
res_param <- map2_dbl(res_param, c(3, 3, 5, 4, rep(3, 7)), round)

CIp       <- glue("[{round(delta_ci[1], digits = 3)}, {round(delta_ci[2], digits = 3)}]")
CIepsilon <- glue("[{round(epsilon_ci[1], digits = 3)}, {round(epsilon_ci[2], digits = 3)}]")
CItheta   <- glue("[{round(theta_ci[1], digits = 3)}, {round(theta_ci[2], digits = 3)}]")
CIsmittsam <- glue("[{round(infectious_report$`2.5%`, digits = 4)}, {round(infectious_report$`97.5%`, digits = 4)}]")

CI_param <- c("", "", CIsmittsam, "", "", CIp, "", CIepsilon, "", CItheta, "")

MAT_para <- matrix(c(res_param,CI_param), ncol = length(res_param), nrow = 2, byrow = TRUE)

df.res <- as.data.frame(MAT_para)
colnames(df.res) <- names(res_param)



XL_file_name <- paste(table.path,"/Res_para_p_non-reported_",p_asymp_use,"_infect_",p_lower_inf_use,".xlsx", sep ="")
write.xlsx(df.res, XL_file_name )




#############################################
## Save results with 31 days forecast ##
#############################################


t <- (Day[1]):(Day[length(Day)]+31) # time in days
fit         <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_par))
fit_S       <- fit[ , 2]
fit_E       <- fit[ , 3]
fit_I_symp  <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_R       <- fit[ , 6]
fit_I       <- fit_I_symp + fit_I_asymp
fit_I_E     <- fit_E + fit_I
fit_cum_inf <- N - fit_S




fitted_incidence  <- p_symp_use * fit_E * eta
fitted_incidence_non_report  <- (1 - p_symp_use) * fit_E * eta


Cum_Inf_inc <- data.frame(Cumulative = fit_cum_inf, 
                          Incidence_reported = fitted_incidence, Incidence_non_reported = fitted_incidence_non_report, Datum = as.Date(t, origin = "2019-12-31"))

fit_Cum_Inf_inc <- cbind(fit, Cum_Inf_inc)


file_name <- paste(table.path,"/Raw_data_fitted_model", "_para_p_asymp", p_asymp_use, "infect", p_lower_inf_use, ".txt", sep ="")

write.table(fit_Cum_Inf_inc, file=file_name, append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)




##################################
## Now calculate bootstrap CI's ##
##################################


fit_S.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_E.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_I_symp.v  <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_I_asymp.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
Fit_I.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_cum_inf.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))

fitted_incidence.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
effective_reprod.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))


for(i in 1:length(p.v)){
  Opt_parDummy = c(delta = p.v[i], epsilon = epsilon.v[i], theta = theta.v[i])
  
  fitDummy          <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_parDummy))
  fit_S.v[i,]       <- fitDummy[ , 2]
  fit_E.v[i,]       <- fitDummy[ , 3]
  fit_I_symp.v[i,]  <- fitDummy[ , 4]
  fit_I_asymp.v[i,] <- fitDummy[ , 5]
  Fit_I.v[i,]       <- fitDummy[ , 4] + fitDummy[ , 5]
  fit_cum_inf.v[i,] <- N - fitDummy[ , 2]
  
  fitted_incidence.v[i,]  <- p_symp_use *  fitDummy[ , 3] * eta
  #fitted_incidence.v[i,] <- beta(fitDummy[,1], delta = Opt_parDummy[1], epsilon = Opt_parDummy[2], theta = Opt_parDummy[3]) * fitDummy[ , 2] *  fitDummy[ , 4]/N 
  effective_reprod.v[i,] <- Basic_repr(fitDummy[,1], delta = Opt_parDummy[1], epsilon = Opt_parDummy[2], theta = Opt_parDummy[3], gamma = gammaD) * fitDummy[ , 2] /N
}


cum_inf_mean       <- apply(fit_cum_inf.v, MARGIN = 2, FUN = mean) 
cum_inf_median     <- apply(fit_cum_inf.v, MARGIN = 2, FUN = median) 
cum_inf_95_up_CRI  <- apply(fit_cum_inf.v, MARGIN = 2, FUN = CRI_95_up) 
cum_inf_95_low_CRI <- apply(fit_cum_inf.v, MARGIN = 2, FUN = CRI_95_low) 


fit_I_mean       <- apply(Fit_I.v, MARGIN = 2, FUN = mean) 
fit_I_median     <- apply(Fit_I.v, MARGIN = 2, FUN = median) 
fit_I_95_up_CRI  <- apply(Fit_I.v, MARGIN = 2, FUN = CRI_95_up) 
fit_I_95_low_CRI <- apply(Fit_I.v, MARGIN = 2, FUN = CRI_95_low) 

fitted_Incidence_mean       <- apply(fitted_incidence.v, MARGIN = 2, FUN = mean) 
fitted_Incidence_median     <- apply(fitted_incidence.v, MARGIN = 2, FUN = median) 
fitted_Incidence_95_up_CRI  <- apply(fitted_incidence.v, MARGIN = 2, FUN = CRI_95_up) 
fitted_Incidence_95_low_CRI <- apply(fitted_incidence.v, MARGIN = 2, FUN = CRI_95_low) 

effective_reprod_mean       <- apply(effective_reprod.v, MARGIN = 2, FUN = mean) 
effective_reprod_median     <- apply(effective_reprod.v, MARGIN = 2, FUN = median) 
effective_reprod_95_up_CRI  <- apply(effective_reprod.v, MARGIN = 2, FUN = CRI_95_up) 
effective_reprod_95_low_CRI <- apply(effective_reprod.v, MARGIN = 2, FUN = CRI_95_low) 




# Cummulative number of infected until 2020-04-11 and until 2020-05-01, with their 95% CI
# 2020-04-11 = dag 102
# 2020-05-01 = dag 122

as.Date(102, origin = "2019-12-31")
as.Date(122, origin = "2019-12-31")


fit_cum_inf[which(t==102)]

cum_inf_95_low_CRI[which(t==102)]
cum_inf_95_up_CRI[which(t==102)]

fit_cum_inf[which(t==102)]/N
cum_inf_95_low_CRI[which(t==102)]/N
cum_inf_95_up_CRI[which(t==102)]/N

fit_cum_inf[which(t==122)]
cum_inf_95_low_CRI[which(t==122)]
cum_inf_95_up_CRI[which(t==122)]

fit_cum_inf[which(t==122)]/N
cum_inf_95_low_CRI[which(t==122)]/N
cum_inf_95_up_CRI[which(t==122)]/N


as.Date(t[which(fit_I==max(fit_I))], origin = "2019-12-31" )
maxdagen <- as.Date(t[which(fit_I==max(fit_I))],origin = "2019-12-31" )

minDag <- as.Date(t[which(fit_I_95_low_CRI==max(fit_I_95_low_CRI))], origin = "2019-12-31" )
maxDag <- as.Date(t[which(fit_I_95_up_CRI==max(fit_I_95_up_CRI))], origin = "2019-12-31" )

as.Date(minDag, origin = "2019-12-31")
as.Date(maxDag, origin = "2019-12-31")




#############################################
## save estimated R0 and their uncertainty ##
#############################################

res_R0        <- round(c(p_asymp_use, p_lower_inf_use, R0_Mean[1], R0_Mean[length(Day)],effective_reprod_mean[length(Day)]), digits = 3)
names(res_R0) <- c("p_0", "q_0", "R0(start)", "R0(end)", "Re(end)")
CI_R01        <- paste("[",round(R0_low[1],digits=3), ", ", round(R0_high[1],digits=3),"]",sep="")
CI_ROend      <- paste("[",round(R0_low[2],digits=3), ", ", round(R0_high[2],digits=3),"]",sep="")
CI_Reend      <- paste("[",round(effective_reprod_95_low_CRI[length(Day)],digits=3), ", ", round(effective_reprod_95_up_CRI[length(Day)],digits=3),"]",sep="")
CIR0          <- c("", "", CI_R01, CI_ROend,CI_Reend)

MAT_R0    <- matrix(c(res_R0,CIR0), ncol =5, nrow=2, byrow=TRUE)
df.resR0  <- as.data.frame(MAT_R0)
colnames(df.resR0) <- names(res_R0)

XL_R0_file_name <- paste(table.path,"/Res_R0_p_non-reported_", p_asymp_use, "_infect_", p_lower_inf_use, ".xlsx", sep ="")
write.xlsx(df.resR0, XL_R0_file_name )



######################################
## save estimated days and their CI ##
######################################




fit_cum_low   <- cum_inf_95_low_CRI
fit_cum_high  <- cum_inf_95_up_CRI

res_days <-c(p_asymp_use, p_lower_inf_use, round(fit_cum_inf[which(t == 102)], digits=0), round(fit_cum_inf[which(t == 102)] / N, digits = 3) ,  
             round(fit_cum_inf[which(t == 122)],digits = 0), round(fit_cum_inf[which(t == 122)] / N, digits = 3))

res_maxprevalens <- round(max(fit_I), digits = 0)
res_days <- c(res_days, as.character(maxdagen), res_maxprevalens)

names(res_days)   <- c("p_0","q_0","2020-04-11", "", "2020-05-01","", "Peakdag", "Prevalens peakdag")
CIdag11mars       <- paste("[",round(fit_cum_low[which(t == 102)], digits=0), ", ", round(fit_cum_high[which(t == 102)], digits = 0), "]", sep="")
CIdag11marsProc   <- paste("[",round(fit_cum_low[which(t == 102)] / N, digits = 3), ", ",round(fit_cum_high[which(t == 102)]/N,digits=3),"]", sep="")
CIdag1maj         <- paste("[",round(fit_cum_low[which(t == 122)], digits=0), ", ", round(fit_cum_high[which(t == 122)],digits=0),"]", sep="")
CIdag1majProc     <- paste("[",round(fit_cum_low[which(t == 122)] / N, digits = 3),", ", round(fit_cum_high[which(t == 122)] / N , digits=3),"]", sep="")
CImaxdag          <- paste("[",as.character(as.Date(minDag, origin = "2019-12-31")), ", ", as.character(as.Date(maxDag, origin = "2019-12-31")), "]", sep="")
CIprevalensMaxdag <- paste("[",round(max(fit_I_95_low_CRI),digits = 0), ", ", round(max(fit_I_95_up_CRI), digits = 0), "]", sep="")



  
CI_dag  <- c("","", CIdag11mars, CIdag11marsProc, CIdag1maj, CIdag1majProc, CImaxdag, CIprevalensMaxdag)
MAT_dag <- matrix(c(res_days, CI_dag), ncol = 8, nrow = 2, byrow = TRUE)
df.dag  <- as.data.frame(MAT_dag)

colnames(df.dag) <- names(res_days)


XL_file_name <- paste(table.path,"/Res_dagar_p_non-reported_", p_asymp_use, "_infect_", p_lower_inf_use, ".xlsx", sep ="")
write.xlsx(df.dag, XL_file_name)


#----------------------------------------------------
# Save figures
#----------------------------------------------------


t   <- (Day[1]):(Day[length(Day)]+14) # time in days
fit <- data.frame(ode(y = init, times = t, func = SEIR_model, parms = Opt_par))

fit_S       <- fit[ , 2]
fit_E       <- fit[ , 3]
fit_I_symp  <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_I       <- fit_I_symp + fit_I_asymp
fit_I_E     <- fit_E + fit_I
fit_cum_inf <- N - fit_S



#fitted_incidence  <- beta(fit[,1], delta = Opt_par[1], epsilon = Opt_par[2], theta = Opt_par[3]) * fit_S *  fit_I_symp/N 
fitted_incidence  <- p_symp_use * fit_E * eta



fitted_incidence_low <- fitted_Incidence_95_low_CRI[1:length(t)]
fitted_incidence_high <- fitted_Incidence_95_up_CRI[1:length(t)]

fitted_I_high <- fit_I_95_up_CRI[1:length(t)]
fitted_I_low <- fit_I_95_low_CRI[1:length(t)]






NameNumber <- paste("/Incidence_number_infected_14Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")



pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



# plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
#      xlab="",ylim=c(0,roundUp(max(fitted_incidence))),xaxt='n')

plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)


# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


## next plot

plot(fit$time, fit[ , 4] + fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)

axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()







NameNumberPNG <- paste("/Incidence_number_infected_14Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)


# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)



axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)



plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)



axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()





NameNumber <- paste("/Incidence_number_infected_14Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")

pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')



grid(nx=NA, ny=NULL)


# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()




NameNumberPNG <- paste("/Incidence_number_infected_14Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')

grid(nx=NA, ny=NULL)

# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_april <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14)
NameDateMarchApril <- as.Date(dayatyear_march_april,origin ="2019-12-31")

abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)

## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')



grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_april,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_april, label = NameDateMarchApril,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()









#----------------------------------------------------
# Save figures 31 day forecast
#----------------------------------------------------


t   <- (Day[1]):(Day[length(Day)] + 31) # time in days
fit <- data.frame(ode(y = init, times = t, func = SEIR_model, parms = Opt_par))

fit_S       <- fit[ , 2]
fit_E       <- fit[ , 3]
fit_I_symp  <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_I       <- fit_I_symp + fit_I_asymp
fit_I_E     <- fit_E + fit_I
fit_cum_inf <- N - fit_S



#fitted_incidence  <- beta(fit[,1], delta = Opt_par[1], epsilon = Opt_par[2], theta = Opt_par[3]) * fit_S *  fit_I_symp/N 
fitted_incidence  <- p_symp_use * fit_E * eta

fitted_incidence_low <- fitted_Incidence_95_low_CRI[1:length(t)]
fitted_incidence_high <- fitted_Incidence_95_up_CRI[1:length(t)]

fitted_I_high <- fit_I_95_up_CRI[1:length(t)]
fitted_I_low <- fit_I_95_low_CRI[1:length(t)]





# 1 mars, 15 mars, 1 april, 15 april
dayatyear_march_may <- c(Day[1], 61, 61 + 14, 61 + 14 + 17, 61 + 14 + 17 + 14 , 61 + 14 + 17 + 14 + 16, 61 + 14 + 17 + 14 + 16 + 14)
NameDateMarchMay <- as.Date(dayatyear_march_may,origin ="2019-12-31")


NameNumber <- paste("/Incidence_number_infected_31Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")



pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



# plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
#      xlab="",ylim=c(0,roundUp(max(fitted_incidence))),xaxt='n')

plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)


abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


## next plot

plot(fit$time, fit[ , 4] + fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)

axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()







NameNumberPNG <- paste("/Incidence_number_infected_31Days_CI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(fitted_incidence_high))),xaxt='n')


grid(nx=NA, ny=NULL)



abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

lines(fit$time, fitted_incidence_low, lty=2)
lines(fit$time, fitted_incidence_high, lty=2)



axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)



plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fitted_I_high)) ) , type="l", col="red", 
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)

lines(fit$time, fitted_I_high, lty=2)
lines(fit$time, fitted_I_low, lty=2)



axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()





NameNumber <- paste("/Incidence_number_infected_31Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".pdf",sep="")

pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 

## Prognosis until May



plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')



grid(nx=NA, ny=NULL)


abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)

axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')

grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()




NameNumberPNG <- paste("/Incidence_number_infected_31Days_noCI_non-reported_", round((1-p_symp_use)*100,digits=1),"perc_less_inf_",p_lower_inf_use*100,"perc",".png",sep="")

#width=9*1.5, height=4.5*1.3 

png(paste(figure.path, NameNumberPNG, sep=""),width = 9*1.5*75, height =  4.5*1.3*75 ) 
par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 


plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
     xlab="",ylim=c(0,roundUp(max(c(fitted_incidence,Observed_incidence)))),xaxt='n')

grid(nx=NA, ny=NULL)


abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))

lines(fit$time, fitted_incidence,col="red",lwd=2)
points(Day, Observed_incidence)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)

## next plot

plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUp(max(fit[ , 4]+fit[ , 5])) ) , type="l", col="red",
     ylab = "Number of infected",lwd=2,
     xlab="",xaxt='n')



grid(nx=NA, ny=NULL)
abline(v=dayatyear_march_may,col = "lightgray", lty = "dotted", lwd = par("lwd"))
lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)


axis(side = 1, at = dayatyear_march_may, label = NameDateMarchMay,las=2)


MainTitle <- paste("Fitted SEIR model covid-19 Stockholm \n Unreported cases ",round((1-p_symp_use)*100,digits=1),
                   "% of infected and their infectivity ",p_lower_inf_use*100,"% compared to reported cases",sep="")

title(MainTitle, outer = TRUE, line = -2.5)


title("Estimated and observed number of reported cases", outer = TRUE, line = -5, adj = 0.11)
title("Number of infected at the same time", outer = TRUE, line = -5, adj = 0.825)


dev.off()








