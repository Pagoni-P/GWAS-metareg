## PART 3: Calculate performance measures
##
## Pagoni P, Higgins JPT, Lawlor DA, Stergiakouli E, Warrington NM, Morris TT, Tilling K.
## Meta-regression of genome-wide association studies to estimate age-varying genetic effects.
## 
## Panagiota Pagoni
## 31/10/2023

#install.packages("haven")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("knitr")
#install.packages("tidyverse")
library("haven")
library("tidyr")
library("dplyr")
library("knitr")
library("tidyverse")

N <- 1000

age <- list( 0 , 25, 50, 75, 100)
  
for ( age_overlap in age ){
  
  for (s in 1:5) {
    
    # load dataset including 10,000 iterations
    
    data <- read_dta( paste( paste ( paste (  "directory/res_all_sc", s, sep = ""),
                                     "_final_", age_overlap,sep=""),
                             N,"1000.dta",sep ="_"))

    na <- sum( is.na(data) )
    
## Set simulated values for each scenario 
    
    if ( s == 1 | s==3 ){
      
      sim_snp <- 1.5 # simulated value for beta_snp 
      sim_int <- 0 # simulated value for beta_int 
      sim_int2 <- 0 # simulated value for beta_int2
      
    }else if ( s ==2 | s==4 ){
      
      sim_snp <- 1.5 + 0.02*10 # simulated value for beta_snp 
      sim_int <- 0.02 # simulated value for beta_int 
      sim_int2 <- 0 # simulated value for beta_int2
      
    }else{
      
      sim_snp <- 1.5 + 0.02*10 + 0.001*10*10 # simulated value for beta_snp 
      sim_int <- 0.02 + 2*0.001*10 # simulated value for beta_int 
      sim_int2 <- 0.001 # simulated value for beta_int2
    
    }
    
## 95% CI of main genetic effect
    
    # Fixed - effect meta - analysis
    
    data$low_ci_beta_snp_fixed <- (data$beta_snp_fixed_metan - 1.96*data$se_snp_fixed_metan) # low CI for the main genetic effect
    data$up_ci_beta_snp_fixed <- (data$beta_snp_fixed_metan + 1.96*data$se_snp_fixed_metan) # high CI for the main genetic effect 
    
    # Random - effects meta - analysis 
    
    data$low_ci_beta_snp_random <- (data$beta_snp_random_metan - 1.96*data$se_snp_random_metan) # low CI for the main genetic effect
    data$up_ci_beta_snp_random <- (data$beta_snp_random_metan + 1.96*data$se_snp_random_metan) # high CI for the main genetic effect 
    
    # Meta - regression ( linear term )
    
    data$low_ci_beta_snp_metareg_A <- (data$beta_snp_metareg_A - 1.96*data$se_beta_snp_metareg_A) # low CI for the main genetic effect
    data$up_ci_beta_snp_metareg_A <- (data$beta_snp_metareg_A + 1.96*data$se_beta_snp_metareg_A) # high CI for the main genetic effect 
    
    data$low_ci_beta_int_metareg_A <- (data$beta_int_metareg_A - 1.96*data$se_int_metareg_A) # low CI for the age dependent effect
    data$up_ci_beta_int_metareg_A <- (data$beta_int_metareg_A + 1.96*data$se_int_metareg_A) # high CI for the age dependent effect 
    
    # Meta - regression ( non linear term )
    
    data$low_ci_beta_snp_metareg_B <- (data$beta_snp_metareg_B - 1.96*data$se_beta_snp_metareg_B) # low CI for the main genetic effect
    data$up_ci_beta_snp_metareg_B <- (data$beta_snp_metareg_B + 1.96*data$se_beta_snp_metareg_B) # high CI for the main genetic effect 
    
    data$low_ci_beta_int_metareg_B <- (data$beta_int_metareg_B - 1.96*data$se_int_metareg_B) # low CI for the age dependent effect int
    data$up_ci_beta_int_metareg_B <- (data$beta_int_metareg_B + 1.96*data$se_int_metareg_B) # high CI for the age dependent effect int
    
    data$low_ci_beta_int2_metareg_B <- (data$beta_int2_metareg_B - 1.96*data$se_int2_metareg_B) # low CI for the age dependent effect int2
    data$up_ci_beta_int2_metareg_B <- (data$beta_int2_metareg_B + 1.96*data$se_int2_metareg_B) # high CI for the age dependent effect int2
    
    
## Create a dataframe to save results per scenario 
    
    simsum_res <- data.frame ( matrix ( NA, nrow = 4 , ncol = 59 ))
    
    colnames(simsum_res) <- c( 'scenario', 'age_overlap','N','method','wmean_all_ages',
                               'sim_snp','sim_int','sim_int2',
                               'mean_beta_snp', 'mean_se_beta_snp', 'low_beta_snp','high_beta_snp','empSE_beta_snp', 
                               'bias_beta_snp', 'bias_mcse_beta_snp', 'low_bias_beta_snp', 'high_bias_beta_snp',
                               'coverage_beta_snp', 'coverage_mcse_beta_snp', 'low_coverage_beta_snp', 'high_coverage_beta_snp',
                               'power_snp', 'power_mcse_beta_snp', 'low_power_beta_snp', 'high_power_beta_snp',
                               
                               
                               'mean_beta_int', 'mean_se_beta_int', 'low_beta_int','high_beta_int','empSE_beta_int', 
                               'bias_beta_int', 'bias_mcse_beta_int', 'low_bias_beta_int', 'high_bias_beta_int',
                               'coverage_beta_int', 'coverage_mcse_beta_int', 'low_coverage_beta_int', 'high_coverage_beta_int',   
                               'power_beta_int', 'power_mcse_beta_int', 'low_power_beta_int', 'high_power_beta_int',   
                               
                               'mean_beta_int2', 'mean_se_beta_int2', 'low_beta_int2','high_beta_int2','empSE_beta_int2', 
                               'bias_beta_int2', 'bias_mcse_beta_int2', 'low_bias_beta_int2', 'high_bias_beta_int2',
                               'coverage_beta_int2', 'coverage_mcse_beta_int2', 'low_coverage_beta_int2', 'high_coverage_beta_int2',
                               'power_beta_int2', 'power_mcse_beta_int2', 'low_power_beta_int2', 'high_power_beta_int2')

    
        
    simsum_res$scenario <- s
    simsum_res$N  <- N
    simsum_res$age_overlap <- age_overlap
    simsum_res$wmean_all_ages <- unique(data$wmean_all_ages)

    simsum_res$sim_snp <- sim_snp
    simsum_res$sim_int <- sim_int
    simsum_res$sim_int2 <- sim_int2
  
    # Fixed - effect meta - analysis
    
    simsum_res$method[[1]] <- 1
    simsum_res$mean_beta_snp[[1]] <- mean ( data$beta_snp_fixed_metan , na.rm = T  ) # mean beta of main genetic effect across simulations
    simsum_res$mean_se_beta_snp[[1]] <- sqrt ( ( 1 / 1000 ) * sum ( ( data$se_snp_fixed_metan )^2 ) ) # se of main genetic effect across simulations
    simsum_res$low_beta_snp[[1]] <- simsum_res$mean_beta_snp[[1]] - 1.96*simsum_res$mean_se_beta_snp[[1]] # low CI band of main genetic effect 
    simsum_res$high_beta_snp[[1]] <- simsum_res$mean_beta_snp[[1]] + 1.96*simsum_res$mean_se_beta_snp[[1]] # hight CI band of main genetic effect
    
    simsum_res$empSE_beta_snp[[1]] <- sd ( data$beta_snp_fixed_metan , na.rm = T ) # empirical SE of main genetic
    
    simsum_res$bias_beta_snp[[1]] <- simsum_res$mean_beta_snp[[1]] - sim_snp # Bias estimate for mean main genetic effect
    simsum_res$bias_mcse_beta_snp[[1]] <-   sqrt( ( 1 / (1000* ( 1000 - 1 ))) * sum( ( data$beta_snp_fixed_metan - simsum_res$mean_beta_snp[[1]] ) ^2 ) ) # Monte Carlo SE of bias
    
    simsum_res$low_bias_beta_snp[[1]] <- simsum_res$bias_beta_snp[[1]] - 1.96*simsum_res$bias_mcse_beta_snp[[1]] # low Ci band of bias 
    simsum_res$high_bias_beta_snp[[1]] <- simsum_res$bias_beta_snp[[1]] + 1.96*simsum_res$bias_mcse_beta_snp[[1]] # high Ci band of bias 
    
    simsum_res$coverage_beta_snp[[1]] <- ( sum( ifelse( (data$up_ci_beta_snp_fixed >= sim_snp & data$low_ci_beta_snp_fixed <= sim_snp), 1 ,0 ) ) / 1000 ) * 100
    simsum_res$coverage_mcse_beta_snp[[1]] <- sqrt ( (simsum_res$coverage_beta_snp[[1]] * 0.01 ) * ( 1 - simsum_res$coverage_beta_snp[[1]] * 0.01 ) / 1000 ) # Monte Carlo SE of coverage
    
    simsum_res$low_coverage_beta_snp[[1]] <- ( simsum_res$coverage_beta_snp[[1]]*0.01 - 1.96*simsum_res$coverage_mcse_beta_snp[[1]] ) * 100
    simsum_res$high_coverage_beta_snp[[1]] <- ( simsum_res$coverage_beta_snp[[1]]*0.01 + 1.96*simsum_res$coverage_mcse_beta_snp[[1]] ) * 100
    
    # Random - effects meta - analysis
    
    simsum_res$method[[2]] <- 2
    simsum_res$mean_beta_snp[[2]] <- mean ( data$beta_snp_random_metan , na.rm = T  ) # mean beta of main genetic effect across simulations
    simsum_res$mean_se_beta_snp[[2]] <- sqrt ( ( 1 / 1000 ) * sum ( ( data$se_snp_random_metan )^2 ) ) # se of main genetic effect across simulations
    simsum_res$low_beta_snp[[2]] <- simsum_res$mean_beta_snp[[2]] - 1.96*simsum_res$mean_se_beta_snp[[2]] # low CI band of main genetic effect 
    simsum_res$high_beta_snp[[2]] <- simsum_res$mean_beta_snp[[2]] + 1.96*simsum_res$mean_se_beta_snp[[2]] # hight CI band of main genetic effect
    
    simsum_res$empSE_beta_snp[[2]] <- sd ( data$beta_snp_random_metan , na.rm = T ) # empirical SE of main genetic
    
    simsum_res$bias_beta_snp[[2]] <- simsum_res$mean_beta_snp[[2]] - sim_snp # Bias estimate for mean main genetic effect
    simsum_res$bias_mcse_beta_snp[[2]] <-   sqrt( ( 1 / (1000* ( 1000 - 1 ))) * sum( ( data$beta_snp_random_metan - simsum_res$mean_beta_snp[[2]] ) ^2 ) ) # Monte Carlo SE of bias
    
    simsum_res$low_bias_beta_snp[[2]] <- simsum_res$bias_beta_snp[[2]] - 1.96*simsum_res$bias_mcse_beta_snp[[2]] # low Ci band of bias 
    simsum_res$high_bias_beta_snp[[2]] <- simsum_res$bias_beta_snp[[2]] + 1.96*simsum_res$bias_mcse_beta_snp[[2]] # high Ci band of bias 
    
    simsum_res$coverage_beta_snp[[2]] <- ( sum( ifelse( (data$up_ci_beta_snp_random >= sim_snp & data$low_ci_beta_snp_random <= sim_snp), 1 ,0 ) ) / 1000 ) * 100
    simsum_res$coverage_mcse_beta_snp[[2]] <- sqrt ( (simsum_res$coverage_beta_snp[[2]] * 0.01 ) * ( 1 - simsum_res$coverage_beta_snp[[2]] * 0.01 ) / 1000 ) # Monte Carlo SE of coverage
    
    simsum_res$low_coverage_beta_snp[[2]] <- ( simsum_res$coverage_beta_snp[[2]]*0.01 - 1.96*simsum_res$coverage_mcse_beta_snp[[2]] ) * 100
    simsum_res$high_coverage_beta_snp[[2]] <- ( simsum_res$coverage_beta_snp[[2]]*0.01 + 1.96*simsum_res$coverage_mcse_beta_snp[[2]] ) * 100

    
# Meta - regression ( linear term ) (A)
    
    # Main genetic effect
    
    simsum_res$method[[3]] <- 3
    simsum_res$mean_beta_snp[[3]] <- mean ( data$beta_snp_metareg_A , na.rm = T  ) # mean beta of main genetic effect across simulations
    simsum_res$mean_se_beta_snp[[3]] <- sqrt ( ( 1 / 1000 ) * sum ( ( data$se_beta_snp_metareg_A )^2 ) ) # se of main genetic effect across simulations
    simsum_res$low_beta_snp[[3]] <- simsum_res$mean_beta_snp[[3]] - 1.96*simsum_res$mean_se_beta_snp[[3]] # low CI band of main genetic effect 
    simsum_res$high_beta_snp[[3]] <- simsum_res$mean_beta_snp[[3]] + 1.96*simsum_res$mean_se_beta_snp[[3]] # hight CI band of main genetic effect
    
    simsum_res$empSE_beta_snp[[3]] <- sd ( data$beta_snp_metareg_A , na.rm = T ) # empirical SE of main genetic
    
    simsum_res$bias_beta_snp[[3]] <- simsum_res$mean_beta_snp[[3]] - sim_snp # Bias estimate for mean main genetic effect
    simsum_res$bias_mcse_beta_snp[[3]] <-   sqrt( ( 1 / (1000* ( 1000 - 1 ))) * sum( ( data$beta_snp_metareg_A - simsum_res$mean_beta_snp[[3]] ) ^2 ) ) # Monte Carlo SE of bias
    
    simsum_res$low_bias_beta_snp[[3]] <- simsum_res$bias_beta_snp[[3]] - 1.96*simsum_res$bias_mcse_beta_snp[[3]] # low Ci band of bias 
    simsum_res$high_bias_beta_snp[[3]] <- simsum_res$bias_beta_snp[[3]] + 1.96*simsum_res$bias_mcse_beta_snp[[3]] # high Ci band of bias 
    
    simsum_res$coverage_beta_snp[[3]] <- ( sum( ifelse( (data$up_ci_beta_snp_metareg_A >= sim_snp & data$low_ci_beta_snp_metareg_A <= sim_snp), 1 ,0 ) ) / 1000 ) * 100
    simsum_res$coverage_mcse_beta_snp[[3]] <- sqrt ( (simsum_res$coverage_beta_snp[[3]] * 0.01 ) * ( 1 - simsum_res$coverage_beta_snp[[3]] * 0.01 ) / 1000 ) # Monte Carlo SE of coverage
    
    simsum_res$low_coverage_beta_snp[[3]] <- ( simsum_res$coverage_beta_snp[[3]]*0.01 - 1.96*simsum_res$coverage_mcse_beta_snp[[3]] ) * 100
    simsum_res$high_coverage_beta_snp[[3]] <- ( simsum_res$coverage_beta_snp[[3]]*0.01 + 1.96*simsum_res$coverage_mcse_beta_snp[[3]] ) * 100
    
    # Age dependent genetic effect (int)
    
    simsum_res$method[[3]] <- 3
    simsum_res$mean_beta_int[[3]] <- mean ( data$beta_int_metareg_A , na.rm = T  ) # mean beta of main genetic effect across simulations
    simsum_res$mean_se_beta_int[[3]] <- sqrt ( ( 1 / 1000 ) * sum ( ( data$se_int_metareg_A )^2 ) ) # se of main genetic effect across simulations
    simsum_res$low_beta_int[[3]] <- simsum_res$mean_beta_int[[3]] - 1.96*simsum_res$mean_se_beta_int[[3]] # low CI band of main genetic effect 
    simsum_res$high_beta_int[[3]] <- simsum_res$mean_beta_int[[3]] + 1.96*simsum_res$mean_se_beta_int[[3]] # hight CI band of main genetic effect
    
    simsum_res$empSE_beta_int[[3]] <- sd ( data$beta_int_metareg_A , na.rm = T ) # empirical SE of main genetic
    
    simsum_res$bias_beta_int[[3]] <- simsum_res$mean_beta_int[[3]] - sim_int # Bias estimate for mean main genetic effect
    simsum_res$bias_mcse_beta_int[[3]] <-   sqrt( ( 1 / (1000* ( 1000 - 1 ))) * sum( ( data$beta_int_metareg_A - simsum_res$mean_beta_int[[3]] ) ^2 ) ) # Monte Carlo SE of bias
    
    simsum_res$low_bias_beta_int[[3]] <- simsum_res$bias_beta_int[[3]] - 1.96*simsum_res$bias_mcse_beta_int[[3]] # low Ci band of bias 
    simsum_res$high_bias_beta_int[[3]] <- simsum_res$bias_beta_int[[3]] + 1.96*simsum_res$bias_mcse_beta_int[[3]] # high Ci band of bias 
    
    simsum_res$coverage_beta_int[[3]] <- ( sum( ifelse( (data$up_ci_beta_int_metareg_A >= sim_int & data$low_ci_beta_int_metareg_A <= sim_int), 1 ,0 ) ) / 1000 ) * 100
    simsum_res$coverage_mcse_beta_int[[3]] <- sqrt ( (simsum_res$coverage_beta_int[[3]] * 0.01 ) * ( 1 - simsum_res$coverage_beta_int[[3]] * 0.01 ) / 1000 ) # Monte Carlo SE of coverage
    
    simsum_res$low_coverage_beta_int[[3]] <- ( simsum_res$coverage_beta_int[[3]]*0.01 - 1.96*simsum_res$coverage_mcse_beta_int[[3]] ) * 100
    simsum_res$high_coverage_beta_int[[3]] <- ( simsum_res$coverage_beta_int[[3]]*0.01 + 1.96*simsum_res$coverage_mcse_beta_int[[3]] ) * 100
    
    simsum_res$power_beta_int[[3]] <- ( sum( ifelse( abs(data$beta_int_metareg_A) >= (1.96*data$se_int_metareg_A), 1 ,0 ) ) / 1000 ) *100 
    simsum_res$power_mcse_beta_int[[3]] <- sqrt ( (simsum_res$power_beta_int[[3]]*0.01) * ( 1 - simsum_res$power_beta_int[[3]]*0.01 ) / 1000 ) # Monte Carlo SE of power

    # Meta - regression ( non linear term )(B)
    
    # Main genetic effect
    
    simsum_res$method[[4]] <- 4
    simsum_res$mean_beta_snp[[4]] <- mean ( data$beta_snp_metareg_B , na.rm = T  ) # mean beta of main genetic effect across simulations
    simsum_res$mean_se_beta_snp[[4]] <- sqrt ( ( 1 / 1000 ) * sum ( ( data$se_beta_snp_metareg_B )^2 ) ) # se of main genetic effect across simulations
    simsum_res$low_beta_snp[[4]] <- simsum_res$mean_beta_snp[[4]] - 1.96*simsum_res$mean_se_beta_snp[[4]] # low CI band of main genetic effect 
    simsum_res$high_beta_snp[[4]] <- simsum_res$mean_beta_snp[[4]] + 1.96*simsum_res$mean_se_beta_snp[[4]] # hight CI band of main genetic effect
    
    simsum_res$empSE_beta_snp[[4]] <- sd ( data$beta_snp_metareg_B , na.rm = T ) # empirical SE of main genetic
    
    simsum_res$bias_beta_snp[[4]] <- simsum_res$mean_beta_snp[[4]] - sim_snp # Bias estimate for mean main genetic effect
    simsum_res$bias_mcse_beta_snp[[4]] <-   sqrt( ( 1 / (1000* ( 1000 - 1 ))) * sum( ( data$beta_snp_metareg_B - simsum_res$mean_beta_snp[[4]] ) ^2 ) ) # Monte Carlo SE of bias
    
    simsum_res$low_bias_beta_snp[[4]] <- simsum_res$bias_beta_snp[[4]] - 1.96*simsum_res$bias_mcse_beta_snp[[4]] # low Ci band of bias 
    simsum_res$high_bias_beta_snp[[4]] <- simsum_res$bias_beta_snp[[4]] + 1.96*simsum_res$bias_mcse_beta_snp[[4]] # high Ci band of bias 
    
    simsum_res$coverage_beta_snp[[4]] <- ( sum( ifelse( (data$up_ci_beta_snp_metareg_B >= sim_snp & data$low_ci_beta_snp_metareg_B <= sim_snp), 1 ,0 ) ) / 1000 ) * 100
    simsum_res$coverage_mcse_beta_snp[[4]] <- sqrt ( (simsum_res$coverage_beta_snp[[4]] * 0.01 ) * ( 1 - simsum_res$coverage_beta_snp[[4]] * 0.01 ) / 1000 ) # Monte Carlo SE of coverage
    
    simsum_res$low_coverage_beta_snp[[4]] <- ( simsum_res$coverage_beta_snp[[4]]*0.01 - 1.96*simsum_res$coverage_mcse_beta_snp[[4]] ) * 100
    simsum_res$high_coverage_beta_snp[[4]] <- ( simsum_res$coverage_beta_snp[[4]]*0.01 + 1.96*simsum_res$coverage_mcse_beta_snp[[4]] ) * 100
    
    # Age dependent genetic effect (int)
    
    simsum_res$method[[4]] <- 4
    simsum_res$mean_beta_int[[4]] <- mean ( data$beta_int_metareg_B , na.rm = T  ) # mean beta of main genetic effect across simulations
    simsum_res$mean_se_beta_int[[4]] <- sqrt ( ( 1 / 1000 ) * sum ( ( data$se_int_metareg_B )^2 ) ) # se of main genetic effect across simulations
    simsum_res$low_beta_int[[4]] <- simsum_res$mean_beta_int[[4]] - 1.96*simsum_res$mean_se_beta_int[[4]] # low CI band of main genetic effect 
    simsum_res$high_beta_int[[4]] <- simsum_res$mean_beta_int[[4]] + 1.96*simsum_res$mean_se_beta_int[[4]] # hight CI band of main genetic effect
    
    simsum_res$empSE_beta_int[[4]] <- sd ( data$beta_int_metareg_B , na.rm = T ) # empirical SE of main genetic
    
    simsum_res$bias_beta_int[[4]] <- simsum_res$mean_beta_int[[4]] - sim_int # Bias estimate for mean main genetic effect
    simsum_res$bias_mcse_beta_int[[4]] <-   sqrt( ( 1 / (1000* ( 1000 - 1 ))) * sum( ( data$beta_int_metareg_B - simsum_res$mean_beta_int[[4]] ) ^2 ) ) # Monte Carlo SE of bias
    
    simsum_res$low_bias_beta_int[[4]] <- simsum_res$bias_beta_int[[4]] - 1.96*simsum_res$bias_mcse_beta_int[[4]] # low Ci band of bias 
    simsum_res$high_bias_beta_int[[4]] <- simsum_res$bias_beta_int[[4]] + 1.96*simsum_res$bias_mcse_beta_int[[4]] # high Ci band of bias 
    
    simsum_res$coverage_beta_int[[4]] <- ( sum( ifelse( (data$up_ci_beta_int_metareg_B >= sim_int & data$low_ci_beta_int_metareg_B <= sim_int), 1 ,0 ) ) / 1000 ) * 100
    simsum_res$coverage_mcse_beta_int[[4]] <- sqrt ( (simsum_res$coverage_beta_int[[4]] * 0.01 ) * ( 1 - simsum_res$coverage_beta_int[[4]] * 0.01 ) / 1000 ) # Monte Carlo SE of coverage
    
    simsum_res$low_coverage_beta_int[[4]] <- ( simsum_res$coverage_beta_int[[4]]*0.01 - 1.96*simsum_res$coverage_mcse_beta_int[[4]] ) * 100
    simsum_res$high_coverage_beta_int[[4]] <- ( simsum_res$coverage_beta_int[[4]]*0.01 + 1.96*simsum_res$coverage_mcse_beta_int[[4]] ) * 100
    
    # Age dependent genetic effect (int2)
    
    simsum_res$method[[4]] <- 4
    simsum_res$mean_beta_int2[[4]] <- mean ( data$beta_int2_metareg_B , na.rm = T  ) # mean beta of main genetic effect across simulations
    simsum_res$mean_se_beta_int2[[4]] <- sqrt ( ( 1 / 1000 ) * sum ( ( data$se_int2_metareg_B )^2 ) ) # se of main genetic effect across simulations
    simsum_res$low_beta_int2[[4]] <- simsum_res$mean_beta_int2[[4]] - 1.96*simsum_res$mean_se_beta_int2[[4]] # low CI band of main genetic effect 
    simsum_res$high_beta_int2[[4]] <- simsum_res$mean_beta_int2[[4]] + 1.96*simsum_res$mean_se_beta_int2[[4]] # hight CI band of main genetic effect
    
    simsum_res$empSE_beta_int2[[4]] <- sd ( data$beta_int2_metareg_B , na.rm = T ) # empirical SE of main genetic
    
    simsum_res$bias_beta_int2[[4]] <- simsum_res$mean_beta_int2[[4]] - sim_int2 # Bias estimate for mean main genetic effect
    simsum_res$bias_mcse_beta_int2[[4]] <-   sqrt( ( 1 / (1000* ( 1000 - 1 ))) * sum( ( data$beta_int2_metareg_B - simsum_res$mean_beta_int2[[4]] ) ^2 ) ) # Monte Carlo SE of bias
    
    simsum_res$low_bias_beta_int2[[4]] <- simsum_res$bias_beta_int2[[4]] - 1.96*simsum_res$bias_mcse_beta_int2[[4]] # low Ci band of bias 
    simsum_res$high_bias_beta_int2[[4]] <- simsum_res$bias_beta_int2[[4]] + 1.96*simsum_res$bias_mcse_beta_int2[[4]] # high Ci band of bias 
    
    simsum_res$coverage_beta_int2[[4]] <- ( sum( ifelse( (data$up_ci_beta_int2_metareg_B >= sim_int2 & data$low_ci_beta_int2_metareg_B <= sim_int2), 1 ,0 ) ) / 1000 ) * 100
    simsum_res$coverage_mcse_beta_int2[[4]] <- sqrt ( (simsum_res$coverage_beta_int2[[4]] * 0.01 ) * ( 1 - simsum_res$coverage_beta_int2[[4]] * 0.01 ) / 1000 ) # Monte Carlo SE of coverage
    
    simsum_res$low_coverage_beta_int2[[4]] <- ( simsum_res$coverage_beta_int2[[4]]*0.01 - 1.96*simsum_res$coverage_mcse_beta_int2[[4]] ) * 100
    simsum_res$high_coverage_beta_int2[[4]] <- ( simsum_res$coverage_beta_int2[[4]]*0.01 + 1.96*simsum_res$coverage_mcse_beta_int2[[4]] ) * 100
    
    simsum_res$power_beta_int2[[4]] <- ( sum( ifelse( abs(data$beta_int2_metareg_B) >= (1.96*data$se_int2_metareg_B), 1 ,0 ) ) / 1000 ) *100 
    simsum_res$power_mcse_beta_int2[[4]] <- sqrt ( (simsum_res$power_beta_int[[4]]*0.01) * ( 1 - simsum_res$power_beta_int[[3]]*0.01 ) / 1000 ) # Monte Carlo SE of power

    assign( paste("simsum_res",s,sep="_"), simsum_res)
    
  }
  
  assign( paste("simsum_res",age_overlap,sep="_"), bind_rows(simsum_res_1,simsum_res_2,simsum_res_3,simsum_res_4,simsum_res_5))
  
}

simsum_A <- bind_rows(simsum_res_0,simsum_res_25,simsum_res_50,simsum_res_75,simsum_res_100)

write_dta (simsum_A,"directory/results.dta", version = 14, label = NULL )

