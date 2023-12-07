
## PART 2: APPLY FIXED-EFFECT AND RANDOM-EFFECTS META-ANALYSE & META - REGRESSION - SCENARIO 1 TO SCENARIO 5 
##
## Pagoni P, Higgins JPT, Lawlor DA, Stergiakouli E, Warrington NM, Morris TT, Tilling K.
## Meta-regression of genome-wide association studies to estimate age-varying genetic effects.
## 
## Panagiota Pagoni
## 31/10/2023

#install.packages("haven")
#install.packages("tidyverse")
#install.packages("Matrix")
#install.packages("meta")
#install.packages("metafor")

library("haven")
library("tidyverse")
library("Matrix")
library("metafor")


# m : sample size of cohorts
# N:Number of studies
# age : age overlap between cohorts

# loop_1 : begining of iteration loop
# loop_2 : end of iteration loop

m <- 1000
n <- 40
age <- list( 0, 25, 50, 75, 100 )

loop_1 <- 1
loop_2 <- 1000

# sc : scenario
# k: number of iterations

for ( age_overlap in age ){

  dir <- paste0("directory/overlap_",age_overlap)
 
  setwd(dir)
  
for( sc in 1:5){
  
  res_final <- as.data.frame(matrix(NA, nrow=1000, ncol= 35))
  
  colnames(res_final) <- c ('iteration','scenario','age_overlap','N', 'centred_age','wmean_all_ages', # Basic info
                            'beta_snp_fixed_metan','se_snp_fixed_metan','i_2_fixed','QE_fixed','QEp_fixed', # FE meta analysis
                            'beta_snp_random_metan','se_snp_random_metan','tau2_random','i_2_random','QE_random','QEp_random', # RE meta analysis
                            
                            'beta_snp_metareg_A','se_beta_snp_metareg_A','beta_int_metareg_A','se_int_metareg_A', # Model_A metareg with x_age
                            'tau2_A','i_2_A', 'QE_A','QEp_A', # Model_A metareg with x_age

                            'beta_snp_metareg_B','se_beta_snp_metareg_B', # Model_B metareg with x_age_centred + x_age2_centred
                            'beta_int_metareg_B','se_int_metareg_B', 
                            'beta_int2_metareg_B','se_int2_metareg_B', 
                            'tau2_B', 'i_2_B', 'QE_B','QEp_B')
  
  for (k in loop_1:loop_2){
    
    data <- read_dta(paste( paste( paste( paste ( paste0(dir,"/sum_stats_sc", sc, sep = ""),
                                          k, sep ="_"),m, sep ="_"),"dta", sep =".")))

    
    data <- data [ 1:n , ]
    
    data$weight <- (1/data$N)
    wmean_all_ages <- weighted.mean(data$x_age,data$weight)
    
    # E(x^2) = Var(x) + [E(x)]^2
    
    data$x_age_centred <- (data$x_age-10) # calculate squared mean age [E(x)]^2 centred in 10 years
    data$var_age <- (data$stdev_age)^2
    data$x_age2_centred <- data$var_age + (data$x_age_centred)^2  # calculate [E(x^2)] to adjust in meta-regression model B
    
    
    res_final$iteration[[k]] <- k
    res_final$scenario[[k]] <- sc
    res_final$age_overlap <- age_overlap
    res_final$N <- m
    res_final$centred_age <- 10
    res_final$wmean_all_ages <- wmean_all_ages

    # Fixed-effect meta analysis
    
    metan_FE <- summary(rma( yi = data$beta_snp, 
                     sei = data$stder_beta_snp,
                     method = "FE"))
    
    
    res_final$beta_snp_fixed_metan[[k]] <- metan_FE$beta
    res_final$se_snp_fixed_metan[[k]] <- metan_FE$se
    res_final$i_2_fixed[[k]] <- metan_FE$I2
    res_final$QE_fixed[[k]] <- metan_FE$QE
    res_final$QEp_fixed[[k]] <- metan_FE$QEp
    
    # Random-effects meta analysis
    
    metan_RE <- summary( rma( yi = data$beta_snp, 
                              sei = data$stder_beta_snp,
                              method = "REML",
                              control=list(stepadj=0.5, maxiter=1000)) )
    
    res_final$beta_snp_random_metan[[k]] <- metan_RE$beta
    res_final$se_snp_random_metan[[k]] <- metan_RE$se
    res_final$tau2_random[[k]] <- metan_RE$tau2
    res_final$i_2_random[[k]] <- metan_RE$I2  
    res_final$QE_random[[k]] <- metan_RE$QE
    res_final$QEp_random[[k]] <- metan_RE$QEp


    # MODEL A: Random-effects meta regression linear term  
    
    mreg_intA <- summary( rma( yi = data$beta_snp, sei = data$stder_beta_snp , mods = ~ x_age_centred, data = data, method = "REML" ) )
    
    res_final$beta_snp_metareg_A[[k]] <- mreg_intA$beta[[1]] # intercept
    res_final$se_beta_snp_metareg_A[[k]] <- mreg_intA$se[[1]]
    
    res_final$beta_int_metareg_A[[k]] <- mreg_intA$beta[[2]] # snp x_age
    res_final$se_int_metareg_A[[k]] <- mreg_intA$se[[2]]
    
    res_final$tau2_A[[k]] <- mreg_intA$tau2
    res_final$i_2_A[[k]] <- mreg_intA$I2
    res_final$QE_A[[k]] <- mreg_intA$QE
    res_final$QEp_A[[k]] <- mreg_intA$QEp
  
    # MODEL B: Random-effects meta regression linear & non - linear term  adjusting for [E(x^2)]
    
    mreg_intB <- summary( rma( yi = data$beta_snp, sei = data$stder_beta_snp , mods= ~ x_age_centred + x_age2_centred, data = data, method = "REML" ) )
    
    res_final$beta_snp_metareg_B[[k]] <- mreg_intB$beta[[1]]
    res_final$se_beta_snp_metareg_B[[k]] <- mreg_intB$se[[1]]
    
    res_final$beta_int_metareg_B[[k]] <- mreg_intB$beta[[2]]
    res_final$se_int_metareg_B[[k]] <- mreg_intB$se[[2]]
    
    res_final$beta_int2_metareg_B[[k]] <- mreg_intB$beta[[3]]
    res_final$se_int2_metareg_B[[k]] <- mreg_intB$se[[3]]
    
    res_final$tau2_B[[k]] <- mreg_intB$tau2
    res_final$i_2_B[[k]] <- mreg_intB$I2  
    res_final$QE_B[[k]] <- mreg_intB$QE
    res_final$QEp_B[[k]] <- mreg_intB$QEp

    print( paste(paste( paste( "Overlap: ",age_overlap, sep = " " ),"Scenario: ",sc, sep = " "), " Iteration: ", k, sep = " ") ) 

  }
  
  write_dta (res_final,paste( paste( paste( paste( paste( paste( dir , "/res_all_sc",sc,sep =""),
                                                          "final",sep = "_"),
                                                   age_overlap,sep ="_"),
                                            m, sep = "_"),
                                     loop_2, sep ="_"),
                              "dta",sep ="."), version = 14, label = NULL )
  
}
}

print ("Script completed")













