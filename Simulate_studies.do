* Part 1: Simulate individual participant data for simulations for: 
*
* Pagoni P, Higgins JPT, Lawlor DA, Stergiakouli E, Warrington NM, Morris TT, Tilling K.
* Meta-regression of genome-wide association studies to estimate age-varying genetic effects.
*
* Panagiota Pagoni
* 31/10/2023


**** Overlap 0%
 
cd "directory"

clear

set seed 1234

* m : indicator for number of participants within each study 1.000

foreach m in 1000 {

* k : indicator for number of iterations

forvalues k = 1/1000 {

clear

set more off

capture set obs 40
generate study = _n

generate beta_0 = 25
generate u_j = rnormal(0,1)

* Create individual level variables (i) *

* number of participants per study *

expand 1000 if `m' == 1000

bysort study: generate id = _n

generate e_ij = rnormal(0,1)

* Create age variable for each study *

generate age_ij = .

replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if  study == 1
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if  study == 2
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if  study == 3
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if  study == 4

replace age_ij = 15 + int(( 19- 15 + 1 )*runiform()) if  study == 5
replace age_ij = 15 + int(( 19- 15 + 1 )*runiform()) if  study == 6
replace age_ij = 15 + int(( 19- 15 + 1 )*runiform()) if  study == 7
replace age_ij = 15 + int(( 19- 15 + 1 )*runiform()) if  study == 8

replace age_ij = 20 + int(( 24 - 20 + 1 )*runiform()) if  study == 9
replace age_ij = 20 + int(( 24 - 20 + 1 )*runiform()) if  study == 10
replace age_ij = 20 + int(( 24 - 20 + 1 )*runiform()) if  study == 11
replace age_ij = 20 + int(( 24 - 20 + 1 )*runiform()) if  study == 12

replace age_ij = 25 + int(( 29 - 25 + 1 )*runiform()) if  study == 13
replace age_ij = 25 + int(( 29 - 25 + 1 )*runiform()) if  study == 14
replace age_ij = 25 + int(( 29 - 25 + 1 )*runiform()) if  study == 15
replace age_ij = 25 + int(( 29 - 25 + 1 )*runiform()) if  study == 16

replace age_ij = 30 + int(( 34 - 30 + 1 )*runiform()) if  study == 17
replace age_ij = 30 + int(( 34 - 30 + 1 )*runiform()) if  study == 18
replace age_ij = 30 + int(( 34 - 30 + 1 )*runiform()) if  study == 19
replace age_ij = 30 + int(( 34 - 30 + 1 )*runiform()) if  study == 20

replace age_ij = 35 + int(( 39 - 35 + 1 )*runiform()) if  study == 21
replace age_ij = 35 + int(( 39 - 35 + 1 )*runiform()) if  study == 22
replace age_ij = 35 + int(( 39 - 35 + 1 )*runiform()) if  study == 23
replace age_ij = 35 + int(( 39 - 35 + 1 )*runiform()) if  study == 24

replace age_ij = 40 + int(( 44 - 40 + 1 )*runiform()) if  study == 25
replace age_ij = 40 + int(( 44 - 40 + 1 )*runiform()) if  study == 26
replace age_ij = 40 + int(( 44 - 40 + 1 )*runiform()) if  study == 27
replace age_ij = 40 + int(( 44 - 40 + 1 )*runiform()) if  study == 28

replace age_ij = 45 + int(( 49 - 45 + 1 )*runiform()) if  study == 29
replace age_ij = 45 + int(( 49 - 45 + 1 )*runiform()) if  study == 30
replace age_ij = 45 + int(( 49 - 45 + 1 )*runiform()) if  study == 31
replace age_ij = 45 + int(( 49 - 45 + 1 )*runiform()) if  study == 32

replace age_ij = 50 + int(( 54 - 50 + 1 )*runiform()) if  study == 33
replace age_ij = 50 + int(( 54 - 50 + 1 )*runiform()) if  study == 34
replace age_ij = 50 + int(( 54 - 50 + 1 )*runiform()) if  study == 35
replace age_ij = 50 + int(( 54 - 50 + 1 )*runiform()) if  study == 36

replace age_ij = 55 + int(( 59 - 55 + 1 )*runiform()) if  study == 37
replace age_ij = 55 + int(( 59 - 55 + 1 )*runiform()) if  study == 38
replace age_ij = 55 + int(( 59 - 55 + 1 )*runiform()) if  study == 39
replace age_ij = 55 + int(( 59 - 55 + 1 )*runiform()) if  study == 40

* Create SNP_ij = {0,1,2} *

generate SNP_ij = rbinomial(2,0.2)
tab SNP_ij

* Check if variables are correlated *

correlate age_ij SNP_ij e_ij u_j

* Set values for simulated betas *

gen beta_SNP = 1.5
gen beta_age = 0.01
gen beta_SNP_age = 0.02
gen beta_age2 = 0.002
gen beta_SNP_age_2 = 0.001

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 1: BMI is dependent on both age and genotoype
***     		  BMI_ij = beta_0 + age_ij + ( SNP_ij + u_j ) + e_ij  

generate bmi_ij_sc1 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc1
     
* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 2: BMI is dependent on both age and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc2 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc2

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 3: BMI is dependent on both age^2 and genotype 
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + e_ij + u_j 

generate bmi_ij_sc3 = beta_0 + beta_age*age_ij + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc3

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 4: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc4 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc4

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 5: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype and an interaction between age^2 and genotype
***     		  BMI_ij = beta_0 + age_ij + age_ij^2 + SNP_ij + age_ij*SNP_ij + age_ij^2*SNP_ij e_ij + u_j 

generate bmi_ij_sc5 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + beta_SNP_age_2*age_ij*age_ij*SNP_ij + e_ij
sum bmi_ij_sc5
	   
* s : indicator for scenarios

forvalues s = 1/5 {

* Save betas and SE by study
 
matrix define sum_stats_sc`s' = J( 40,12,.)

 * i : indicator for study
  
forvalues i = 1/40 {


    matrix sum_stats_sc`s'[`i',1] = `k'
	matrix sum_stats_sc`s'[`i',2] = `s'

	* Exctract N, mean age, sd age from each study*
	
	sum age_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',3] = `i'
	matrix sum_stats_sc`s'[`i',4] = r(N)
	matrix sum_stats_sc`s'[`i',5] = r(mean)
	matrix sum_stats_sc`s'[`i',6] = r(sd)

	* Fit a linear model in each study 
	
	regress bmi_ij_sc`s' age_ij SNP_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',7] = _b[_cons]
	matrix sum_stats_sc`s'[`i',8] = _se[_cons]

	matrix sum_stats_sc`s'[`i',9] = _b[age_ij]
	matrix sum_stats_sc`s'[`i',10] = _se[age_ij]
	
	matrix sum_stats_sc`s'[`i',11] = _b[SNP_ij]
	matrix sum_stats_sc`s'[`i',12] = _se[SNP_ij]
		
	}
	
	matrix list sum_stats_sc`s'

preserve

svmat sum_stats_sc`s'

keep sum_stats_sc`s'*

rename sum_stats_sc`s'1 iteration
rename sum_stats_sc`s'2 scenario
rename sum_stats_sc`s'3 study
rename sum_stats_sc`s'4 N
rename sum_stats_sc`s'5 x_age
rename sum_stats_sc`s'6 stdev_age
rename sum_stats_sc`s'7 beta_0
rename sum_stats_sc`s'8 stder_beta_0
rename sum_stats_sc`s'9 beta_age
rename sum_stats_sc`s'10 stder_beta_age
rename sum_stats_sc`s'11 beta_snp
rename sum_stats_sc`s'12 stder_beta_snp

save sum_stats_sc`s'_`k'_`m',replace

restore

}
}
}


**** Overlap 25%

* m : indicator for number of participants within each study 1.000


foreach m in 1000 {

* k : indicator for number of iterations

forvalues k = 1/1000 {

clear

set more off

capture set obs 40
generate study = _n

generate beta_0 = 25
generate u_j = rnormal(0,1)

* Create individual level variables (i) *

* number of participants per study *

expand 1000 if `m' == 1000

bysort study: generate id = _n

generate e_ij = rnormal(0,1)

* Create age variable for each study *

generate age_ij = .
 
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 1
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 2
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 3
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 4

sum age_ij if ( study > = 1 & study < = 4) ,d


replace age_ij = r(p75) + int(( 19- r(p75) + 1 )*runiform()) if study == 5
replace age_ij = r(p75) + int(( 19- r(p75) + 1 )*runiform()) if study == 6
replace age_ij = r(p75) + int(( 19- r(p75) + 1 )*runiform()) if study == 7
replace age_ij = r(p75) + int(( 19- r(p75) + 1 )*runiform()) if study == 8

sum age_ij if  ( study > 4 & study < = 8 ) ,d

replace age_ij = r(p75) + int(( 24 - r(p75) + 1 )*runiform()) if study == 9
replace age_ij = r(p75) + int(( 24 - r(p75) + 1 )*runiform()) if study == 10
replace age_ij = r(p75) + int(( 24 - r(p75) + 1 )*runiform()) if study == 11
replace age_ij = r(p75) + int(( 24 - r(p75) + 1 )*runiform()) if study == 12

sum age_ij if ( study > 8 & study < = 12 ) ,d

replace age_ij = r(p75) + int(( 29 - r(p75) + 1 )*runiform()) if study == 13
replace age_ij = r(p75) + int(( 29 - r(p75) + 1 )*runiform()) if study == 14
replace age_ij = r(p75) + int(( 29 - r(p75) + 1 )*runiform()) if study == 15
replace age_ij = r(p75) + int(( 29 - r(p75) + 1 )*runiform()) if study == 16

sum age_ij if ( study > 12 & study < = 16 ) ,d

replace age_ij = r(p75) + int(( 34 - r(p75) + 1 )*runiform()) if study == 17
replace age_ij = r(p75) + int(( 34 - r(p75) + 1 )*runiform()) if study == 18
replace age_ij = r(p75) + int(( 34 - r(p75) + 1 )*runiform()) if study == 19
replace age_ij = r(p75) + int(( 34 - r(p75) + 1 )*runiform()) if study == 20

sum age_ij if ( study > 16 & study < = 20 ) ,d

replace age_ij = r(p75) + int(( 39 - r(p75) + 1 )*runiform()) if study == 21
replace age_ij = r(p75) + int(( 39 - r(p75) + 1 )*runiform()) if study == 22
replace age_ij = r(p75) + int(( 39 - r(p75) + 1 )*runiform()) if study == 23
replace age_ij = r(p75) + int(( 39 - r(p75) + 1 )*runiform()) if study == 24

sum age_ij if ( study > 20 & study < = 24 ) ,d

replace age_ij = r(p75) + int(( 44 - r(p75) + 1 )*runiform()) if study == 25
replace age_ij = r(p75) + int(( 44 - r(p75) + 1 )*runiform()) if study == 26
replace age_ij = r(p75) + int(( 44 - r(p75) + 1 )*runiform()) if study == 27
replace age_ij = r(p75) + int(( 44 - r(p75) + 1 )*runiform()) if study == 28

sum age_ij if ( study > 24 & study < = 28 ) ,d

replace age_ij = r(p75) + int(( 49 - r(p75) + 1 )*runiform()) if study == 29
replace age_ij = r(p75) + int(( 49 - r(p75) + 1 )*runiform()) if study == 30
replace age_ij = r(p75) + int(( 49 - r(p75) + 1 )*runiform()) if study == 31
replace age_ij = r(p75) + int(( 49 - r(p75) + 1 )*runiform()) if study == 32

sum age_ij if ( study > 28 & study < = 32 ) ,d

replace age_ij = r(p75) + int(( 54 - r(p75) + 1 )*runiform()) if study == 33
replace age_ij = r(p75) + int(( 54 - r(p75) + 1 )*runiform()) if study == 34
replace age_ij = r(p75) + int(( 54 - r(p75) + 1 )*runiform()) if study == 35
replace age_ij = r(p75) + int(( 54 - r(p75) + 1 )*runiform()) if study == 36

sum age_ij if ( study > 32 & study < = 36 ) ,d

replace age_ij = r(p75) + int(( 59 - r(p75) + 1 )*runiform()) if study == 37
replace age_ij = r(p75) + int(( 59 - r(p75) + 1 )*runiform()) if study == 38
replace age_ij = r(p75) + int(( 59 - r(p75) + 1 )*runiform()) if study == 39
replace age_ij = r(p75) + int(( 59 - r(p75) + 1 )*runiform()) if study == 40

sum age_ij if ( study > 36 & study < = 40 ) ,d


*scatter age_ij study

* Create SNP_ij = {0,1,2} *

generate SNP_ij = rbinomial(2,0.2)
tab SNP_ij

* Check if variables are correlated *

correlate age_ij SNP_ij e_ij u_j

* Set values for simulated betas *

gen beta_SNP = 1.5
gen beta_age = 0.01
gen beta_SNP_age = 0.02
gen beta_age2 = 0.002
gen beta_SNP_age_2 = 0.001

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 1: BMI is dependent on both age and genotoype
***     		  BMI_ij = beta_0 + age_ij + ( SNP_ij + u_j ) + e_ij  

generate bmi_ij_sc1 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc1

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 2: BMI is dependent on both age and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc2 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc2

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 3: BMI is dependent on both age^2 and genotype 
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + e_ij + u_j 

generate bmi_ij_sc3 = beta_0 + beta_age*age_ij + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc3

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 4: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc4 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc4

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 5: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype and an interaction between age^2 and genotype
***     		  BMI_ij = beta_0 + age_ij + age_ij^2 + SNP_ij + age_ij*SNP_ij + age_ij^2*SNP_ij e_ij + u_j 

generate bmi_ij_sc5 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + beta_SNP_age_2*age_ij*age_ij*SNP_ij + e_ij
sum bmi_ij_sc5
	   
* s : indicator for scenarios

forvalues s = 1/5 {

* Save betas and SE by study
 
matrix define sum_stats_sc`s' = J( 40,12,.)

 * i : indicator for study
  
forvalues i = 1/40 {


    matrix sum_stats_sc`s'[`i',1] = `k'
	matrix sum_stats_sc`s'[`i',2] = `s'

	* Exctract N, mean age, sd age from each study*
	
	sum age_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',3] = `i'
	matrix sum_stats_sc`s'[`i',4] = r(N)
	matrix sum_stats_sc`s'[`i',5] = r(mean)
	matrix sum_stats_sc`s'[`i',6] = r(sd)

	* Fit a linear model in each study 
	
	regress bmi_ij_sc`s' age_ij SNP_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',7] = _b[_cons]
	matrix sum_stats_sc`s'[`i',8] = _se[_cons]

	matrix sum_stats_sc`s'[`i',9] = _b[age_ij]
	matrix sum_stats_sc`s'[`i',10] = _se[age_ij]
	
	matrix sum_stats_sc`s'[`i',11] = _b[SNP_ij]
	matrix sum_stats_sc`s'[`i',12] = _se[SNP_ij]
		
	}
	
	matrix list sum_stats_sc`s'

preserve

svmat sum_stats_sc`s'

keep sum_stats_sc`s'*

rename sum_stats_sc`s'1 iteration
rename sum_stats_sc`s'2 scenario
rename sum_stats_sc`s'3 study
rename sum_stats_sc`s'4 N
rename sum_stats_sc`s'5 x_age
rename sum_stats_sc`s'6 stdev_age
rename sum_stats_sc`s'7 beta_0
rename sum_stats_sc`s'8 stder_beta_0
rename sum_stats_sc`s'9 beta_age
rename sum_stats_sc`s'10 stder_beta_age
rename sum_stats_sc`s'11 beta_snp
rename sum_stats_sc`s'12 stder_beta_snp

save sum_stats_sc`s'_`k'_`m',replace

restore

}
}
}


**** Overlap 50%

* m : indicator for number of participants within each study 1.000

foreach m in 1000 {

* k : indicator for number of iterations

forvalues k = 1/1000 {

clear

set more off

capture set obs 40
generate study = _n

generate beta_0 = 25
generate u_j = rnormal(0,1)

* Create individual level variables (i) *

* number of participants per study *

expand 1000 if `m' == 1000

bysort study: generate id = _n

generate e_ij = rnormal(0,1)

* Create age variable for each study *

generate age_ij = .
 
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 1
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 2
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 3
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 4

sum age_ij if ( study > = 1 & study < = 4), d

replace age_ij = r(p50) + int(( 19- r(p50) + 1 )*runiform()) if study == 5
replace age_ij = r(p50) + int(( 19- r(p50) + 1 )*runiform()) if study == 6
replace age_ij = r(p50) + int(( 19- r(p50) + 1 )*runiform()) if study == 7
replace age_ij = r(p50) + int(( 19- r(p50) + 1 )*runiform()) if study == 8


sum age_ij if  ( study > 4 & study < = 8 ), d

replace age_ij = r(p50) + int(( 24 - r(p50) + 1 )*runiform()) if study == 9
replace age_ij = r(p50) + int(( 24 - r(p50) + 1 )*runiform()) if study == 10
replace age_ij = r(p50) + int(( 24 - r(p50) + 1 )*runiform()) if study == 11
replace age_ij = r(p50) + int(( 24 - r(p50) + 1 )*runiform()) if study == 12

sum age_ij if ( study > 8 & study < = 12 ), d

replace age_ij = r(p50) + int(( 29 - r(p50) + 1 )*runiform()) if study == 13
replace age_ij = r(p50) + int(( 29 - r(p50) + 1 )*runiform()) if study == 14
replace age_ij = r(p50) + int(( 29 - r(p50) + 1 )*runiform()) if study == 15
replace age_ij = r(p50) + int(( 29 - r(p50) + 1 )*runiform()) if study == 16

sum age_ij if ( study > 12 & study < = 16 ), d

replace age_ij = r(p50) + int(( 34 - r(p50) + 1 )*runiform()) if study == 17
replace age_ij = r(p50) + int(( 34 - r(p50) + 1 )*runiform()) if study == 18
replace age_ij = r(p50) + int(( 34 - r(p50) + 1 )*runiform()) if study == 19
replace age_ij = r(p50) + int(( 34 - r(p50) + 1 )*runiform()) if study == 20

sum age_ij if ( study > 16 & study < = 20 ), d

replace age_ij = r(p50) + int(( 39 - r(p50) + 1 )*runiform()) if study == 21
replace age_ij = r(p50) + int(( 39 - r(p50) + 1 )*runiform()) if study == 22
replace age_ij = r(p50) + int(( 39 - r(p50) + 1 )*runiform()) if study == 23
replace age_ij = r(p50) + int(( 39 - r(p50) + 1 )*runiform()) if study == 24

sum age_ij if ( study > 20 & study < = 24 ), d

replace age_ij = r(p50) + int(( 44 - r(p50) + 1 )*runiform()) if study == 25
replace age_ij = r(p50) + int(( 44 - r(p50) + 1 )*runiform()) if study == 26
replace age_ij = r(p50) + int(( 44 - r(p50) + 1 )*runiform()) if study == 27
replace age_ij = r(p50) + int(( 44 - r(p50) + 1 )*runiform()) if study == 28

sum age_ij if ( study > 24 & study < = 28 ), d

replace age_ij = r(p50) + int(( 49 - r(p50) + 1 )*runiform()) if study == 29
replace age_ij = r(p50) + int(( 49 - r(p50) + 1 )*runiform()) if study == 30
replace age_ij = r(p50) + int(( 49 - r(p50) + 1 )*runiform()) if study == 31
replace age_ij = r(p50) + int(( 49 - r(p50) + 1 )*runiform()) if study == 32

sum age_ij if ( study > 28 & study < = 32 ), d

replace age_ij = r(p50) + int(( 54 - r(p50) + 1 )*runiform()) if study == 33
replace age_ij = r(p50) + int(( 54 - r(p50) + 1 )*runiform()) if study == 34
replace age_ij = r(p50) + int(( 54 - r(p50) + 1 )*runiform()) if study == 35
replace age_ij = r(p50) + int(( 54 - r(p50) + 1 )*runiform()) if study == 36

sum age_ij if ( study > 32 & study < = 36 ), d

replace age_ij = r(p50) + int(( 59 - r(p50) + 1 )*runiform()) if study == 37
replace age_ij = r(p50) + int(( 59 - r(p50) + 1 )*runiform()) if study == 38
replace age_ij = r(p50) + int(( 59 - r(p50) + 1 )*runiform()) if study == 39
replace age_ij = r(p50) + int(( 59 - r(p50) + 1 )*runiform()) if study == 40

sum age_ij if ( study > 36 & study < = 40 ), d


*scatter age_ij study 

* Create SNP_ij = {0,1,2} *

generate SNP_ij = rbinomial(2,0.2)
tab SNP_ij

* Check if variables are correlated *

correlate age_ij SNP_ij e_ij u_j

* Set values for simulated betas *

gen beta_SNP = 1.5
gen beta_age = 0.01
gen beta_SNP_age = 0.02
gen beta_age2 = 0.002
gen beta_SNP_age_2 = 0.001

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 1: BMI is dependent on both age and genotoype
***     		  BMI_ij = beta_0 + age_ij + ( SNP_ij + u_j ) + e_ij  

generate bmi_ij_sc1 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc1
   
* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 2: BMI is dependent on both age and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc2 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc2

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 3: BMI is dependent on both age^2 and genotype 
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + e_ij + u_j 

generate bmi_ij_sc3 = beta_0 + beta_age*age_ij + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc3

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 4: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc4 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc4

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 5: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype and an interaction between age^2 and genotype
***     		  BMI_ij = beta_0 + age_ij + age_ij^2 + SNP_ij + age_ij*SNP_ij + age_ij^2*SNP_ij e_ij + u_j 

generate bmi_ij_sc5 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + beta_SNP_age_2*age_ij*age_ij*SNP_ij + e_ij
sum bmi_ij_sc5

*scatter age_ij study
	   
* s : indicator for scenarios

forvalues s = 1/5 {

* Save betas and SE by study
 
matrix define sum_stats_sc`s' = J( 40,12,.)

 * i : indicator for study
  
forvalues i = 1/40 {


    matrix sum_stats_sc`s'[`i',1] = `k'
	matrix sum_stats_sc`s'[`i',2] = `s'

	* Exctract N, mean age, sd age from each study*
	
	sum age_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',3] = `i'
	matrix sum_stats_sc`s'[`i',4] = r(N)
	matrix sum_stats_sc`s'[`i',5] = r(mean)
	matrix sum_stats_sc`s'[`i',6] = r(sd)

	* Fit a linear model in each study 
	
	regress bmi_ij_sc`s' age_ij SNP_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',7] = _b[_cons]
	matrix sum_stats_sc`s'[`i',8] = _se[_cons]

	matrix sum_stats_sc`s'[`i',9] = _b[age_ij]
	matrix sum_stats_sc`s'[`i',10] = _se[age_ij]
	
	matrix sum_stats_sc`s'[`i',11] = _b[SNP_ij]
	matrix sum_stats_sc`s'[`i',12] = _se[SNP_ij]
		
	}
	
	matrix list sum_stats_sc`s'

preserve

svmat sum_stats_sc`s'

keep sum_stats_sc`s'*

rename sum_stats_sc`s'1 iteration
rename sum_stats_sc`s'2 scenario
rename sum_stats_sc`s'3 study
rename sum_stats_sc`s'4 N
rename sum_stats_sc`s'5 x_age
rename sum_stats_sc`s'6 stdev_age
rename sum_stats_sc`s'7 beta_0
rename sum_stats_sc`s'8 stder_beta_0
rename sum_stats_sc`s'9 beta_age
rename sum_stats_sc`s'10 stder_beta_age
rename sum_stats_sc`s'11 beta_snp
rename sum_stats_sc`s'12 stder_beta_snp

save sum_stats_sc`s'_`k'_`m',replace

restore

}
}
}



**** Overlap 75%

* m : indicator for number of participants within each study 1.000


foreach m in 1000 {

* k : indicator for number of iterations

forvalues k = 1/1000 {

clear

set more off

capture set obs 40
generate study = _n

generate beta_0 = 25
generate u_j = rnormal(0,1)

* Create individual level variables (i) *

* number of participants per study *

expand 1000 if `m' == 1000

bysort study: generate id = _n

generate e_ij = rnormal(0,1)

* Create age variable for each study *

generate age_ij = .
 
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 1
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 2
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 3
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 4

sum age_ij if ( study > = 1 & study < = 4), d 

replace age_ij = r(p25) + int(( 19- r(p25) + 1 )*runiform()) if study == 5
replace age_ij = r(p25) + int(( 19- r(p25) + 1 )*runiform()) if study == 6
replace age_ij = r(p25) + int(( 19- r(p25) + 1 )*runiform()) if study == 7
replace age_ij = r(p25) + int(( 19- r(p25) + 1 )*runiform()) if study == 8

sum age_ij if  ( study > 4 & study < = 8 ), d 

replace age_ij = r(p25) + int(( 24 - r(p25) + 1 )*runiform()) if study == 9
replace age_ij = r(p25) + int(( 24 - r(p25) + 1 )*runiform()) if study == 10
replace age_ij = r(p25) + int(( 24 - r(p25) + 1 )*runiform()) if study == 11
replace age_ij = r(p25) + int(( 24 - r(p25) + 1 )*runiform()) if study == 12

sum age_ij if ( study > 8 & study < = 12 ), d 

replace age_ij = r(p25) + int(( 29 - r(p25) + 1 )*runiform()) if study == 13
replace age_ij = r(p25) + int(( 29 - r(p25) + 1 )*runiform()) if study == 14
replace age_ij = r(p25) + int(( 29 - r(p25) + 1 )*runiform()) if study == 15
replace age_ij = r(p25) + int(( 29 - r(p25) + 1 )*runiform()) if study == 16

sum age_ij if ( study > 12 & study < = 16 ), d 

replace age_ij = r(p25) + int(( 34 - r(p25) + 1 )*runiform()) if study == 17
replace age_ij = r(p25) + int(( 34 - r(p25) + 1 )*runiform()) if study == 18
replace age_ij = r(p25) + int(( 34 - r(p25) + 1 )*runiform()) if study == 19
replace age_ij = r(p25) + int(( 34 - r(p25) + 1 )*runiform()) if study == 20

sum age_ij if ( study > 16 & study < = 20 ), d 

replace age_ij = r(p25) + int(( 39 - r(p25) + 1 )*runiform()) if study == 21
replace age_ij = r(p25) + int(( 39 - r(p25) + 1 )*runiform()) if study == 22
replace age_ij = r(p25) + int(( 39 - r(p25) + 1 )*runiform()) if study == 23
replace age_ij = r(p25) + int(( 39 - r(p25) + 1 )*runiform()) if study == 24

sum age_ij if ( study > 20 & study < = 24 ), d 

replace age_ij = r(p25) + int(( 44 - r(p25) + 1 )*runiform()) if study == 25
replace age_ij = r(p25) + int(( 44 - r(p25) + 1 )*runiform()) if study == 26
replace age_ij = r(p25) + int(( 44 - r(p25) + 1 )*runiform()) if study == 27
replace age_ij = r(p25) + int(( 44 - r(p25) + 1 )*runiform()) if study == 28

sum age_ij if ( study > 24 & study < = 28 ), d 

replace age_ij = r(p25) + int(( 49 - r(p25) + 1 )*runiform()) if study == 29
replace age_ij = r(p25) + int(( 49 - r(p25) + 1 )*runiform()) if study == 30
replace age_ij = r(p25) + int(( 49 - r(p25) + 1 )*runiform()) if study == 31
replace age_ij = r(p25) + int(( 49 - r(p25) + 1 )*runiform()) if study == 32

sum age_ij if ( study > 28 & study < = 32 ), d 

replace age_ij = r(p25) + int(( 54 - r(p25) + 1 )*runiform()) if study == 33
replace age_ij = r(p25) + int(( 54 - r(p25) + 1 )*runiform()) if study == 34
replace age_ij = r(p25) + int(( 54 - r(p25) + 1 )*runiform()) if study == 35
replace age_ij = r(p25) + int(( 54 - r(p25) + 1 )*runiform()) if study == 36

sum age_ij if ( study > 32 & study < = 36 ), d 

replace age_ij = r(p25) + int(( 59 - r(p25) + 1 )*runiform()) if study == 37
replace age_ij = r(p25) + int(( 59 - r(p25) + 1 )*runiform()) if study == 38
replace age_ij = r(p25) + int(( 59 - r(p25) + 1 )*runiform()) if study == 39
replace age_ij = r(p25) + int(( 59 - r(p25) + 1 )*runiform()) if study == 40

sum age_ij if ( study > 36 & study < = 40 ), d 

* Create SNP_ij = {0,1,2} *

generate SNP_ij = rbinomial(2,0.2)
tab SNP_ij

* Check if variables are correlated *

correlate age_ij SNP_ij e_ij u_j

* Set values for simulated betas *

gen beta_SNP = 1.5
gen beta_age = 0.01
gen beta_SNP_age = 0.02
gen beta_age2 = 0.002
gen beta_SNP_age_2 = 0.001


* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 1: BMI is dependent on both age and genotoype
***     		  BMI_ij = beta_0 + age_ij + ( SNP_ij + u_j ) + e_ij  

generate bmi_ij_sc1 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc1
   
* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 2: BMI is dependent on both age and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc2 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc2

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 3: BMI is dependent on both age^2 and genotype 
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + e_ij + u_j 

generate bmi_ij_sc3 = beta_0 + beta_age*age_ij + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc3

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 4: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc4 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc4

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 5: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype and an interaction between age^2 and genotype
***     		  BMI_ij = beta_0 + age_ij + age_ij^2 + SNP_ij + age_ij*SNP_ij + age_ij^2*SNP_ij e_ij + u_j 

generate bmi_ij_sc5 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + beta_SNP_age_2*age_ij*age_ij*SNP_ij + e_ij
sum bmi_ij_sc5
	   
* s : indicator for scenarios

forvalues s = 1/5 {

* Save betas and SE by study
 
matrix define sum_stats_sc`s' = J( 40,12,.)

 * i : indicator for study
  
forvalues i = 1/40 {


    matrix sum_stats_sc`s'[`i',1] = `k'
	matrix sum_stats_sc`s'[`i',2] = `s'

	* Exctract N, mean age, sd age from each study*
	
	sum age_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',3] = `i'
	matrix sum_stats_sc`s'[`i',4] = r(N)
	matrix sum_stats_sc`s'[`i',5] = r(mean)
	matrix sum_stats_sc`s'[`i',6] = r(sd)

	* Fit a linear model in each study 
	
	regress bmi_ij_sc`s' age_ij SNP_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',7] = _b[_cons]
	matrix sum_stats_sc`s'[`i',8] = _se[_cons]

	matrix sum_stats_sc`s'[`i',9] = _b[age_ij]
	matrix sum_stats_sc`s'[`i',10] = _se[age_ij]
	
	matrix sum_stats_sc`s'[`i',11] = _b[SNP_ij]
	matrix sum_stats_sc`s'[`i',12] = _se[SNP_ij]
		
	}
	
	matrix list sum_stats_sc`s'

preserve

svmat sum_stats_sc`s'

keep sum_stats_sc`s'*

rename sum_stats_sc`s'1 iteration
rename sum_stats_sc`s'2 scenario
rename sum_stats_sc`s'3 study
rename sum_stats_sc`s'4 N
rename sum_stats_sc`s'5 x_age
rename sum_stats_sc`s'6 stdev_age
rename sum_stats_sc`s'7 beta_0
rename sum_stats_sc`s'8 stder_beta_0
rename sum_stats_sc`s'9 beta_age
rename sum_stats_sc`s'10 stder_beta_age
rename sum_stats_sc`s'11 beta_snp
rename sum_stats_sc`s'12 stder_beta_snp

save sum_stats_sc`s'_`k'_`m',replace

restore

}
}
}


**** Overlap 100%

* m : indicator for number of participants within each study 1.000


foreach m in 1000 {

* k : indicator for number of iterations

forvalues k = 1/1000 {

clear

set more off

capture set obs 40
generate study = _n

generate beta_0 = 25
generate u_j = rnormal(0,1)

* Create individual level variables (i) *

* number of participants per study *

expand 1000 if `m' == 1000

bysort study: generate id = _n

generate e_ij = rnormal(0,1)

* Create age variable for each study *

generate age_ij = .
 
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 1
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 2
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 3
replace age_ij = 10 + int(( 14 - 10 + 1 )*runiform()) if study == 4
sum age_ij if study == 1, d 

replace age_ij = 10 + int(( 19- 10 + 1 )*runiform()) if  study == 5
replace age_ij = 10 + int(( 19- 10 + 1 )*runiform()) if  study == 6
replace age_ij = 10 + int(( 19- 10 + 1 )*runiform()) if  study == 7
replace age_ij = 10 + int(( 19- 10 + 1 )*runiform()) if  study == 8

replace age_ij = 10 + int(( 24 - 10 + 1 )*runiform()) if  study == 9
replace age_ij = 10 + int(( 24 - 10 + 1 )*runiform()) if  study == 10
replace age_ij = 10 + int(( 24 - 10 + 1 )*runiform()) if  study == 11
replace age_ij = 10 + int(( 24 - 10 + 1 )*runiform()) if  study == 12

replace age_ij = 10 + int(( 29 - 10 + 1 )*runiform()) if  study == 13
replace age_ij = 10 + int(( 29 - 10 + 1 )*runiform()) if  study == 14
replace age_ij = 10 + int(( 29 - 10 + 1 )*runiform()) if  study == 15
replace age_ij = 10 + int(( 29 - 10 + 1 )*runiform()) if  study == 16

replace age_ij = 10 + int(( 34 - 10 + 1 )*runiform()) if  study == 17
replace age_ij = 10 + int(( 34 - 10 + 1 )*runiform()) if  study == 18
replace age_ij = 10 + int(( 34 - 10 + 1 )*runiform()) if  study == 19
replace age_ij = 10 + int(( 34 - 10 + 1 )*runiform()) if  study == 20

replace age_ij = 10 + int(( 39 - 10 + 1 )*runiform()) if  study == 21
replace age_ij = 10 + int(( 39 - 10 + 1 )*runiform()) if  study == 22
replace age_ij = 10 + int(( 39 - 10 + 1 )*runiform()) if  study == 23
replace age_ij = 10 + int(( 39 - 10 + 1 )*runiform()) if  study == 24

replace age_ij = 10 + int(( 44 - 10 + 1 )*runiform()) if  study == 25
replace age_ij = 10 + int(( 44 - 10 + 1 )*runiform()) if  study == 26
replace age_ij = 10 + int(( 44 - 10 + 1 )*runiform()) if  study == 27
replace age_ij = 10 + int(( 44 - 10 + 1 )*runiform()) if  study == 28

replace age_ij = 10 + int(( 49 - 10 + 1 )*runiform()) if  study == 29
replace age_ij = 10 + int(( 49 - 10 + 1 )*runiform()) if  study == 30
replace age_ij = 10 + int(( 49 - 10 + 1 )*runiform()) if  study == 31
replace age_ij = 10 + int(( 49 - 10 + 1 )*runiform()) if  study == 32

replace age_ij = 10 + int(( 54 - 10 + 1 )*runiform()) if  study == 33
replace age_ij = 10 + int(( 54 - 10 + 1 )*runiform()) if  study == 34
replace age_ij = 10 + int(( 54 - 10 + 1 )*runiform()) if  study == 35
replace age_ij = 10 + int(( 54 - 10 + 1 )*runiform()) if  study == 36

replace age_ij = 10 + int(( 59 - 10 + 1 )*runiform()) if  study == 37
replace age_ij = 10 + int(( 59 - 10 + 1 )*runiform()) if  study == 38
replace age_ij = 10 + int(( 59 - 10 + 1 )*runiform()) if  study == 39
replace age_ij = 10 + int(( 59 - 10 + 1 )*runiform()) if  study == 40

*scatter age_ij study

* Create SNP_ij = {0,1,2} *

generate SNP_ij = rbinomial(2,0.2)
tab SNP_ij

* Check if variables are correlated *

correlate age_ij SNP_ij e_ij u_j

* Set values for simulated betas *

gen beta_SNP = 1.5
gen beta_age = 0.01
gen beta_SNP_age = 0.02
gen beta_age2 = 0.002
gen beta_SNP_age_2 = 0.001

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 1: BMI is dependent on both age and genotoype
***     		  BMI_ij = beta_0 + age_ij + ( SNP_ij + u_j ) + e_ij  

generate bmi_ij_sc1 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc1
   
* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 2: BMI is dependent on both age and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc2 = beta_0 + beta_age*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc2

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 3: BMI is dependent on both age^2 and genotype 
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + e_ij + u_j 

generate bmi_ij_sc3 = beta_0 + beta_age*age_ij + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + e_ij
sum bmi_ij_sc3

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 4: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype
***     		  BMI_ij = beta_0 + age_ij^2 + SNP_ij + age_ij*SNP_ij + e_ij + u_j 

generate bmi_ij_sc4 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + e_ij
sum bmi_ij_sc4

* ------------------------------------------------------------------------------------------------------------------------------------------------ *

*** SCENARIO 5: BMI is dependent on both age^2 and genotype and there is an interaction between age and genotype and an interaction between age^2 and genotype
***     		  BMI_ij = beta_0 + age_ij + age_ij^2 + SNP_ij + age_ij*SNP_ij + age_ij^2*SNP_ij e_ij + u_j 

generate bmi_ij_sc5 = beta_0 + beta_age*age + beta_age2*age_ij*age_ij + (beta_SNP + u_j)*SNP_ij + beta_SNP_age*age_ij*SNP_ij + beta_SNP_age_2*age_ij*age_ij*SNP_ij + e_ij
sum bmi_ij_sc5
	   
* s : indicator for scenarios

forvalues s = 1/5 {

* Save betas and SE by study
 
matrix define sum_stats_sc`s' = J( 40,12,.)

 * i : indicator for study
  
forvalues i = 1/40 {


    matrix sum_stats_sc`s'[`i',1] = `k'
	matrix sum_stats_sc`s'[`i',2] = `s'

	* Exctract N, mean age, sd age from each study*
	
	sum age_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',3] = `i'
	matrix sum_stats_sc`s'[`i',4] = r(N)
	matrix sum_stats_sc`s'[`i',5] = r(mean)
	matrix sum_stats_sc`s'[`i',6] = r(sd)

	* Fit a linear model in each study 
	
	regress bmi_ij_sc`s' age_ij SNP_ij if study == `i'
	
	matrix sum_stats_sc`s'[`i',7] = _b[_cons]
	matrix sum_stats_sc`s'[`i',8] = _se[_cons]

	matrix sum_stats_sc`s'[`i',9] = _b[age_ij]
	matrix sum_stats_sc`s'[`i',10] = _se[age_ij]
	
	matrix sum_stats_sc`s'[`i',11] = _b[SNP_ij]
	matrix sum_stats_sc`s'[`i',12] = _se[SNP_ij]
		
	}
	
	matrix list sum_stats_sc`s'

preserve

svmat sum_stats_sc`s'

keep sum_stats_sc`s'*

rename sum_stats_sc`s'1 iteration
rename sum_stats_sc`s'2 scenario
rename sum_stats_sc`s'3 study
rename sum_stats_sc`s'4 N
rename sum_stats_sc`s'5 x_age
rename sum_stats_sc`s'6 stdev_age
rename sum_stats_sc`s'7 beta_0
rename sum_stats_sc`s'8 stder_beta_0
rename sum_stats_sc`s'9 beta_age
rename sum_stats_sc`s'10 stder_beta_age
rename sum_stats_sc`s'11 beta_snp
rename sum_stats_sc`s'12 stder_beta_snp

save sum_stats_sc`s'_`k'_`m',replace

restore

}
}
}


*** END OF SCRIPT ***

