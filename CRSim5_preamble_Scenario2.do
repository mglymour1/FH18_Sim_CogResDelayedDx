/******************************************************************************/
/***	Scenario 2/B (U influences stroke but and mortality).				***/
/******************************************************************************/

/*Step 1: Set parameters*/

*specify prevalence of exposure
global pexp = 0.5 	

*parameters for death
//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
global g1_0to1 =	0.30
global g1_1to5 = 	0.37
global g1_5to10 = 	0.23
global g1_10to15 = 	0.56
global g1_15to20 = 	0.94
global g1_20to25 = 	0.92
global g1_25to30 = 	0.79
global g1_30to35 = 	0.78
global g1_35to40 = 	0.77
global g1_40to45 = 	0.74
*values for ages 45+ identified from search loop because strokes begin occurring, stroke increase risk of death
global g1_45to50 = 	.642
global g1_50to55 = 	.63599998
global g1_55to60 = 	.51799998
global g1_60to65 = 	.34800001
global g1_65to70 = 	.29
global g1_70to75 = 	.128
global g1_75to80 = 	.002
global g1_80to85 = 	-.066
global g1_85to90 = 	-.184
global g1_90to95 = 	-.236

*effects of covariates on mortality risk
global g2 = ln(1.5)	//log(HR) for effect of U on death 
global g3 = 0 		  //log(HR) for interaction effect of exposure & U on death
global g4 = ln(2)	  //log(HR) for effect of stroke history on death

*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
*values identified from search loop 
global lambda_0to1 = 	.0702
global lambda_1to5 = 	.0066
global lambda_5to10 = 	.003
global lambda_10to15 = 	.0024
global lambda_15to20 = 	.0038
global lambda_20to25 = 	.005
global lambda_25to30 = 	.0058
global lambda_30to35 = 	.0068
global lambda_35to40 = 	.0072
global lambda_40to45 = 	.0078
global lambda_45to50 =  .0096
global lambda_50to55 = 	.0118 
global lambda_55to60 =	.0164 
global lambda_60to65 = 	.0238 
global lambda_65to70 =  .0364 
global lambda_70to75 = 	.0606 
global lambda_75to80 = 	.1058 
global lambda_80to85 = 	.17600001 
global lambda_85to90 = 	.294 
global lambda_90to95 = 	.47

/*baseline hazard of stroke (exp=0 whites)*/ 
*baseline hazard of stroke identified from search loop because U influences risk of death from birth in this scenario
global stk_lambda_exp0_45to50 = .00054
global stk_lambda_exp0_50to55 = .00116
global stk_lambda_exp0_55to60 = .00194
global stk_lambda_exp0_60to65 = .0035
global stk_lambda_exp0_65to70 = .0051
global stk_lambda_exp0_70to75 = .0081
global stk_lambda_exp0_75to80 = .01226
global stk_lambda_exp0_80to85 = .0164
global stk_lambda_exp0_85to90 = .0212
global stk_lambda_exp0_90to95 = .029

/*baseline hazard of stroke (exp=1 blacks)*/
global stk_lambda_delta = 0.002 //age-constant incidence rate difference for exposure=1 vs. exposure=0
global stk_lambda_exp1_45to50 = 	$stk_lambda_exp0_45to50 + $stk_lambda_delta
global stk_lambda_exp1_50to55 = 	$stk_lambda_exp0_50to55 + $stk_lambda_delta
global stk_lambda_exp1_55to60 =		$stk_lambda_exp0_55to60 + $stk_lambda_delta
global stk_lambda_exp1_60to65 = 	$stk_lambda_exp0_60to65 + $stk_lambda_delta
global stk_lambda_exp1_65to70 =  	$stk_lambda_exp0_65to70 + $stk_lambda_delta
global stk_lambda_exp1_70to75 = 	$stk_lambda_exp0_70to75 + $stk_lambda_delta
global stk_lambda_exp1_75to80 = 	$stk_lambda_exp0_75to80 + $stk_lambda_delta
global stk_lambda_exp1_80to85 = 	$stk_lambda_exp0_80to85 + $stk_lambda_delta
global stk_lambda_exp1_85to90 = 	$stk_lambda_exp0_85to90 + $stk_lambda_delta 
global stk_lambda_exp1_90to95 = 	$stk_lambda_exp0_90to95 + $stk_lambda_delta

*effects of covariates on stroke risk
global b1 = ln(1.5)		//log (HR) for U on stroke

*probability of death at stroke
global pstrokedeath = 0.25 

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

