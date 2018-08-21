/******************************************************************************/
/***	Scenario 1/A (U influences stroke but not mortality).				***/
/******************************************************************************/

/*Step 1: Set parameters*/

*specify prevalence of exposure
global pexp = 0.5 	

*parameters for death
//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
*values for ages 0-45 are calculated directly from US life tables
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
*values for ages 45+ identified from search loop because strokes begin occurring, stroke increase risk of death, and exposure increases stroke risk
global g1_45to50 = 	.622
global g1_50to55 = 	.58
global g1_55to60 = 	.426
global g1_60to65 = 	.29
global g1_65to70 = 	.186
global g1_70to75 = 	.054
global g1_75to80 = 	-.052
global g1_80to85 = 	-.132
global g1_85to90 = 	-.204
global g1_90to95 = 	-.264

*effects of covariates on mortality risk
global g2 = 0		  //log(HR) for effect of U on death 
global g3 = 0 		//log(HR) for interaction effect of exposure & U on death
global g4 = ln(2)	//log(HR) for effect of stroke history on death

*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
global lambda_0to1 = 	0.0749
global lambda_1to5 = 	0.0082
global lambda_5to10 = 	0.0028
global lambda_10to15 = 	0.0021
global lambda_15to20 = 	0.0034
global lambda_20to25 = 	0.0048
global lambda_25to30 = 	0.0056
global lambda_30to35 = 	0.0062
global lambda_35to40 = 	0.0068
global lambda_40to45 = 	0.0078
*values for ages 45+ identified from search loop because strokes begin occuring, stroke increase risk of death
global lambda_45to50 = 	.0102
global lambda_50to55 = 	.0132
global lambda_55to60 =	.0188
global lambda_60to65 = 	.0262
global lambda_65to70 =  .0388
global lambda_70to75 = 	.0606
global lambda_75to80 = 	.0918
global lambda_80to85 = 	.142
global lambda_85to90 = 	.214
global lambda_90to95 = 	.3

/*baseline hazard of stroke (exp=0 whites)*/
*values identified from search loop 	
global stk_lambda_exp0_45to50 = .00054
global stk_lambda_exp0_50to55 = .00108
global stk_lambda_exp0_55to60 = .0018
global stk_lambda_exp0_60to65 = .00334
global stk_lambda_exp0_65to70 = .00456
global stk_lambda_exp0_70to75 = .00724
global stk_lambda_exp0_75to80 = .00968
global stk_lambda_exp0_80to85 = .0118
global stk_lambda_exp0_85to90 = .014
global stk_lambda_exp0_90to95 = .0178

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

