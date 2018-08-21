/**********************************************************************************/
/***	Scenario $causalscenario (1=A, 2=B, 3=C)								***/
/***	Modify the local parameters by changing the preamble file called		***/
/***	to modify the scenario.													***/
/**********************************************************************************/

set more off
clear
*set seed 67208113

/*Create blank data set*/
set obs 40000 //creates blank dataset with XXXXX observations
gen id = _n


/*Step 1: Set parameters*/
do preamble_$causalscenario.do

*specify prevalence of exposure
local pexp = $pexp 	

*parameters for mortality  
//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
local g1_0to1 =		$g1_0to1
local g1_1to5 = 	$g1_1to5
foreach x in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 {
	local g1_`x'to`=`x'+5' = ${g1_`x'to`=`x'+5'}
}

*effects of covariates on mortality risk
local g2 = $g2 //log(HR) for effect of U on death 
local g3 = $g3 //log(HR) for interaction effect of exposure & U on death
local g4 = $g4 //log(HR) for effect of stroke on death		

*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
local lambda_0to1 = 	$lambda_0to1
local lambda_1to5 = 	$lambda_1to5
foreach x in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 {
	local lambda_`x'to`=`x'+5' = ${lambda_`x'to`=`x'+5'}
}

/*baseline hazard of stroke (exp=0 whites), based on race- and age-specific stroke incidence rates in REGARDS 
(rates provided by Dr. George Howard in Dec 2016 that are updes of rates published in Howard Ann Neurol 2011)*/
foreach x in 45 50 55 60 65 70 75 80 85 90 {
	local stk_lambda_exp0_`x'to`=`x'+5' = ${stk_lambda_exp0_`x'to`=`x'+5'}
}

/*baseline hazard of stroke (exp=0 whites), based on race- and age-specific stroke incidence rates in REGARDS 
(rates provided by Dr. George Howard in Dec 2016 that are updes of rates published in Howard Ann Neurol 2011)*/
foreach x in 45 50 55 60 65 70 75 80 85 90 {
	local stk_lambda_exp1_`x'to`=`x'+5' = ${stk_lambda_exp1_`x'to`=`x'+5'}
}

*parameter for stroke hazard
local b1 = $b1 			//log (HR) for U on stroke

*probability of death at stroke
local pstrokedeath = $pstrokedeath 


/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
gen exposure = runiform()<`pexp' 
/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
gen U = rnormal(0,1)


/*Step 4: Generate survival time (Tij) for each person and strokes for people alive
at each interval. 
a. Each person's underlying time to death is generated for each age interval, 
conditional on the past provided the person has not died in a previous interval, 
under an exponential survival distribtion. If the person's generated survival 
time exceeds the length of the interval between study visits j and j+1, 
she is considered alive at study visit j+1 and a new survival time is 
generated for the next interval conditional on history up to the start of the 
interval. The process is repeated until the person's survival time falls 
within a given interval or the end of the study, whichever comes first. Each 
person's mortality hazard function is defined as:
h(tij|x) = lambda*exp(g1*exposurei + g2*Uij + g3*exposurei*Uij + g4*stroke_historyij)
A person's survival time for a given time interval at risk is generated using 
the inverse cumulative hazard function transformation formula described by 
Bender et al. (Stat Med 2011)
b. Stroke code is adapted from survival time code. The effect of exposure on stroke
risk is set to a constant incidence rate difference (see preamble files): 
stk_lambda_exp1_agejtoagej+1 = stk_lambda_exp0_agejtoagej+1 + stk_lambda_delta*/

*ia. Generate uniform random variable for generating survival time
gen RV_0to1 = runiform()
gen RV_1to5 = runiform()
foreach x in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 {
	gen RV_`x'to`=`x'+5' = runiform()
}

*ib. Generate uniform random variable for generating stroke time
foreach x in 45 50 55 60 65 70 75 80 85 90 {
	gen RV2_`x'to`=`x'+5' = runiform()
}


*ii. Generate survival time and stroke time for each interval

gen survage = .
gen strokeage = .

/***Ages 0-45: no strokes, so only need to generate survival***/

/*Interval 0-1*/
*Generate survival time from time 0
gen survtime0to1 = -ln(RV_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
			+`g3'*exposure*U +`g4'*0))
*Generate death indicator for interval 0-1 
gen death0to1 = 0
replace death0to1 = 1 if (survtime0to1 < 1) 
replace survage = survtime0to1 if (death0to1==1)
*Generate indicator for death before age 1
gen death1 = death0to1


/*Interval 1-5*/
*Generate survival time from time 1
gen survtime1to5 = -ln(RV_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
			+`g3'*exposure*U +`g4'*0)) ///
				if (death0to1==0) 
*Generate death indicator for interval 1-5 
gen death1to5 = 0 if (death1==0)
replace death1to5 = 1 if (survtime1to5 < 5)
replace survage = 1 + survtime1to5 if (death1to5==1)
*Generate indicator for death before age 5
gen death5 = 0
replace death5 = 1 if survage < 5


/*Use loop to generate survival for 5-45*/
foreach x in 5 10 15 20 25 30 35 40 {
	/*Generate survival time fom age x to x+5*/
	gen survtime`x'to`=`x'+5' = -ln(RV_`x'to`=`x'+5')/(`lambda_`x'to`=`x'+5''*exp(`g1_`x'to`=`x'+5''*exposure +`g2'*U ///
				+`g3'*exposure*U +`g4'*0)) if (death`x'==0)
	/*Generate death indicator for interval x to x+5*/
	gen death`x'to`=`x'+5' = 0 if (death`x'==0)
	replace death`x'to`=`x'+5'= 1 if (survtime`x'to`=`x'+5' < 5) 
	replace survage = `x' + survtime`x'to`=`x'+5' if (death`x'to`=`x'+5'==1)
	/*Generate indicator for death before age x+5*/
	gen death`=`x'+5' = 0
	replace death`=`x'+5' = 1 if survage < `=`x'+5'
}


/***Starting at age 45--people are at risk of stroke, and prevalent stroke 
increases mortality risk, so we need to start iteratively generating deaths and strokes.
Prevalent stroke for age interval age j to age j+1 is a stroke that occured prior to 
age j. Since strokes first occurs starting at age 45, there are no prevalent strokes at
age 45. Thus, the effect of prevalent stroke on mortality is set to 0 for this interval:
	survtime45to50: `g4'*0
	For subsequent age intervals: `g4'*stroke_history***/

/*Generate prevalent stroke variable (prior to age 45, stroke_history=0 for everyone),
values replaced as strokes occur after age 45*/
gen stroke_history = 0

/*Generate variable for stroke up to age 45 (set to 0 for everyone, but creating 
variable because stroke up to age x is used for generating stroke between ages x and x+5*/
gen stroke45 = 0


/*Use loop to generate survival and stroke for 45-95*/
foreach x in 45 50 55 60 65 70 75 80 85 90 {
***a. Survival
	/*Generate survival time fom age x to x+5*/
	gen survtime`x'to`=`x'+5' = -ln(RV_`x'to`=`x'+5')/(`lambda_`x'to`=`x'+5''*exp(`g1_`x'to`=`x'+5''*exposure +`g2'*U ///
				+`g3'*exposure*U +`g4'*stroke_history)) if (death`x'==0)
	/*Generate death indicator for interval x to x+5*/
	gen death`x'to`=`x'+5' = 0 if (death`x'==0)
	replace death`x'to`=`x'+5'= 1 if (survtime`x'to`=`x'+5' < 5) 
	replace survage = `x' + survtime`x'to`=`x'+5' if (death`x'to`=`x'+5'==1)
	/*Generate indicator for death before age x+5*/
	gen death`=`x'+5' = 0
	replace death`=`x'+5' = 1 if survage < `=`x'+5'
	
***b. Stroke
	/*Generate stroke time from age x for exp==0*/
	gen stroketime`x'to`=`x'+5' = -ln(RV2_`x'to`=`x'+5')/(`stk_lambda_exp0_`x'to`=`x'+5''*exp(`b1'*U)) ///
					if (exp==0 & death`x'==0 & stroke`x'==0)
	/*Generate stroke time from age x for exp==1*/
	replace stroketime`x'to`=`x'+5' = -ln(RV2_`x'to`=`x'+5')/(`stk_lambda_exp1_`x'to`=`x'+5''*exp(`b1'*U)) ///
					if (exp==1 & death`x'==0 & stroke`x'==0)
	/*Generate stroke indicator for interval x to x+5*/
	gen stroke`x'to`=`x'+5' = 0 if (death`x'==0 & stroke`x'==0)
	replace stroke`x'to`=`x'+5' = 1 if (stroketime`x'to`=`x'+5' < 5) 
	replace stroke`x'to`=`x'+5' = 0 if (death`x'to`=`x'+5'==1 & stroketime`x'to`=`x'+5' != . & stroketime`x'to`=`x'+5' > survtime`x'to`=`x'+5')
	/*Update prevalent stroke variable (stroke up to age x+5)*/
	replace stroke_history = 1 if (stroke`x'to`=`x'+5'==1)
	replace strokeage = `x' + stroketime`x'to`=`x'+5' if (stroke`x'to`=`x'+5'==1)
	/*Generate indicator for stroke before age x+5*/						
	gen stroke`=`x'+5' = 0
	replace stroke`=`x'+5' = 1 if strokeage < `=`x'+5'

	/*Generate stroke deaths for interval x to x+5, and replace death and survage accordingly if stroke death=1*/																
	gen strokedeath`x'to`=`x'+5' = runiform()<`pstrokedeath' if stroke`x'to`=`x'+5'==1
	replace death`x'to`=`x'+5' = 1 if (strokedeath`x'to`=`x'+5'==1) 
	replace survage = strokeage if (strokedeath`x'to`=`x'+5'==1) 
	/*Generate indicator for death before age x+5*/
	replace death`=`x'+5' = 1 if strokedeath`x'to`=`x'+5'==1	
}

replace survage = 90 + survtime90to95 if (survtime90to95!=. & death95==0) 	//this line is what differs from survtimes at younger intervals
	

*everyone dies
gen death = 1


*iv. Generate variables for stroke and stroke death between ages 45 to 95
gen stroke = 0
replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | ///
stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80==1 | ///
stroke80to85==1 | stroke85to90==1 | stroke90to95==1)
gen strokedeath = 0
replace strokedeath = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
strokedeath80to85==1 | strokedeath85to90==1 | strokedeath90to95==1)


*v. Generate strokeage for people who didn't develop a stroke
replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45


*vi. Generate age-stratified stroke variables
foreach x in 45 55 65 75 85 {
	gen strokeage`x'to`=`x'+10'_start = `x' if (strokeage !=. & strokeage>=`x')
	gen strokeage`x'to`=`x'+10'_end = `=`x'+9.99999' if (strokeage !=. & strokeage>=`=`x'+10')
	replace strokeage`x'to`=`x'+10'_end = strokeage if (strokeage !=. & strokeage>=`x' & strokeage<`=`x'+10')
	gen strokeage`x'to`=`x'+10'_contributor = 1 if (strokeage !=. & strokeage>=`x')
	replace strokeage`x'to`=`x'+10'_contributor = 0 if (strokeage !=. & strokeage<`x')
	gen stroke`x'to`=`x'+10' = .
	replace stroke`x'to`=`x'+10' = 0 
	replace stroke`x'to`=`x'+10' = 1 if (stroke`x'to`=`x'+5'==1 | stroke`=`x'+5'to`=`x'+10'==1)
}
/**************************************************************/

/**********************************************************************************/
/***	End data generation														***/
/**********************************************************************************/


/**********************************************************************************/
/***	Select a sample of observations still alive at start of study			***/
/**********************************************************************************/


/**********************************************************************************/
/***	MODELS																	***/
/**********************************************************************************/

/******************************************/
*pull N
scalar N= _N

*pull N in each exposure group
summarize exposure
scalar N_exp1 = r(mean)*r(N)
scalar N_exp0 = r(N) - r(mean)*r(N)
/******************************************/

/******************************************/
*pull proportion of deaths at start and end of stroke f/u
*age 45, 50, 55, ..., 95
foreach x in 45 50 55 60 65 70 75 80 85 90 95 {
	summarize death`x', meanonly
	scalar p_death`x' = r(mean)
	summarize death`x' if (exposure==1), meanonly
	scalar p_death`x'_exp1 = r(mean)
	summarize death`x' if (exposure==0), meanonly
	scalar p_death`x'_exp0 = r(mean)
}


*pull median survival time
summarize survage, detail 
scalar med_survage = r(p50)
qui summarize survage if (exposure==1), detail
scalar med_survage_exp1 = r(p50)
qui summarize survage if (exposure==0), detail
scalar med_survage_exp0 = r(p50)


*pull percentage of strokes
summarize stroke, meanonly 
scalar p_stroke = r(mean)
summarize stroke if (exposure==1), meanonly 
scalar p_stroke_exp1 = r(mean)
summarize stroke if (exposure==0), meanonly 
scalar p_stroke_exp0 = r(mean)
/******************************************/

/******************************************/
*distribution of U among people at risk for stroke in each age group at birth and in 5-year intervals starting at age 45
*U and race should be independent at birth
*birth
qui sum U if (exposure==1)
scalar meanUatrisk0_exp1= r(mean)
qui sum U if (exposure==0)
scalar meanUatrisk0_exp0= r(mean)

*age 45, 50, 55, ..., 95
foreach x in 45 50 55 60 65 70 75 80 85 90 95 {
	qui sum U if (exposure==1 & strokeage !=. & strokeage>=`x')
	scalar meanUatrisk`x'_exp1= r(mean)
	scalar Natrisk`x'_exp1= r(N)
	qui sum U if (exposure==0 & strokeage !=. & strokeage>=`x')
	scalar meanUatrisk`x'_exp0= r(mean)
	scalar Natrisk`x'_exp0= r(N)
}
/******************************************/

/******************************************/
*age-stratified incidence rates by exposure, incidence rate differences, and incidence rate ratios
foreach x in 45to55 55to65 65to75 75to85 85to95 {
	qui stset strokeage, failure(stroke`x'==1) id(id) enter(strokeage`x'_start) exit(strokeage`x'_end)
	*stroke IR for blacks
	qui stptime if (exposure==1), title(person-years) per(10000)
	scalar nstrokes`x'_exp1 = r(failures)
	scalar ptime`x'_exp1 = r(ptime)
	scalar strokerate`x'_exp1 = r(rate)
	*stroke IR for whites
	qui stptime if (exposure==0), title(person-years) per(10000)
	scalar nstrokes`x'_exp0 = r(failures)
	scalar ptime`x'_exp0 = r(ptime)
	scalar strokerate`x'_exp0 = r(rate)
	*stroke IRR and IRD for blacks vs. whites
	qui stir exposure
	scalar strokeIRR`x' = r(irr)
	scalar strokelnIRR`x' = ln(r(irr))
	scalar strokelnIRR`x'_SE = sqrt((1/nstrokes`x'_exp1)+(1/nstrokes`x'_exp0)) 
	scalar strokeIRR`x'_ub = r(ub_irr)
	scalar strokeIRR`x'_lb = r(lb_irr)
	scalar strokeIRD`x' = r(ird)*10000
	scalar strokeIRD`x'_SE = sqrt((nstrokes`x'_exp1/((ptime`x'_exp1)^2))+(1/nstrokes`x'_exp0/((ptime`x'_exp0)^2)))*10000 
	scalar strokeIRD`x'_ub = r(ub_ird)*10000
	scalar strokeIRD`x'_lb = r(lb_ird)*10000
}
/******************************************/

/******************************************/
*pull n strokes per 5-year intervals
*age 45, 50, 55, ..., 95
foreach x in 45to50 50to55 55to60 60to65 65to70 70to75 75to80 80to85 85to90 90to95 {
	sum stroke`x' if exposure==0
	scalar p_stroke`x'_exp0 = r(mean)
	scalar nstrokes`x'_exp0 = r(mean)*r(N)

	sum stroke`x' if exposure==1
	scalar p_stroke`x'_exp1 = r(mean)
	scalar nstrokes`x'_exp1 = r(mean)*r(N)
}
/******************************************/

