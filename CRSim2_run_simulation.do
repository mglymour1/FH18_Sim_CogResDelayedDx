* This code sets up all the local macro variables for the mortality rates
* calls the actual simulate code with the requested input macros
* and summarizes simulation findings.

set more off
clear all

di "$S_TIME"

timer clear 1

timer on 1

capture log close


local B = $B //number of simulations
/*other global vars called in this program:
outputrow //output row for summarized simulation results
causalscenario //where 1=A, 2=B, 3=C*/


/*Pull in scalars from data generation and analysis file*/
local simlist ""
foreach x in N N_exp1 N_exp0 ///
	p_death45 p_death45_exp0 p_death45_exp1 ///
	p_death50 p_death50_exp0 p_death50_exp1 ///
	p_death55 p_death55_exp0 p_death55_exp1 ///
	p_death60 p_death60_exp0 p_death60_exp1 ///
	p_death65 p_death65_exp0 p_death65_exp1 ///
	p_death70 p_death70_exp0 p_death70_exp1 ///
	p_death75 p_death75_exp0 p_death75_exp1 ///
	p_death80 p_death80_exp0 p_death80_exp1 ///
	p_death85 p_death85_exp0 p_death85_exp1 ///
	p_death90 p_death90_exp0 p_death90_exp1 ///
	p_death95 p_death95_exp0 p_death95_exp1 ///
	med_survage med_survage_exp0 med_survage_exp1 ///
	p_stroke p_stroke_exp0 p_stroke_exp1 ///
	strokerate45to55_exp0 strokerate45to55_exp1 ///
	strokerate55to65_exp0 strokerate55to65_exp1 ///
	strokerate65to75_exp0 strokerate65to75_exp1 ///
	strokerate75to85_exp0 strokerate75to85_exp1 ///
	strokerate85to95_exp0 strokerate85to95_exp1 ///
	strokeIRR45to55 strokeIRR45to55_lb strokeIRR45to55_ub strokelnIRR45to55 strokelnIRR45to55_SE ///
	strokeIRR55to65 strokeIRR55to65_lb strokeIRR55to65_ub strokelnIRR55to65 strokelnIRR55to65_SE ///
	strokeIRR65to75 strokeIRR65to75_lb strokeIRR65to75_ub strokelnIRR65to75 strokelnIRR65to75_SE ///
	strokeIRR75to85 strokeIRR75to85_lb strokeIRR75to85_ub strokelnIRR75to85 strokelnIRR75to85_SE ///
	strokeIRR85to95 strokeIRR85to95_lb strokeIRR85to95_ub strokelnIRR85to95 strokelnIRR85to95_SE ///
	strokeIRD45to55 strokeIRD45to55_lb strokeIRD45to55_ub strokeIRD45to55_SE ///
	strokeIRD55to65 strokeIRD55to65_lb strokeIRD55to65_ub strokeIRD55to65_SE ///
	strokeIRD65to75 strokeIRD65to75_lb strokeIRD65to75_ub strokeIRD65to75_SE ///
	strokeIRD75to85 strokeIRD75to85_lb strokeIRD75to85_ub strokeIRD75to85_SE ///
	strokeIRD85to95 strokeIRD85to95_lb strokeIRD85to95_ub strokeIRD85to95_SE ///
	nstrokes45to55_exp1 ptime45to55_exp1 ///
	nstrokes45to55_exp0 ptime45to55_exp0 ///
	nstrokes55to65_exp1 ptime55to65_exp1 ///
	nstrokes55to65_exp0 ptime55to65_exp0 ///
	nstrokes65to75_exp1 ptime65to75_exp1 ///
	nstrokes65to75_exp0 ptime65to75_exp0 ///
	nstrokes75to85_exp1 ptime75to85_exp1 ///
	nstrokes75to85_exp0 ptime75to85_exp0 ///
	nstrokes85to95_exp1 ptime85to95_exp1 ///
	nstrokes85to95_exp0 ptime85to95_exp0 ///
	meanUatrisk0_exp1 ///
	meanUatrisk0_exp0 ///
	Natrisk45_exp1 meanUatrisk45_exp1 ///
	Natrisk45_exp0 meanUatrisk45_exp0 ///
	Natrisk50_exp1 meanUatrisk50_exp1 ///
	Natrisk50_exp0 meanUatrisk50_exp0 ///
	Natrisk55_exp1 meanUatrisk55_exp1 ///
	Natrisk55_exp0 meanUatrisk55_exp0 ///
	Natrisk60_exp1 meanUatrisk60_exp1 ///
	Natrisk60_exp0 meanUatrisk60_exp0 ///
	Natrisk65_exp1 meanUatrisk65_exp1 ///
	Natrisk65_exp0 meanUatrisk65_exp0 ///
	Natrisk70_exp1 meanUatrisk70_exp1 ///
	Natrisk70_exp0 meanUatrisk70_exp0 ///
	Natrisk75_exp1 meanUatrisk75_exp1 ///
	Natrisk75_exp0 meanUatrisk75_exp0 ///
	Natrisk80_exp1 meanUatrisk80_exp1 ///
	Natrisk80_exp0 meanUatrisk80_exp0 ///
	Natrisk85_exp1 meanUatrisk85_exp1 ///
	Natrisk85_exp0 meanUatrisk85_exp0 ///
	Natrisk90_exp1 meanUatrisk90_exp1 ///
	Natrisk90_exp0 meanUatrisk90_exp0 ///
	nstrokes45to50_exp0 nstrokes50to55_exp0 ///
	nstrokes55to60_exp0 nstrokes60to65_exp0 ///
	nstrokes65to70_exp0 nstrokes70to75_exp0 ///
	nstrokes75to80_exp0 nstrokes80to85_exp0 ///
	nstrokes85to90_exp0 nstrokes90to95_exp0 ///
	nstrokes45to50_exp1 nstrokes50to55_exp1 ///
	nstrokes55to60_exp1 nstrokes60to65_exp1 ///
	nstrokes65to70_exp1 nstrokes70to75_exp1 ///
	nstrokes75to80_exp1 nstrokes80to85_exp1 ///
	nstrokes85to90_exp1 nstrokes90to95_exp1 ///
	p_stroke45to50_exp0 p_stroke50to55_exp0 ///
	p_stroke55to60_exp0 p_stroke60to65_exp0 ///
	p_stroke65to70_exp0 p_stroke70to75_exp0 ///
	p_stroke75to80_exp0 p_stroke80to85_exp0 ///
	p_stroke85to90_exp0 p_stroke90to95_exp0 ///
	p_stroke45to50_exp1 p_stroke50to55_exp1 ///
	p_stroke55to60_exp1 p_stroke60to65_exp1 ///
	p_stroke65to70_exp1 p_stroke70to75_exp1 ///
	p_stroke75to80_exp1 p_stroke80to85_exp1 ///
	p_stroke85to90_exp1 p_stroke90to95_exp1 {  
   local simlist "`simlist' `x'=`x'"
}

/*Run simulation*/
simulate `simlist', ///
reps(`B') seed(67208113): do data_generation_analysis

*number of individuals
local N = N
local N_rounded = round(N,100)

*generate p alive variables
gen p_alive45 = 1 - p_death45 
gen p_alive45_exp0 = 1 - p_death45_exp0 
gen p_alive45_exp1 = 1 - p_death45_exp1

gen p_alive50 = 1 - p_death50 
gen p_alive50_exp0 = 1 - p_death50_exp0 
gen p_alive50_exp1 = 1 - p_death50_exp1

gen p_alive55 = 1 - p_death55 
gen p_alive55_exp0 = 1 - p_death55_exp0 
gen p_alive55_exp1 = 1 - p_death55_exp1

gen p_alive60 = 1 - p_death60 
gen p_alive60_exp0 = 1 - p_death60_exp0 
gen p_alive60_exp1 = 1 - p_death60_exp1

gen p_alive65 = 1 - p_death65 
gen p_alive65_exp0 = 1 - p_death65_exp0 
gen p_alive65_exp1 = 1 - p_death65_exp1

gen p_alive70 = 1 - p_death70 
gen p_alive70_exp0 = 1 - p_death70_exp0 
gen p_alive70_exp1 = 1 - p_death70_exp1

gen p_alive75 = 1 - p_death75 
gen p_alive75_exp0 = 1 - p_death75_exp0 
gen p_alive75_exp1 = 1 - p_death75_exp1

gen p_alive80 = 1 - p_death80 
gen p_alive80_exp0 = 1 - p_death80_exp0 
gen p_alive80_exp1 = 1 - p_death80_exp1

gen p_alive85 = 1 - p_death85 
gen p_alive85_exp0 = 1 - p_death85_exp0 
gen p_alive85_exp1 = 1 - p_death85_exp1

gen p_alive90 = 1 - p_death90 
gen p_alive90_exp0 = 1 - p_death90_exp0 
gen p_alive90_exp1 = 1 - p_death90_exp1

gen p_alive95 = 1 - p_death95 
gen p_alive95_exp0 = 1 - p_death95_exp0 
gen p_alive95_exp1 = 1 - p_death95_exp1


*scalars for means
foreach x in N_exp1 N_exp0 ///
	p_death45 p_death45_exp0 p_death45_exp1 ///
	p_alive45 p_alive45_exp0 p_alive45_exp1 ///
	p_alive50 p_alive50_exp0 p_alive50_exp1 ///
	p_alive55 p_alive55_exp0 p_alive55_exp1 ///
	p_alive60 p_alive60_exp0 p_alive60_exp1 ///
	p_alive65 p_alive65_exp0 p_alive65_exp1 ///
	p_alive70 p_alive70_exp0 p_alive70_exp1 ///
	p_alive75 p_alive75_exp0 p_alive75_exp1 ///
	p_alive80 p_alive80_exp0 p_alive80_exp1 ///
	p_alive85 p_alive85_exp0 p_alive85_exp1 ///
	p_alive90 p_alive90_exp0 p_alive90_exp1 ///
	p_alive95 p_alive95_exp0 p_alive95_exp1 ///
	med_survage med_survage_exp0 med_survage_exp1 ///
	p_stroke p_stroke_exp0 p_stroke_exp1 ///
	strokerate45to55_exp0 strokerate45to55_exp1 ///
	strokerate55to65_exp0 strokerate55to65_exp1 ///
	strokerate65to75_exp0 strokerate65to75_exp1 ///
	strokerate75to85_exp0 strokerate75to85_exp1 ///
	strokerate85to95_exp0 strokerate85to95_exp1 ///
	strokeIRR45to55 ///
	strokeIRR55to65 ///
	strokeIRR65to75 ///
	strokeIRR75to85 ///
	strokeIRR85to95 ///
	strokeIRD45to55 ///
	strokeIRD55to65 ///
	strokeIRD65to75 ///
	strokeIRD75to85 ///
	strokeIRD85to95 ///
	nstrokes45to55_exp0 ptime45to55_exp0 ///
	nstrokes45to55_exp1 ptime45to55_exp1 ///
	nstrokes55to65_exp0 ptime55to65_exp0 ///
	nstrokes55to65_exp1 ptime55to65_exp1 ///
	nstrokes65to75_exp0 ptime65to75_exp0 ///
	nstrokes65to75_exp1 ptime65to75_exp1 ///
	nstrokes75to85_exp0 ptime75to85_exp0 ///
	nstrokes75to85_exp1 ptime75to85_exp1 ///
	nstrokes85to95_exp0 ptime85to95_exp0 ///
	nstrokes85to95_exp1 ptime85to95_exp1 ///
	meanUatrisk0_exp1 ///
	meanUatrisk0_exp0 ///
	Natrisk45_exp1 meanUatrisk45_exp1 ///
	Natrisk45_exp0 meanUatrisk45_exp0 ///
	Natrisk50_exp1 meanUatrisk50_exp1 ///
	Natrisk50_exp0 meanUatrisk50_exp0 ///
	Natrisk55_exp1 meanUatrisk55_exp1 ///
	Natrisk55_exp0 meanUatrisk55_exp0 ///
	Natrisk60_exp1 meanUatrisk60_exp1 ///
	Natrisk60_exp0 meanUatrisk60_exp0 ///
	Natrisk65_exp1 meanUatrisk65_exp1 ///
	Natrisk65_exp0 meanUatrisk65_exp0 ///
	Natrisk70_exp1 meanUatrisk70_exp1 ///
	Natrisk70_exp0 meanUatrisk70_exp0 ///
	Natrisk75_exp1 meanUatrisk75_exp1 ///
	Natrisk75_exp0 meanUatrisk75_exp0 ///
	Natrisk80_exp1 meanUatrisk80_exp1 ///
	Natrisk80_exp0 meanUatrisk80_exp0 ///
	Natrisk85_exp1 meanUatrisk85_exp1 ///
	Natrisk85_exp0 meanUatrisk85_exp0 ///
	Natrisk90_exp1 meanUatrisk90_exp1 ///
	Natrisk90_exp0 meanUatrisk90_exp0 ///
	nstrokes45to50_exp0 nstrokes50to55_exp0 ///
	nstrokes55to60_exp0 nstrokes60to65_exp0 ///
	nstrokes65to70_exp0 nstrokes70to75_exp0 ///
	nstrokes75to80_exp0 nstrokes80to85_exp0 ///
	nstrokes85to90_exp0 nstrokes90to95_exp0 ///
	nstrokes45to50_exp1 nstrokes50to55_exp1 ///
	nstrokes55to60_exp1 nstrokes60to65_exp1 ///
	nstrokes65to70_exp1 nstrokes70to75_exp1 ///
	nstrokes75to80_exp1 nstrokes80to85_exp1 ///
	nstrokes85to90_exp1 nstrokes90to95_exp1 ///
	p_stroke45to50_exp0 p_stroke50to55_exp0 ///
	p_stroke55to60_exp0 p_stroke60to65_exp0 ///
	p_stroke65to70_exp0 p_stroke70to75_exp0 ///
	p_stroke75to80_exp0 p_stroke80to85_exp0 ///
	p_stroke85to90_exp0 p_stroke90to95_exp0 ///
	p_stroke45to50_exp1 p_stroke50to55_exp1 ///
	p_stroke55to60_exp1 p_stroke60to65_exp1 ///
	p_stroke65to70_exp1 p_stroke70to75_exp1 ///
	p_stroke75to80_exp1 p_stroke80to85_exp1 ///
	p_stroke85to90_exp1 p_stroke90to95_exp1 ///
	strokelnIRR45to55 strokelnIRR55to65 ///
	strokelnIRR65to75 strokelnIRR75to85 ///
	strokelnIRR85to95 ///
	strokelnIRR45to55_SE strokeIRD45to55_SE ///
	strokelnIRR55to65_SE strokeIRD55to65_SE ///
	strokelnIRR65to75_SE strokeIRD65to75_SE ///
	strokelnIRR75to85_SE strokeIRD75to85_SE ///
	strokelnIRR85to95_SE strokeIRD85to95_SE {
summarize `x', meanonly
scalar mean_`x' = round(r(mean),0.001)
}


*scalars for empirical SEs
foreach x in strokerate45to55_exp0 strokerate45to55_exp1 ///
	strokelnIRR45to55  strokeIRD45to55 ///
	strokerate55to65_exp0 strokerate55to65_exp1 ///
	strokelnIRR55to65 strokeIRD55to65 ///
	strokerate65to75_exp0 strokerate65to75_exp1 ///
	strokelnIRR65to75 strokeIRD65to75  ///
	strokerate75to85_exp0 strokerate75to85_exp1 ///
	strokelnIRR75to85 strokeIRD75to85  ///
	strokerate85to95_exp0 strokerate85to95_exp1 ///
	strokelnIRR85to95 strokeIRD85to95  {
summarize `x'
scalar empSE_`x' = round(r(sd),0.001)
}


/**************************************************************************************************/
/***	summarize simulation results															***/
/***	results summarized across different Excel sheets										***/
/**************************************************************************************************/
/***sheet 1: IRR results***/
putexcel A1=("Scenario") ///
B1=("IRR 45to55") C1=("strokelnIRR45to55") D1=("empSE(lnIRR 45to55)") E1=("avgSE(lnIRR 45to55)") ///
F1=("IRR 55to65") G1=("strokelnIRR55to65") H1=("empSE(lnIRR 55to65)") I1=("avgSE(lnIRR 55to65)") ///
J1=("IRR 65to75") K1=("strokelnIRR65to75") L1=("empSE(lnIRR 65to75)") M1=("avgSE(lnIRR 65to75)") ///
N1=("IRR 75to85") O1=("strokelnIRR75to85") P1=("empSE(lnIRR 75to85)") Q1=("avgSE(lnIRR 75to85)") ///
R1=("IRR 85to95") S1=("strokelnIRR85to95") T1=("empSE(lnIRR 85to95)") U1=("avgSE(lnIRR 85to95)") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_strokeIRR45to55) C$outputrow=(mean_strokelnIRR45to55) D$outputrow=(empSE_strokelnIRR45to55) E$outputrow=(mean_strokelnIRR45to55_SE) ///
F$outputrow=(mean_strokeIRR55to65) G$outputrow=(mean_strokelnIRR55to65) H$outputrow=(empSE_strokelnIRR55to65) I$outputrow=(mean_strokelnIRR55to65_SE) ///
J$outputrow=(mean_strokeIRR65to75) K$outputrow=(mean_strokelnIRR65to75) L$outputrow=(empSE_strokelnIRR65to75) M$outputrow=(mean_strokelnIRR65to75_SE) ///
N$outputrow=(mean_strokeIRR75to85) O$outputrow=(mean_strokelnIRR75to85) P$outputrow=(empSE_strokelnIRR75to85) Q$outputrow=(mean_strokelnIRR75to85_SE) ///
R$outputrow=(mean_strokeIRR85to95) S$outputrow=(mean_strokelnIRR85to95) T$outputrow=(empSE_strokelnIRR85to95) U$outputrow=(mean_strokelnIRR85to95_SE) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Stroke_IRR") modify

/***sheet 2: IRD results***/
putexcel A1=("Scenario") ///
B1=("IRD 45to55") C1=("empSE(IRD 45to55)") D1=("avgSE(IRD 45to55)") ///
E1=("IRD 55to65") F1=("empSE(IRD 55to65)") G1=("avgSE(IRD 55to65)") ///
H1=("IRD 65to75") I1=("empSE(IRD 65to75)") J1=("avgSE(IRD 65to75)") ///
K1=("IRD 75to85") L1=("empSE(IRD 75to85)") M1=("avgSE(IRD 75to85)") ///
N1=("IRD 85to95") O1=("empSE(IRD 85to95)") P1=("avgSE(IRD 85to95)") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_strokeIRD45to55) C$outputrow=(empSE_strokeIRD45to55) D$outputrow=(mean_strokeIRD45to55_SE) ///
E$outputrow=(mean_strokeIRD55to65) F$outputrow=(empSE_strokeIRD55to65) G$outputrow=(mean_strokeIRD55to65_SE) ///
H$outputrow=(mean_strokeIRD65to75) I$outputrow=(empSE_strokeIRD65to75) J$outputrow=(mean_strokeIRD65to75_SE) ///
K$outputrow=(mean_strokeIRD75to85) L$outputrow=(empSE_strokeIRD75to85) M$outputrow=(mean_strokeIRD75to85_SE) ///
N$outputrow=(mean_strokeIRD85to95) O$outputrow=(empSE_strokeIRD85to95) P$outputrow=(mean_strokeIRD85to95_SE) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Stroke_IRD") modify

/***sheet 3: IR results per 1,000 person-years***/
putexcel A1=("Scenario") ///
B1=("IR/10000 PYs 45to55 whites") C1=("IR/10000 PYs 55to65 whites") ///
D1=("IR/10000 PYs 65to75 whites") E1=("IR/10000 PYs 75to85 whites") ///
F1=("IR/10000 PYs 85to95 whites") G1=("IR/10000 PYs 45to55 blacks") ///
H1=("IR/10000 PYs 55to65 blacks") I1=("IR/10000 PYs 65to75 blacks") ///
J1=("IR/10000 PYs 75to85 blacks") K1=("IR/1000 PYs 85to95 blacks") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_strokerate45to55_exp0) C$outputrow=(mean_strokerate55to65_exp0) ///
D$outputrow=(mean_strokerate65to75_exp0) E$outputrow=(mean_strokerate75to85_exp0) ///
F$outputrow=(mean_strokerate85to95_exp0) G$outputrow=(mean_strokerate45to55_exp1) ///
H$outputrow=(mean_strokerate55to65_exp1) I$outputrow=(mean_strokerate65to75_exp1) ///
J$outputrow=(mean_strokerate75to85_exp1) K$outputrow=(mean_strokerate85to95_exp1) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Stroke_IR10000py") modify

/***sheet 4: # people at risk for stroke and PYs***/
putexcel A1=("Scenario") ///
B1=("# strokes 45to55 whites") C1=("PYs 45to55 whites") ///
D1=("# strokes 55to65 whites") E1=("PYs 55to65 whites") ///
F1=("# strokes 65to75 whites") G1=("PYs 65to75 whites") ///
H1=("# strokes 75to85 whites") I1=("PYs 75to85 whites") ///
J1=("# strokes 85to95 whites") K1=("PYs 85to95 whites") ///
L1=("# strokes 45to55 blacks") M1=("PYs 45to55 blacks") ///
N1=("# strokes 55to65 blacks") O1=("PYs 55to65 blacks") ///
P1=("# strokes 65to75 blacks") Q1=("PYs 65to75 blacks") ///
R1=("# strokes 75to85 blacks") S1=("PYs 75to85 blacks") ///
T1=("# strokes 85to95 blacks") U1=("PYs 85to95 blacks") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_nstrokes45to55_exp0) C$outputrow=(mean_ptime45to55_exp0) ///
D$outputrow=(mean_nstrokes55to65_exp0) E$outputrow=(mean_ptime55to65_exp0) ///
F$outputrow=(mean_nstrokes65to75_exp0) G$outputrow=(mean_ptime65to75_exp0) ///
H$outputrow=(mean_nstrokes75to85_exp0) I$outputrow=(mean_ptime75to85_exp0) ///
J$outputrow=(mean_nstrokes85to95_exp0) K$outputrow=(mean_ptime85to95_exp0) ///
L$outputrow=(mean_nstrokes45to55_exp1) M$outputrow=(mean_ptime45to55_exp1) ///
N$outputrow=(mean_nstrokes55to65_exp1) O$outputrow=(mean_ptime55to65_exp1) ///
P$outputrow=(mean_nstrokes65to75_exp1) Q$outputrow=(mean_ptime65to75_exp1) ///
R$outputrow=(mean_nstrokes75to85_exp1) S$outputrow=(mean_ptime75to85_exp1) ///
T$outputrow=(mean_nstrokes85to95_exp1) U$outputrow=(mean_ptime85to95_exp1) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Nstks_PYs") modify

/***sheet 5: Pull death info***/
putexcel A1=("Scenario") ///
B1=("N whites at birth") ///
C1=("Median surv age whites") D1=("P alive age 45 whites") ///
E1=("P alive age 50 whites") F1=("P alive age 55 whites") ///
G1=("P alive age 60 whites") H1=("P alive age 65 whites") ///
I1=("P alive age 70 whites") J1=("P alive age 75 whites") ///
K1=("P alive age 80 whites") L1=("P alive age 85 whites") ///
M1=("P alive age 90 whites") N1=("P alive age 95 whites") ///
O1=("N blacks at birth") ///
P1=("Median surv age blacks") Q1=("P alive age 45 blacks") ///
R1=("P alive age 50 blacks") S1=("P alive age 55 blacks") ///
T1=("P alive age 60 blacks") U1=("P alive age 65 blacks") ///
V1=("P alive age 70 blacks") W1=("P alive age 75 blacks") ///
X1=("P alive age 80 blacks") Y1=("P alive age 85 blacks") ///
Z1=("P alive age 90 blacks") AA1=("P alive age 95 blacks") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_N_exp0) ///
C$outputrow=(mean_med_survage_exp0) D$outputrow=(mean_p_alive45_exp0) ///
E$outputrow=(mean_p_alive50_exp0) F$outputrow=(mean_p_alive55_exp0) ///
G$outputrow=(mean_p_alive60_exp0) H$outputrow=(mean_p_alive65_exp0) ///
I$outputrow=(mean_p_alive70_exp0) J$outputrow=(mean_p_alive75_exp0) ///
K$outputrow=(mean_p_alive80_exp0) L$outputrow=(mean_p_alive85_exp0) ///
M$outputrow=(mean_p_alive90_exp0) N$outputrow=(mean_p_alive95_exp0) ///
O$outputrow=(mean_N_exp1) ///
P$outputrow=(mean_med_survage_exp1) Q$outputrow=(mean_p_alive45_exp1) ///
R$outputrow=(mean_p_alive50_exp1) S$outputrow=(mean_p_alive55_exp1) ///
T$outputrow=(mean_p_alive60_exp1) U$outputrow=(mean_p_alive65_exp1) ///
V$outputrow=(mean_p_alive70_exp1) W$outputrow=(mean_p_alive75_exp1) ///
X$outputrow=(mean_p_alive80_exp1) Y$outputrow=(mean_p_alive85_exp1) ///
Z$outputrow=(mean_p_alive90_exp1) AA$outputrow=(mean_p_alive95_exp1) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Surv") modify

/***sheet 6: distributions of U and N at risk at each age***/
putexcel A1=("Scenario") ///
B1=("mean U at risk birth white") C1=("mean U at risk 45 white") D1=("mean U at risk 50 white") E1=("mean U at risk 55 white") ///
F1=("mean U at risk 60 white") G1=("mean U at risk 65 white") H1=("mean U at risk 70 white") ///
I1=("mean U at risk 75 white") J1=("mean U at risk 80 white") K1=("mean U at risk 85 white") ///
L1=("mean U at risk 90 white") ///
M1=("mean U at risk birth black") N1=("mean U at risk 45 black") O1=("mean U at risk 50 black") P1=("mean U at risk 55 black") ///
Q1=("mean U at risk 60 black") R1=("mean U at risk 65 black") S1=("mean U at risk 70 black") ///
T1=("mean U at risk 75 black") U1=("mean U at risk 80 black") V1=("mean U at risk 85 black") ///
W1=("mean U at risk 90 black") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_meanUatrisk0_exp0)  C$outputrow=(mean_meanUatrisk45_exp0) D$outputrow=(mean_meanUatrisk50_exp0) E$outputrow=(mean_meanUatrisk55_exp0) ///
F$outputrow=(mean_meanUatrisk60_exp0) G$outputrow=(mean_meanUatrisk65_exp0) H$outputrow=(mean_meanUatrisk70_exp0) ///
I$outputrow=(mean_meanUatrisk75_exp0) J$outputrow=(mean_meanUatrisk80_exp0) K$outputrow=(mean_meanUatrisk85_exp0) ///
L$outputrow=(mean_meanUatrisk90_exp0) ///
M$outputrow=(mean_meanUatrisk0_exp1) N$outputrow=(mean_meanUatrisk45_exp1) O$outputrow=(mean_meanUatrisk50_exp1) P$outputrow=(mean_meanUatrisk55_exp1) ///
Q$outputrow=(mean_meanUatrisk60_exp1) R$outputrow=(mean_meanUatrisk65_exp1) S$outputrow=(mean_meanUatrisk70_exp1) ///
T$outputrow=(mean_meanUatrisk75_exp1) U$outputrow=(mean_meanUatrisk80_exp1) V$outputrow=(mean_meanUatrisk85_exp1) ///
W$outputrow=(mean_meanUatrisk90_exp1) ///
using SimulationResults_N`N_rounded'_B`B', sheet("meanUatrisk") modify

/***sheet 7: N at risk at each age***/
putexcel A1=("Scenario") ///
B1=("N at risk 45 white") C1=("N at risk 50 white") D1=("N at risk 55 white") ///
E1=("N at risk 60 white") F1=("N at risk 65 white") G1=("N at risk 70 white") ///
H1=("N at risk 75 white") I1=("N at risk 80 white") J1=("N at risk 85 white") ///
K1=("N at risk 90 white") ///
L1=("N at risk 45 black") M1=("N at risk 50 black") N1=("N at risk 55 black") ///
O1=("N at risk 60 black") P1=("N at risk 65 black") Q1=("N at risk 70 black") ///
R1=("N at risk 75 black") S1=("N at risk 80 black") T1=("N at risk 85 black") ///
U1=("N at risk 90 black") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_Natrisk45_exp0) C$outputrow=(mean_Natrisk50_exp0) D$outputrow=(mean_Natrisk55_exp0) ///
E$outputrow=(mean_Natrisk60_exp0) F$outputrow=(mean_Natrisk65_exp0) G$outputrow=(mean_Natrisk70_exp0) ///
H$outputrow=(mean_Natrisk75_exp0) I$outputrow=(mean_Natrisk80_exp0) J$outputrow=(mean_Natrisk85_exp0) ///
K$outputrow=(mean_Natrisk90_exp0) ///
L$outputrow=(mean_Natrisk45_exp1) M$outputrow=(mean_Natrisk50_exp1) N$outputrow=(mean_Natrisk55_exp1) ///
O$outputrow=(mean_Natrisk60_exp1) P$outputrow=(mean_Natrisk65_exp1) Q$outputrow=(mean_Natrisk70_exp1) ///
R$outputrow=(mean_Natrisk75_exp1) S$outputrow=(mean_Natrisk80_exp1) T$outputrow=(mean_Natrisk85_exp1) ///
U$outputrow=(mean_Natrisk90_exp1) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Natrisk") modify

/***sheet 8: N strokes at risk at each age***/
putexcel A1=("Scenario") ///
B1=("N strokes 45-50 white") C1=("N strokes 50-55 white") D1=("N strokes 55-60 white") ///
E1=("N strokes 60-65 white") F1=("N strokes 65-70 white") G1=("N strokes 70-75 white") ///
H1=("N strokes 75-80 white") I1=("N strokes 80-85 white") J1=("N strokes 85-90 white") ///
K1=("N strokes 90-95 white") ///
L1=("N strokes 45-50 black") M1=("N strokes 50-55 black") N1=("N strokes 55-60 black") ///
O1=("N strokes 60-65 black") P1=("N strokes 65-70 black") Q1=("N strokes 70-75 black") ///
R1=("N strokes 75-80 black") S1=("N strokes 80-85 black") T1=("N strokes 85-90 black") ///
U1=("N strokes 90-95 black") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_nstrokes45to50_exp0) C$outputrow=(mean_nstrokes50to55_exp0) D$outputrow=(mean_nstrokes55to60_exp0) ///
E$outputrow=(mean_nstrokes60to65_exp0) F$outputrow=(mean_nstrokes65to70_exp0) G$outputrow=(mean_nstrokes70to75_exp0) ///
H$outputrow=(mean_nstrokes75to80_exp0) I$outputrow=(mean_nstrokes80to85_exp0) J$outputrow=(mean_nstrokes85to90_exp0) ///
K$outputrow=(mean_nstrokes90to95_exp0) ///
L$outputrow=(mean_nstrokes45to50_exp1) M$outputrow=(mean_nstrokes50to55_exp1) N$outputrow=(mean_nstrokes55to60_exp1) ///
O$outputrow=(mean_nstrokes60to65_exp1) P$outputrow=(mean_nstrokes65to70_exp1) Q$outputrow=(mean_nstrokes70to75_exp1) ///
R$outputrow=(mean_nstrokes75to80_exp1) S$outputrow=(mean_nstrokes80to85_exp1) T$outputrow=(mean_nstrokes85to90_exp1) ///
U$outputrow=(mean_nstrokes90to95_exp1) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Nstrokes") modify

/***sheet 9: P strokes at risk at each age***/
putexcel A1=("Scenario") ///
B1=("P strokes 45-50 white") C1=("P strokes 50-55 white") D1=("P strokes 55-60 white") ///
E1=("P strokes 60-65 white") F1=("P strokes 65-70 white") G1=("P strokes 70-75 white") ///
H1=("P strokes 75-80 white") I1=("P strokes 80-85 white") J1=("P strokes 85-90 white") ///
K1=("P strokes 90-95 white") ///
L1=("P strokes 45-50 black") M1=("P strokes 50-55 black") N1=("P strokes 55-60 black") ///
O1=("P strokes 60-65 black") P1=("P strokes 65-70 black") Q1=("P strokes 70-75 black") ///
R1=("P strokes 75-80 black") S1=("P strokes 80-85 black") T1=("P strokes 85-90 black") ///
U1=("P strokes 90-95 black") ///
A$outputrow=("$causalscenario") ///
B$outputrow=(mean_p_stroke45to50_exp0) C$outputrow=(mean_p_stroke50to55_exp0) D$outputrow=(mean_p_stroke55to60_exp0) ///
E$outputrow=(mean_p_stroke60to65_exp0) F$outputrow=(mean_p_stroke65to70_exp0) G$outputrow=(mean_p_stroke70to75_exp0) ///
H$outputrow=(mean_p_stroke75to80_exp0) I$outputrow=(mean_p_stroke80to85_exp0) J$outputrow=(mean_p_stroke85to90_exp0) ///
K$outputrow=(mean_p_stroke90to95_exp0) ///
L$outputrow=(mean_p_stroke45to50_exp1) M$outputrow=(mean_p_stroke50to55_exp1) N$outputrow=(mean_p_stroke55to60_exp1) ///
O$outputrow=(mean_p_stroke60to65_exp1) P$outputrow=(mean_p_stroke65to70_exp1) Q$outputrow=(mean_p_stroke70to75_exp1) ///
R$outputrow=(mean_p_stroke75to80_exp1) S$outputrow=(mean_p_stroke80to85_exp1) T$outputrow=(mean_p_stroke85to90_exp1) ///
U$outputrow=(mean_p_stroke90to95_exp1) ///
using SimulationResults_N`N_rounded'_B`B', sheet("Pstroke") modify


/**************************************************************************************************/
/***	export data to Excel																	***/
/***	one row = one simulated sample															***/
/**************************************************************************************************/
export excel using EachSimulationResults_N`N_rounded'_B`B', sheet("$causalscenario") sheetmodify firstrow(variables)

/**************************************************************************************************/
/***	save data in Stata																		***/
/***	one row = one simulated sample															***/
/**************************************************************************************************/
save "each_sim_results_$causalscenario.dta", replace


/**************************************************************************************************/
/***	histograms of mean U at risk-these figures are not included in paper-Figure 4			***/
/***	plots the average mean U by age and race for each scenario; these histograms			***/
/***	show the distribution of mean U by race for each age group and scenario  				***/
/***	(one plot for each age group and each scenario).										***/
/**************************************************************************************************/
twoway (histogram meanUatrisk0_exp0, fcolor(none) lcolor(blue)) ///
       (histogram meanUatrisk0_exp1, fcolor(none) lcolor(red)), ///
	   xlabel(-1.25(0.25)0.25) legend(order(1 "white" 2 "black")) graphregion(fcolor(white)) title("$causalscenario: mean U at risk, birth", color(black) size(med))
	   graph save Graph "meanUhistogram0_$causalscenario.gph", replace

twoway (histogram meanUatrisk45_exp0, fcolor(none) lcolor(blue)) ///
       (histogram meanUatrisk45_exp1, fcolor(none) lcolor(red)), ///
	   xlabel(-1.25(0.25)0.25) legend(order(1 "white" 2 "black")) graphregion(fcolor(white)) title("$causalscenario: mean U at risk, age 45", color(black) size(med))
	   graph save Graph "meanUhistogram45_$causalscenario.gph", replace
    
twoway (histogram meanUatrisk55_exp0, fcolor(none) lcolor(blue)) ///
       (histogram meanUatrisk55_exp1, fcolor(none) lcolor(red)), ///
	   xlabel(-1.25(0.25)0.25) legend(order(1 "white" 2 "black")) graphregion(fcolor(white)) title("$causalscenario: mean U at risk, age 55", color(black) size(med))
	   graph save Graph "meanUhistogram55_$causalscenario.gph", replace	   
	   
twoway (histogram meanUatrisk65_exp0, fcolor(none) lcolor(blue)) ///
       (histogram meanUatrisk65_exp1, fcolor(none) lcolor(red)), ///
	   xlabel(-1.25(0.25)0.25) legend(order(1 "white" 2 "black")) graphregion(fcolor(white)) title("$causalscenario: mean U at risk, age 65", color(black) size(med))
	   graph save Graph "meanUhistogram65_$causalscenario.gph", replace
	   
twoway (histogram meanUatrisk75_exp0, fcolor(none) lcolor(blue)) ///
       (histogram meanUatrisk75_exp1, fcolor(none) lcolor(red)), ///
	   xlabel(-1.25(0.25)0.25) legend(order(1 "white" 2 "black")) graphregion(fcolor(white)) title("$causalscenario: mean U at risk, age 75", color(black) size(med))
	   graph save Graph "meanUhistogram75_$causalscenario.gph", replace
	   
twoway (histogram meanUatrisk85_exp0, fcolor(none) lcolor(blue)) ///
       (histogram meanUatrisk85_exp1, fcolor(none) lcolor(red)), ///
	   xlabel(-1.25(0.25)0.25) legend(order(1 "white" 2 "black")) graphregion(fcolor(white)) title("$causalscenario: mean U at risk, age 85", color(black) size(med))
	   graph save Graph "meanUhistogram85_$causalscenario.gph", replace
 
 
di "$S_TIME"

 timer off 1
 
 timer list 1

