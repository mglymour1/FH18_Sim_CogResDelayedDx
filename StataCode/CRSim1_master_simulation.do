* This code calls all other code for the Cognitive Reserve simulation
* Use this code to specify which causal scenarios are used and the
* number of iterations completed (B)

set more off
clear all

di "$S_TIME"

timer clear 1

timer on 1


foreach causalscenario in 1 2 3 {
	global causalscenario = "Scenario`causalscenario'"
	
	global outputrow = `causalscenario'+1 
	
	global B = 5000 //desired number of iterations of sample generation 
	
	include CRSim_run_simulation.do
	
	}
	
di "$S_TIME"

 timer off 1
 
 timer list 1
