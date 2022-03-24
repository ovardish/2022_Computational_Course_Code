Endogenous Grid Points
2022 Yale Computational Macro Course
Instructor Oliko Vardishvili

Adopted from original code by Lukas Nord, EUI.
--------------------------------------------------------------------------------------

MAIN.m: script, sets parameters and calls functions to solve transitions

solve_aiyagari.m: function, contains GE loop over the steady state interest rate, calls solution to HH problem and distribution

solveHH_EGM.m: function, solves HH problem by endogenous grid method
	- iterates backward on Euler Equation, continuous asset choice
	- reference: Carroll (2006, Economic Letters)

getDist_continuous.m: function, computes distribution over state-space based on continuous asset policy
	- reference: Young (2010, Journal of Economic Dynamics and Control)

solve_trans.m: function, computes transitional dynamics to an MIT productivity shock
	- updating of interest rate path based on extended path approach
