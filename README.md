pq-sim
======

R and ADMB code for dealing with pesky zeros

Instructions:
1. source(pqBatch.r)
2. Set fixed "data" values n, F and Rbar on the GUI
3. Set "True" parameter values M ,ah, gh, qtil, b, N on the GUI. These are the parameters to be estimated in ADMB
M: Natural mortality
ah: Age at 50% harvest in the logistic selectivity function
gh: standard deviation of the logistic selectivity function
qtil: parameter in the relationship between pi and qi (0 < qtil < 1)
b:  parameter in the relationship between pi and qi (b > 0) 
N: Effective sample size for Dirichlet distribution

4. Estimation options:
i. Simulate nsamples datasets and run MPD estimation nsamples times
ii. Simulate one dataset and run an MCMC estimation: you can set the number of MCMC samples on the GUI

CHANGED in this version
April 29 2014: 
1. Reduced upper bound on b in PARAMETER_SECTION of pq.tpl. As pi tends to zero, logit(pi) tends to a large negative number, and (Eq T2.2) tends to a large positive number	and qi tends to 1. When Eq T2.2 gets to be around 44, qi = 1 and the likelihood function crashes due to log(1-qi)
2. Corrected objective function so that the second term is always calculated.
