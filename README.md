# Multilevel-Poisson-simulation
Here is the study design:
Background:
The estimation of power in two-level models used to analyze data that are hierarchically structured is particularly complex because the outcome contains variance at two levels that is regressed on predictors at two levels. Methods for the estimation of power in two-level models have been based on formulas and Monte Carlo simulation. However, there is scant literature or Monte Carlo simulations focused on ordinal/count data outcomes estimation. 

Study design: 
Test of statistical power for a two-level model with the proportional odds model
Yij ~ (poisson (lamda))
ln(lamda) = b0 + b1(X)
b0 = g00 + g01(W) + u0j
b1 = g10 + g11(W)
Where you ~ N(0, t00)
Constants: g00, g01, g10,
Nj (number of cases per cluster) = 20
Tau matrix (variance of the slope) = half size of the intercept, covariance is zero
Variance of u0j = t00 = 0.5
Y is count data with poisson distribution
X & W are continuous, normally distributed and homoscedastic (X & W, standardized normal variables)
Define the values of 
g00 = 2 
g01 = .5 - a medium effects
g10 = .5 
Number of MC replications per condition = 1000
Variance of u0j = t00 = (0.004, 0.0015, 0.001)
IV1 (3-levels): ratio of  ICC =  (0.25, 0.09, 0.06)
IV2 (2-levels): g11 = (.2, 0.45)
IV3 (3-levels): #-clusters = 30, 50, 100
DV for the simulation
Power ( Type-I error rate)
Bias/relative bias
g00, g01, g10, g11 
