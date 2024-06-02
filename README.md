# Mortality

Last Updated: May 2023
Replication Instructions for Lleras-Muney
and Moreau (2022) Demography


1-Data
The data for this project comes from the Human Mortality Database and is publicly available but
requires registration. The data can be downloaded from https://www.mortality.org/


2-Codes
There are several codes to produce the results of the paper.

a. STATA codes

1-hmd_read.do: this file reads in data downloaded from the Human Mortality Database and
outputs a data set called hmd_drate.dta
the input files are labeled Mx_1x1 in the website. The program assumes these have been read into
STATA and the first few lines of text have been removed.
For simplicity these figures use the mortality rates that the human mortality database provides, that
is we do not re-compute these rates for these figures.

2-Figures_AERI.do
*uses data from the Human Mortality Database (saved as .dta by hmd.do) to produce Figure 1 and
Appendix Figure 11
*simulates data to produce Figures 4b, 4d, and 5d and Appendix Figures 7, 8, 9, 14, 15, 16 and 17.
The seed and parameters used in the simulations are in the code.

3. GeneratetargetMRv4.do
This code generare mortality rate series in correct format that can be used for matlab
estimation. Outputs one csv file per cohort, with format ‘MR_FRA_YYYY_Fq.csv’where YYYY is the
year

b. Matlab codes:

1- main_figuresv5_hoff_V21.m produces all the figures based on estimated series

2- hof_timeseriesv6_L.m is a wrapper that runs the estimation for all series using the
Hoffman server and saves the data in a table, launching the estimation recursively

TS_fun_kII.m function that estimates for a given cohort based on a guess it
uses the following ancillary functions. From the year the function will automatically pick up the
corresponding data file saved locally.

Criterion_s.m wrapper, forms the squared errors with the empirical survival
curve to minimize

Crit_fun.m wrapper, forms the squared errors with the empirical mortality curve to
minimize

Trunc_data.m truncates the mortality rate data to avoid the last years that are
extrapolated (see Human Mortality Database Documentation)

Survivalcurve.m compute the survival curve from the log mortality rates

Powell.m a routine to improve on Matlab’s default maximization
algorithms

Life_exp4.m compute the life expectancy

logMR_IIk.m mortality model with 2 investment shocks and one accident
shock through lifteime

logMR_add.m used for the monkeys

logMR.m is the most basic version of the 5-parameter model.

Example: To estimate a single series
year = 1816
init_guess = [0.3665 0.00047854 0.83838 1.7934 0.91652 0.3665 0.0093653 0.3665]
[res,fit_S,LMR1,S1,LE1,fit0]=TSfun_kII(year,init_guess)
% figure(),
% subplot(1,2,1)
% plot(LMR1)
% subplot(1,2,2)
% plot(S1)

c. Excel code: Figure 13 is simulated in excel using the estimated parameters.
