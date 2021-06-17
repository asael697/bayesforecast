**bayesforecast NEWS**
============

**Bayesforecast 1.0.1. Date: 16/06/2021**
----------------------------------

### Features:

- Local level models using local_level() and stan_local_level() methods

- Adding introductory vignette.

- Add mcmc_plot for posterior check

### Improvements:

-  adding seasonal logical function to auto.sarima() function for avoiding Seasonal models.


### Changes:

-   No current changes

### Fixes:

- Correcting the default number of iterations per chain to 2000.

- Correcting ets Stan code, for seasonality. 

- correcting bug in predictive_error() function.

-fix doi:10.2307/1392556 in GARCH vignette

-fix URl in GARCH vignette

**Bayesforecast 0.1.1. Date: 29/01/2021**
----------------------------------

### Features:

- bayesforecast package posted on CRAN


### Improvements:

-  No current changes


### Changes:

-   No current changes

### Fixes:

- Correcting the default number of iterations per chain to 2000.

- Correcting ets stan code, for seasonality. 


**Bayesforecast 0.0.1. Date: 07/12/2020**
----------------------------------

### Features:

- bayesforecast package created.

- Less stan models to compile than the previous *varstan* package.

- simplified interface.


### Improvements:

-  *stan_model()* functions created for similar syntax to *brms* and *rstanarm* packages.

- *forecast* EST methods implemented for Bayesian estimation.

- local global trend methods implemented from the *Orbit* package.

- automatic forecast methods for similar syntax to *forecast* package.


### Changes:

-   No current changes

### Fixes:

-   No current changes
