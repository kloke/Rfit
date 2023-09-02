Rfit: Rank-based estimation for linear models 
---------------------------------------------

### Installation ###
* Rfit may be installed from 
	* [CRAN](https://cran.r-project.org/) e.g. `install.packages('Rfit')`
	* [github](https://github.com/) e.g. `install_github('kloke/Rfit')` 
		* requires [devtools](https://cran.r-project.org/package=devtools)
* Rfit no longer requires the package [quantreg](https://cran.r-project.org/package=quantreg)
  The inital fit is now based on least squares to avoid additional dependencies.  One may still use quantreg to obtain the inital fit which will ensure the robustnesses of the result.

### Background ###
Rank-based (R) estimation for statistical models is a robust nonparametric alternative to classical estimation procedures such as least squares. R methods have been developed for models ranging from linear models, to linear mixed models, to timeseries, to nonlinear models. Advantages of these R methods over traditional methods such as maximum-likelihood or least squares is that they require fewer assumptions, are robust to gross outliers, and are highly efficient at a wide range of distributions. Rfit uses standard linear model syntax and includes commonly used functions for inference and diagnostic procedures.  Wilcoxon scores, the default, are robust and highly efficent relative to the least squares estimator when the errors are normally distributed.

### Fine Print ###
As with R itself, Rfit is free software and comes with ABSOLUTELY NO WARRANTY.
