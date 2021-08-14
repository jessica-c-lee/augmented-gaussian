# augmented-gaussian

Citation:
Lee, J. C., Mills, L., Hayes, B. K., & Livesey, E. J. (2021). Modelling generalisation gradients as augmented Gaussian functions. *Quarterly Journal of Experimental Psychology, 74(1)*, 106-121. https://doi.org/10.1177/1747021820949470 

R code to fit augmented Gaussians to individual generalization gradients in a hierarchical Bayesian model. Uses the [rstan](https://mc-stan.org/users/interfaces/rstan) and [bayestestR](https://github.com/easystats/bayestestR) packages.

The augmented Gaussian has 4 parameters that allow it to fit asymmetrical gradients:
* mean: the location of the gradient peak
* width-: the width (SD) of the left side of the gradient
* width+: the width (SD) of the right side of the gradient
* height: the height (peak) of the gradient

## How to use
* Open the .Rproj file, add your data file to the "data" folder, input parameters in the index.R file, and run the model following the example.
* The data file must be in long format, with the grouping variable labelled as "group", the stimulus dimension variable labelled as "x" and the dependent variable labelled as "y" (see demo1.csv)

## Note
* The example code fits gradients with 21 test stimuli, with the CS+ at the midpoint of the dimension (coded as 0). All functions use the fixed stimulus values -0.5:0.1:+0.5. These values are arbitrary, but the model specification is dependent on these values.
* The code is set up to model responses ranging from 0-100. The model must be re-specified to work with a different range.

## Output generated
### For each group
* summary.csv: summary statistics for the posterior samples for each parameter
* waic.csv: Widely Applicable Information Criterion computed with the [loo package](https://cran.r-project.org/web/packages/loo/index.html) (top row) and SE (bottom row)
* loo.csv: loo computed with the [loo package](https://cran.r-project.org/web/packages/loo/index.html)(top row) and SE (bottom row)
* hdis.csv 
  - HDI lim: Highest Density Interval calculated (%)
  - HDI low: Highest Density Interval lower limit
  - HDI high: Highest Density Interval upper limit
  - p(direction): proportion of the posterior that is positive or negative (whatever is most probable). Note that this is meaningless for the width and height parameters since they can only be positive
  - rope_low: lower limit of the user-defined Region of Practical Equivalence
  - rope_high: upper limit of the user-defined Region of Practical Equivalence
  - prop_rope: proportion of the posterior contained within the ROPE limits
* gradients.jpeg: plot of the empirical gradients facetted by subject
* postpreds.jpeg: plot of the empirical gradients facetted by subject with posterior predictives overlayed

### For each dataset
* density.jpeg: 5 panelled figure with the mean generalisation gradients and posterior density plots of the 4 group-level augmented Gaussian parameters

## Contact
Contact jessica.lee@unsw.edu.au. 
