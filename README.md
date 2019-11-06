# DEAUtils
Utilities for Data Envelopment Analysis (DEA)

![GitHub All Releases](https://img.shields.io/github/downloads/kepiej/DEAUtils/total)

This repository will host various utilities that one may need when doing a DEA analysis.

If you use any of these functions in your own work, then please acknowledge my work whenever possible.

## Nonparametric kernel density test for testing equality of densities
Currently, it contains a MATLAB implementation of the nonparametric kernel density test of Li et al. (2009) in the function "litest2009". This function can be used as follows:

<code>
[h,p,Tn,exitflag,bw] = litest2009(X,Y,alpha,nboot,bwmethod)
</code>

where <code>X</code> and <code>Y</code> are i.i.d samples of some dimension p, but possibly with different sample sizes. <code>alpha</code> is the desired significance level (default = 0.05), <code>nboot</code> the number of bootstrap samples (default = 2000) and <code>bwmethod</code> the bandwidth selection method (default = 'lscv'). The null hypothesis states that both densities of samples <code>X</code> and <code>Y</code> are equal. Bandwidth tuning is performed using either least-squares ('lscv') or maximum likelihood ('mlcv') cross-validation.
The function returns 5 arguments:
- <code>h</code> the result of the hypothesis test: h = 0 means the null hypothesis cannot be rejected at the given significance level <code>alpha</code>. h = 1 means the null hypothesis is rejected at the given significance level.
- <code>p</code> the p-value obtained by a bootstrapping procedure using <code>nboot</code> bootstrap iterations. See Li et al.(2009) for details.
- <code>Tn</code> the computed test statistic.
- <code>exitflag</code> the exitflag from the optimization of the bandwidth <code>bw</code>.
- <code>bw</code> the bandwidth used in the kernel estimation.

The function might print this warning message during computation while still finishing successfully:

>Exiting: Maximum number of function evaluations has been exceeded
>         - increase MaxFunEvals option.
>         Current function value: XXXXX

The warning message originates from optimizing the bandwidth smoothing parameter <code>bw</code> over an unbounded interval. The warning message indicates that this search does not converge fast enough, so the optimization is restarted from the rule-of-thumb bandwidth but now limiting the search interval. The result should still be reliable though. The warning message is printed to give full transparency to the user.

Simar & Zelenyuk (2006) propose a modification of this test for use with DEA efficiency scores. The function "smootheffscorebeforelitest" implements their "Algorithm II". See Simar & Zelenyuk (2006) for details. The function needs to be called before using the above mentioned <code>litest2009</code> as follows:

<code>
  TES = smootheffscorebeforelitest(TE,ninputs,noutputs,alpha,effval)
</code>

where the arguments are:
- TE: (k x 1) vector of DEA efficiency scores;
- ninputs: number of inputs;
- noutputs: number of outputs;
- alpha: significance level of the Li test (optional, default = 5%);
- effval: value that indicates an efficient DMU (optional, default = 1.0);

The function returns the "smoothed" efficiency scores in the vector <code>TES</code>. The inefficient scores remain untouched and only the efficient scores are "smoothed". These "smoothed" efficiency scores can then be used as arguments <code>X</code> and <code>Y</code> in <code>litest2009</code>.

Note that productivity scores need not be smoothed, because these are not bounded!

### Limitations
The <code>litest2009</code> implementation can only handle continuous data at the moment. Furthermore, the test is for unconditional densities only. If there is any interest then I might implement support for categorical data and conditional densities in the future.

### Acknowledgements
I gratefully acknowledge Qianying Jin for reporting strange results in some cases when using the least-squares cross-validation method. This led me to add the maximum likelihood cross-validation as a bandwidth selection method.

## References
* Qi Li, Esfandiar Maasoumi & Jeffrey S. Racine (2009) "A nonparametric test for equality of distributions with mixed categorical and continuous data", Journal of Econometrics, 148, 186-200, DOI: http://dx.doi.org/10.1016/j.jeconom.2008.10.007
*  Léopold Simar & Valentin Zelenyuk (2006) "On Testing Equality of Distributions of Technical Efficiency Scores", Econometric Reviews, 25:4, 497-522, DOI:10.1080/07474930600972582
