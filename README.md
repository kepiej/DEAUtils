# DEAUtils
Utility functions for Data Envelopment Analysis (DEA)

This repository will host various utility functions that one may need when doing a DEA analysis.

## Nonparametric kernel density test for testing equality of densities
Currently, it contains a MATLAB implementation of the nonparametric kernel density test of Li et al. (2009) in the function "litest2009". This function can be used as follows:

<code>
[h,p,Tn] = litest2009(X,Y,alpha,nboot)
</code>

where <code>X</code> and <code>Y</code> are i.i.d samples of some dimension p, but possibly with different sample sizes. <code>alpha</code> is the desired significance level (default = 0.05) and <code>nboot</code> the number of bootstrap samples (default = 2000). The null hypothesis states that both densities of samples <code>X</code> and <code>Y</code> are equal.
The function returns 3 arguments:
- <code>h</code> the result of the hypothesis test: h = 0 means the null hypothesis cannot be rejected at the given significance level <code>alpha</code>. h = 1 means the null hypothesis is rejected at the given significance level.
- <code>p</code> the p-value obtained by a bootstrapping procedure using <code>nboot</code> bootstrap iterations. See Li et al.(2009) for details.
- <code>Tn</code> the computed test statistic.

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
The <code>litest2009</code> implementation can only handle continuous data at the moment. If there is any interest then I might implement support for categorical data in the future.


## References
* Qi Li, Esfandiar Maasoumi & Jeffrey S. Racine (2009) "A nonparametric test for equality of distributions with mixed categorical and continuous data", Journal of Econometrics, 148, 186-200, DOI: http://dx.doi.org/10.1016/j.jeconom.2008.10.007
*  Léopold Simar & Valentin Zelenyuk (2006) "On Testing Equality of Distributions of Technical Efficiency Scores", Econometric Reviews, 25:4, 497-522, DOI:10.1080/07474930600972582
