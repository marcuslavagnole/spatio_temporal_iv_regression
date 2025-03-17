In this repository, one can find the R routines used in the article
[Spatio-temporal instrumental variables regression with missing data: A Bayesian approach](https://doi.org/10.1007/s10614-022-10269-z). It is a joint work with [Kelly C. M. Gonçalves](https://sites.google.com/dme.ufrj.br/kelly/) and Mario Jorge Mendonça, published in Computational Economics. 

The paper proposes an extension of the Bayesian instrumental variables regression, which allows spatial and temporal correlation among observations. For that, we introduce a double separable covariance matrix, adopting a Conditional Autoregressive structure for the spatial component and a first-order autoregressive process for the temporal component. We also introduce a Bayesian multiple imputation to handle missing data considering uncertainty. The methodology is applied for the estimation of how broadband affects the Gross Domestic Product of municipalities in the state of Mato Grosso do Sul from 2010 to 2017.

This repo includes:

- **MCMC_Code_IVST.R** : Main file containing the MCMC routine; 
- **Full_Conditionals_IVST.R** : Auxiliary file with all full conditional distributions;
- **Missing_Imputation.R** : File with the multiple imputation function;
- **aux_func.R** : File with the auxiliary functions;
- **data.RData** : real data set;
- **W_MS.RData** : adjacency matrix for the municipalities of Mato Grosso do Sul state.
