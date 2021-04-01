# SBC
## Efficient simulation-based calibration for Bayesian models
SBC provides tools to easily validate and offer corrections on prior, likelihood, and computation algorithms based on the self-recovering property of Bayesian models. This package contains tools such as SBC rank histograms, ECDF plots, and their summary statistics which can be used to assess computational faithfulness. Modelers could, and are recommended to, diagnose their model based on multiple test results.

Please refer to the included vignette for detailed usage.

Some usecases include the following:
1. Test prior and likelihood on the basis of computational consistency  
    Any canonical Bayesian Model has the *self recovering property*, in which averaging over posterior distributions fitted with samples from the prior predictive distribution will always be equal to the prior distribution.
SBC uses the above principle and evaluates the combination of the prior and likelihood model under a fixed computation algorithm. Users should choose one computation algorithm in advance, such as full HMC, ADVI, Laplace approximation.
2. Test approximation algorithms  
    Approximation based Bayesian computation is very promising but one limitation is that it can be hard to diagnose its reliability. For example, full HMC benchmark is needed to measure its error. SBC which evaluates how well an algorithm samples from the posterior distribution, given a model and a prior could be an alternative tool for measuring reliability.
---
### Currently supports:
* Rank Histogram
* ECDF plot
* Uniformity checks
* ~~Centered Histogram plot~~
* Support for automated predictive sampling of Stan models in canonical form yhat ~ P(y | theta)
* Multi-parameter SBC plots
---
### TODO:
* ECDF\_diff plot
* More envelope metrics
* Add verbose diagnostics, akin to stan's get\_hmc\_diagnostics
* Inferential Calibration
---
### References:
* [Validating Bayesian Inference
Algorithms with Simulation-Based
Calibration](https://arxiv.org/pdf/1804.06788.pdf), Sean Talts, Michael Betancourt, et al.
* [Rank-Normalization, Folding, and Localization: An Improved R-hat for Assessing Convergence of MCMC](https://arxiv.org/abs/1903.08008), Aki Vehtari, Andrew Gelman, et al.
* [Graphical Test for Discrete Uniformity and its Applications in Goodness of Fit Evaluation and Multiple Sample Comparison](https://arxiv.org/abs/2103.10522) Teemu Säilynoja, Paul-Christian Bürkner, Aki Vehtari, 2021
* [Bayes factor workflow](https://arxiv.org/pdf/2103.08744.pdf) Schad, Nicenboim, Bürkner, Betancourt, Vasishth, 2021
* [ECDF with codes](https://avehtari.github.io/rhat_ess/rhat_ess.html)
