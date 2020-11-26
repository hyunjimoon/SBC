# SBC

## Efficient simulation-based calibration for Bayesian models

SBC is an R package that provide tools to easily validate and offer correction on prior, likelihood, and computation algorithm. Tools include SBC rank plot, ECDF, and ECDF\_diff plot and their summary statistics. Modelers could, and recommended to, upgrade their model based on multiple test results.

---
### Currently supports:

* Rank Histogram
* ECDF plot
* Centered Histogram plot
* Automated sampling of P(\theta),P(y | \theta) for cmdstan and rstan models

---
### TODO:
* Uniformity checks
* More envelope metrics
* Add verbose diagnostics(somewhat like get\_hmc\_diagnostics)
* Refine SBCData for better support
* **Automate prior and posterior sampling** (WIP)
* Inferential Calibration
* ECDF_diff plot

---
### References:
[Validating Bayesian Inference
Algorithms with Simulation-Based
Calibration](https://arxiv.org/pdf/1804.06788.pdf), Sean Talts, Michael Betancourt, et al.

[Rank-Normalization, Folding, and Localization: An Improved $\hat{R}$ for Assessing Convergence of MCMC](https://arxiv.org/abs/1903.08008), Aki Vehtari, Andrew Gelman, et al.

