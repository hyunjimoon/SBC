# Simulation-based Calibration: SBC

SBC provides tools to validate your Bayesian model and/or a sampling algorithm via the self-recovering property of Bayesian models. This package lets you run SBC easily and perform postprocessing and visualisations of the results to assess computational faithfulness. 

## Installation

To install the development version of SBC, run

```r
devtools::install_github("hyunjimoon/SBC")
```


## Quick tour

To use SBC, you need a piece of code that generates simulated data that should
match your model (a _generator_) and a statistical model + algorithm + 
algorithm parameters that can fit the model to data (a _backend_). SBC then lets you
discover when the backend and generator don't encode the same data generating process 
(up to [certain limitations](https://hyunjimoon.github.io/SBC/articles/limits_of_SBC.html)).

For a quick example, we'll use a simple generator producing normally-distributed
data (basically `y <- rnorm(N, mu, sigma)`) with a backend in Stan that mismatches
the generator by wrongly assuming Stan parametrizes the normal distribution via
precision (i.e. it has `y ~ normal(mu, 1 / sigma ^ 2)`).

```r
library(SBC)
gen <- SBC_example_generator("normal")
# interface = "cmdstanr" or "rjags" is also supported
backend_bad <- SBC_example_backend("normal_bad", interface = "rstan")
```

_Note: Using the `cmdstanr` interface, a small number of rejected steps will be reported. Those are false positives and do not threaten validity (they happen during warmup). This is a result of difficulties in parsing the output of `cmdstanr`. We are working on a resolution._

You can use `SBC_print_example_model("normal_bad")` to inspect the model used.

We generate 50 simulated datasets and perform SBC:

```r
ds <- generate_datasets(gen, n_sims = 50)
results_bad <- compute_SBC(ds, backend_bad)
```

The results then give us diagnostic plots that immediately show a problem:
the distribution of SBC ranks is not uniform as witnessed by both the rank histogram
and the difference between sample ECDF and the expected deviations from theoretical CDF.

```r
plot_rank_hist(results_bad)
plot_ecdf_diff(results_bad)
```

We can then run SBC with a backend that uses the correct parametrization 
(i.e. with `y ~ normal(mu, sigma)`):

```r
backend_sd <- SBC_example_backend("normal_sd", interface = "rstan")
results_sd <- compute_SBC(ds, backend_sd)

plot_rank_hist(results_sd)
plot_ecdf_diff(results_sd)
```

The diagnostic plots show no problems in this case. As with any other
software test, we can observe clear failures, but absence of failures does not imply
correctness. We can however make the SBC check more thorough by using a lot of
simulations and including suitable derived quantities to guard against
[known limitations of vanilla SBC](https://hyunjimoon.github.io/SBC/articles/limits_of_SBC.html).

## Paralellization

The examples above are very fast to compute, but in real use cases, 
you almost certainly want to let the computation run in parallel via the
[`future`](https://future.futureverse.org/) package.

```r
library(future)
plan(multisession)
```

## More information

The [package vignettes](https://hyunjimoon.github.io/SBC/articles/) provide 
additional context and examples. Notably:

- [The main vignette](https://hyunjimoon.github.io/SBC/articles/SBC.html)
has more theoretical background and instructions how to integrate your own simulation code and 
models with SBC.
- [Small model workflow](https://hyunjimoon.github.io/SBC/articles/small_model_workflow.html)
discusses how SBC integrates with model implementation workflow and how you can
use SBC to safely develop complex models step-by-step.

Currently `SBC` supports `cmdstanr`, `rstan`, and `brms` models out of the box. 
With a little additional work, you can integrate SBC with any exact or approximate fitting method as shown in the [Implementing backends vignette](https://hyunjimoon.github.io/SBC/articles/implementing_backends.html).

## Citing SBC and related software
Developing and maintaining open source software is an important yet often underappreciated contribution to scientific progress. Thus, whenever you are using open source software, please make sure to cite it appropriately so that developers get credit for their work.

When using SBC, please cite the following publication:

Modrák, M., Moon, A. H., Kim, S., Bürkner, P., Huurre, N., Faltejsková, K., ... & Vehtari, A. (2023). Simulation-based calibration checking for Bayesian computation: The choice of test quantities shapes sensitivity. Bayesian Analysis, 1(1), 1-28.

Further, SBC relies on several other R packages and, of course, on R itself. To find out how to cite R and its packages, use `citation` function. SBC specifically rely on `posterior` for manipulating posterior draws and `future` for parallel processing. Also, do not forget to cite the probabilistic inference tool you use as backend (e.g. Stan, JAGS, brms, ...)

## References

* Theoretical support
   * [Simulation-Based Calibration Checking for Bayesian Computation: The Choice of Test Quantities Shapes Sensitivity](https://arxiv.org/abs/2211.02383v1) Modrák, Moon, Kim, Bürkner, Huurre, Faltejsková, Gelman, Vehtari, 2022
   * [Validating Bayesian Inference Algorithms with Simulation-Based Calibration](http://www.stat.columbia.edu/~gelman/research/unpublished/sbc.pdf) Talts, Betancourt, Simpson, Vehtari, Gelman, 2018
   * [Graphical Test for Discrete Uniformity and its Applications in Goodness of Fit Evaluation and Multiple Sample Comparison](https://arxiv.org/abs/2103.10522)  Säilynoja, Bürkner, Vehtari, 2021
   * [Bayesian Workflow](https://arxiv.org/abs/2011.01808), Gelman et al., 2020
   * [Toward a principled Bayesian workflow in cognitive science](https://psycnet.apa.org/record/2020-43606-001) Schad, Betancourt, Vasishth, 2021
   * [Bayes factor workflow](https://arxiv.org/pdf/2103.08744.pdf) Schad, Nicenboim, Bürkner, Betancourt, Vasishth, 2021

* Application support
   * [Cognitive science, response time fitting](https://link.springer.com/content/pdf/10.3758/s13428-019-01318-x.pdf)
   * [Bioinformatics, effect of mutation prediction](https://www.biorxiv.org/content/10.1101/2020.10.27.356758v1.full.pdf)
   * [Earth science, earthquake prediction](https://gmd.copernicus.org/articles/11/4383/2018/gmd-11-4383-2018.pdf )
   * [Sequential Neural Likelihood](http://proceedings.mlr.press/v89/papamakarios19a/papamakarios19a.pdf) 

## FAQ

> How does calibration relate to prediction accuracy?

Comparing the ground truth and the simulated result is a backbone of calibration and comparison target greatly affects the calibrated (i.e. trained) result, similar to reward in reinforcement learning. In this sense, if the U(a(y), theta) term is designed for prediction, the model will be calibrated to have best predictive result as possible.

## Acknowledgements

Development of this package was supported by [ELIXIR CZ](https://www.elixir-czech.cz/) research infrastructure project (Ministry of Youth, Education and Sports of the Czech Republic, Grant No: LM2018131) including access to computing and storage facilities.
