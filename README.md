# Simulation-based Calibration: SBC
## Efficient simulation-based calibration for Bayesian models
SBC provides tools to easily validate and offer corrections on prior, likelihood, and computation algorithms based on the self-recovering property of Bayesian models. This package contains tools such as SBC rank histograms, ECDF plots, and their summary statistics which can be used to assess computational faithfulness. 

### Varieties of calibrations: scope of this package
Calibration (i.e. reliability) is not a sufficient condition for a good forecast but a minimal property that any forecast should satisfy (FV1998). It serves as a bootstrap for model development and its method and target varies. Target is chosen as modeler's quantity of interest and directly affects the calibrated result as reward in reinforcement learning. Method depends on how much you marginalized or conditioned the full joint space to test coverage. Scope of this package is checked below.

## Interface and Usage

SBC is designed to be primarily used with [Stan](https://mc-stan.org/) models, offering a highly customizable interface to integrate Simulation Based Calibration into existing Bayesian workflows with minimal effort. Its main feature is the `api` interface, which defines a fully-blown SBC pipeline starting from dataset generation to posterior sampling. Once a user has a valid Stan model and a minor R function defining the data generating process(referred to as `Generator`), running SBC becomes as simple as:

```
n_datasets <- 100  # Number of SBC iterations to run

sbc_generator <- SBC::function_SBC_generator(Generator)
sbc_dataset <- SBC::generate_datasets(
  sbc_generator, 
  n_datasets)

cmdstan_backend <- SBC::cmdstan_sample_SBC_backend(
    cmdstan_model, iter_warmup = 1000, iter_sampling = 1000)
    
results <- SBC::compute_results(sbc_dataset, cmdstan_backend)
plot_rank_hist(results)
```

For detailed usage, please refer to the included vignettes.

### Compatibility
Currently `SBC` supports `cmdstan`, `rstan`, and `brms` models out of the box. However, adding Backends for other platforms is supported.

## Installation
To install the development version of SBC, run
```
devtools::install_github("hyunjimoon/SBC")
```
from your R console.

### References:
Theoretical support
* [Validating Bayesian Inference Algorithms with Simulation-Based Calibration](https://arxiv.org/pdf/1804.06788.pdf) Talts, Betancourt, Simpson, Vehtari, Gelman, 2018
* [Graphical Test for Discrete Uniformity and its Applications in Goodness of Fit Evaluation and Multiple Sample Comparison](https://arxiv.org/abs/2103.10522)  Säilynoja, Bürkner, Vehtari, 2021
* [Bayesian Workflow](https://arxiv.org/abs/2011.01808), Gelman et al., 2020
* [Toward a principled Bayesian workflow in cognitive science](https://psycnet.apa.org/record/2020-43606-001) Schad, Betancourt, Vasishth, 2021
* [Bayes factor workflow](https://arxiv.org/pdf/2103.08744.pdf) Schad, Nicenboim, Bürkner, Betancourt, Vasishth, 2021

Application support
* [Cognitive science, response time fitting](https://link.springer.com/content/pdf/10.3758/s13428-019-01318-x.pdf)
* [Bioinformatics, effect of mutation prediction](https://www.biorxiv.org/content/10.1101/2020.10.27.356758v1.full.pdf)
* [Earth science, earthquake prediction](https://gmd.copernicus.org/articles/11/4383/2018/gmd-11-4383-2018.pdf )
* [Sequential Neural Likelihood](http://proceedings.mlr.press/v89/papamakarios19a/papamakarios19a.pdf) 

Vignette
* [ECDF with codes](https://avehtari.github.io/rhat_ess/rhat_ess.html) (new implementation by Teemu Säilynoja will be available in `bayesplot` and `SBC` package soon)

## FAQ
> How does calibration relate to prediction accuracy?

Comparing the ground truth and the simulated result is a backbone of calibration and comparison target greatly affects the calibrated (i.e. trained) result, similar to reward in reinforcement learning. In this sense, if the U(a(y), theta) term is designed for prediction, the model will be calibrated to have best predictive result as possible.
