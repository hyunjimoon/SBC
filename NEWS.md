# SBC 0.5.0

## Breaking changes

* We changed how backend diagnostics are reported, which may affect your custom backends, see https://hyunjimoon.github.io/SBC/articles/implementing_backends.html
  for example on how to implement custom diagnostics for the new standard. See help `diagnostic_types()` for some more details.
* Some results cache files may be invalidated and may need to be recomputed to work properly (this sould generally happen automatically)

## What's changed

* Support for Bayes Factor workflow via `bridgesampling`, `BayesFactor` - see https://hyunjimoon.github.io/SBC/articles/bayes_factor.html 
  the code is used to support our preprint on SBC for Bayes factors https://arxiv.org/abs/2508.11814 
  and some further example usage can be seen in the accompanying repo https://github.com/martinmodrak/SBCBayesFactors/
  This includes methods to check calibration for binary variables that are more efficient than the eCDF plots.
* `SBC_datasets` now support setting "variable attributes" (`var_attributes()`) that provide further info for correct handling
  of variables (e.g. allow NA, treat variable as binary)
* Backends may implement `SBC_posterior_cdf()` to provide exact continuous rank with respect to CDFs instead of ranks in samples.
  This is currently used only for posterior model probabilities in Bayes factors but has some potential for wider use with e.g. variational
  inference or Laplace approximations where it would let us bypass samples completely and gain a bit more precision.
* `SBC_backend_cached()` lets you wrap any backend and cache the results of the specific fits (instead of caching the whole results)
* `SBC_generator_function()` can now be paralellized
* Support `fullrank` algorithm for `brms` backends.

**Full Changelog**: https://github.com/hyunjimoon/SBC/compare/v0.3.0...v0.5.0

# SBC 0.3.0

## What's Changed

* Add support for multivariate brms models by @bmfazio in https://github.com/hyunjimoon/SBC/pull/81
* Skip data validation inside `brms::posterior_predict` to speed up generation by @bmfazio in https://github.com/hyunjimoon/SBC/pull/82
* Create CITATION by @hyunjimoon in https://github.com/hyunjimoon/SBC/pull/92
* Updating stan code to the most recent syntax by @timostenz in https://github.com/hyunjimoon/SBC/pull/94
* Add `combine` argument to ecdf plots by @bmfazio in https://github.com/hyunjimoon/SBC/pull/93
* separate ecdf_alpha logic from combine by @bmfazio in https://github.com/hyunjimoon/SBC/pull/96
* out_stan_file for brms backends to avoid recompilation
* init_factory as a way to support data-dependent inits (currently for cmdstan sampling only)
* Other minor fixes and quality of life improvements

## New Contributors
* @bmfazio made their first contribution in https://github.com/hyunjimoon/SBC/pull/81
* @timostenz made their first contribution in https://github.com/hyunjimoon/SBC/pull/94

**Full Changelog**: https://github.com/hyunjimoon/SBC/compare/v0.2.0...v0.3.0

# SBC 0.2.0

## What's changed

- Minor improvements in performance
- Renaming and deprecating a bunch of functions and improving documentation to be generally more sensible and better match the [preprint describing our main ideas](https://arxiv.org/abs/2211.02383v1) - all older code should still work, but may issue some deprecation warnings
