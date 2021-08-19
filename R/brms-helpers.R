stanmodel_for_brms <- function(...) {
  model_code <- brms::make_stancode(...)

  args <- list(...)
  if(!is.null(args$sample_prior)) {
    stop("Do not specify `sample_prior`")
  }
  if(!is.null(args$empty)) {
    stop("Do not specify `empty`")
  }
  backend <- args$backend
  if(is.null(backend)) {
    backend <- getOption("brms.backend", "rstan")
  }
  if(backend == "cmdstanr") {
    compiled_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(model_code))
  } else if(backend == "rstan") {
    compiled_model <- rstan::stan_model(model_code = model_code)
  } else {
    stop("Unsupported backend: ", backend)

  }

  compiled_model
}

translate_rstan_args_to_cmdstan <- function(args, include_unrecognized = TRUE) {
  ignored_args <- c("cores", "data")
  recognized_but_unchanged <- c("thin", "refresh")
  translated_args <- list()
  for(old in names(args)) {
    if(old == "control") {
      if(!is.null(args$control$adapt_delta)) {
        translated_args$adapt_delta = args$control$adapt_delta
      }
      if(!is.null(args$control$max_treedepth)) {
        translated_args$max_depth = args$control$max_treedepth
      }
    } else if(old == "iter") {
      if(is.null(args$warmup)) {
        translated_args$iter_warmup = args$iter / 2
        translated_args$iter_sampling = args$iter/ 2
      } else {
        translated_args$iter_warmup = args$warmup
        translated_args$iter_sampling = args$iter - args$warmup
      }
    } else if(old == "warmup") {
      if(is.null(args$iter)) {
        translated_args$iter_warmup = args$warmup
      } #If iter is present, warmup will be handled when handling iter
    } else if(!(old %in% ignored_args)) {
      if(include_unrecognized || old %in% recognized_but_unchanged) {
        translated_args[[old]] = args[[old]]
      }
    }
  }
  translated_args
}

sampling_backend_from_stanmodel <- function(stanmodel, args) {

  if(inherits(stanmodel, "CmdStanModel")) {
    translated_args <- translate_rstan_args_to_cmdstan(args)

    do.call(SBC_backend_cmdstan_sample, combine_args(translated_args, list(model = stanmodel)))
  } else if(inherits(stanmodel, "stanmodel")) {
    do.call(SBC_backend_rstan_sample, combine_args(args,list(model = stanmodel)))
  } else {
    stop("stanmodel does not inherit from `stanmodel` or `CmdStanModel`")
  }
}

brmsfit_from_stanfit <- function(fit, brmsargs) {
  fit_brms <- do.call(brms::brm, combine_args(brmsargs, list(empty = TRUE)))
  if(inherits(fit, "CmdStanMCMC")) {
    fit_brms$fit <- rstan::read_stan_csv(fit$output_files())
  } else {
    fit_brms$fit <- fit
  }
  fit_brms <- brms::rename_pars(fit_brms)

  fit_brms
}
