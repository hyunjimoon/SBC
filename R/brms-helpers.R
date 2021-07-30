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
  } else {
    stop("Only cmdstanr backend currently supported")
  }

  compiled_model
}

sampling_backend_from_stanmodel <- function(stanmodel, args) {

  if(inherits(stanmodel, "CmdStanModel")) {
    ignored_args <- c("cores", "data")
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
        translated_args[[old]] = args[[old]]
      }
    }

    do.call(cmdstan_sample_SBC_backend, c(list(model = stanmodel), translated_args))
  } else if(inherits(stanomodel, "stanmodel")) {
    stop("rstan backend not supported yet")
  }
}

brmsfit_from_stanfit <- function(fit, brmsargs) {
  fit_brms <- do.call(brms::brm, c(list(empty = TRUE), brmsargs))
  if(inherits(fit, "CmdStanMCMC")) {
    fit_brms$fit <- rstan::read_stan_csv(fit$output_files())
  } else {
    fit_brms$fit <- fit
  }
  fit_brms <- brms::rename_pars(fit_brms)

  fit_brms
}
