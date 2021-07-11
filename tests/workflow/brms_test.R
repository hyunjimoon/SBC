library(brms)
dat2 <- data.frame(
  x = as.vector(1:20),
  y = x + rnorm(20)
)

warmup <- 150
iter <- 200
chains <- 1

prior <- c(
  set_prior("normal(0, 0.5)", class="b", coef="Intercept"),
  set_prior("normal(1, 0.1", class="b", coef="x")
)

brmsfit_example4 <- brm(
  bf(y ~ 0 + Intercept + x),
  data = dat2,
  warmup = warmup, iter = iter, chains = chains, prior=prior,
  sample_prior="only"
)
update(brmsfit_example4)
