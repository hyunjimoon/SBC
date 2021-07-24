# workflow test code
library(SBC)
generator <- function(){
  function(){
    mu <- rnorm(1, 0, 5)
    tau <- rcauchy(1, 0, 5)
    theta_trans <- rnorm(8, 0, 1)
    theta <-theta_trans * tau + mu
    list(
      generated = rnorm(8, mu, sigma),
      parameters = list(
        mu = mu,
        theta_trans=theta_trans,
        tau = tau,
        theta = theta
      )
    )
  }
}


J <- 8
y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)

data = list("J"=J, "y"=y, "sigma"=sigma)
ncp_model = cmdstanr::cmdstan_model("tests/workflow/eightschools_ncp.stan")


workflow <- SBC::SBCWorkflow$new(ncp_model, generator())

workflow$simulate(100)
d <- workflow$fit_model(200, 200, data, thin_factor = 1)
prior <- workflow$prior_samples
post <- workflow$posterior_samples
ranks <- workflow$calculate_rank()


plot_ecdf(workflow, var="mu")
plot_ecdf_diff(workflow, var="theta")
