# workflow test code

generator <- function(){
  function(J){
    mu <- rnorm(1, 0, 5)
    tau <- rcauchy(1, 0, 5)
    theta_trans <- rnorm(J, 0, 1)
    theta <-theta_trans * tau + mu
    list(
      generated = rnorm(J, mu, sigma),
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
ncp_model = cmdstanr::cmdstan_model("tests/eightschools_ncp.stan")


workflow <- SBC::SBCWorkflow$new(ncp_model, generator())

workflow$simulate(10, J)
d <- workflow$fit_model(4000, 4000, data)
ranks <- workflow$calculate_rank()
