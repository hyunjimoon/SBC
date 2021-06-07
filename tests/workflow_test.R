# workflow test code

generator <- function(){
  function(N){
    mu = 5
    sigma = 1
    list(
      generated = rnorm(N, mu, sigma),
      parameters = list(
        mu = 5,
        sigma=1,
        tau <- c(1,2,3,4,5)
      )
    )
  }
}

test_stan_model <- 1
workflow <- SBC::SBCWorkflow$new(test_stan_model, generator())

prior <- workflow$simulate(100, 100)
