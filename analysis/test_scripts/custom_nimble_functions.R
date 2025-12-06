dBetaBinom_One <- nimbleFunction(
  run = function(x = double(0),
                 N = double(0),
                 shape1 = double(0),
                 shape2 = double(0),
                 log = integer(0, default = 0)) {
    
    logprob <- 0
    logprob <- logprob +
      nimBetaFun(a = x + shape1, b = N - x + shape2, log = TRUE) -
      nimBetaFun(a = shape1, b = shape2, log = TRUE) +
      lgamma(N+1) - (lgamma(x+1) + lgamma(N - x + 1))
    if (log) return(logprob)
    return(exp(logprob))
    returnType(double(0))
  }
)

rBetaBinom_One <- nimbleFunction(
  run = function(n = double(0),
                 N = double(0),
                 shape1 = double(0),
                 shape2 = double(0)) {
    p <- rbeta(1, shape1, shape2)
    x <- rbinom(1, N, p)
    return(x)
    returnType(double())
  })