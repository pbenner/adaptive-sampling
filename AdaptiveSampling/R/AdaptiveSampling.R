make.options <- function(
	n_moments=2, 
	model_posterior=1, 
	bprob = 1, utility = 1, 
	kl_component=1, 
	kl_multibin=0, 
	effective_counts = 0,
	marginal = 1,
	marginal_step = 0.01,
	marginal_range = c(0, 1),
	epsilon = 0.00001,
	threads = 1,
	stacksize = 256*1024,
	algorithm = 0,
	which = 0,
	samples = c(100, 2000)) {
  env <- environment()
  env$n_moments             <- n_moments
  env$model_posterior       <- model_posterior
  env$bprob                 <- bprob
  env$utility               <- utility
  env$kl_component          <- kl_component
  env$kl_multibin           <- kl_multibin
  env$effective_counts      <- effective_counts
  env$marginal              <- marginal
  env$marginal_step         <- marginal_step
  env$marginal_range        <- marginal_range
  env$epsilon               <- epsilon
  env$threads               <- threads
  env$stacksize             <- stacksize
  env$algorithm             <- algorithm
  env$which                 <- which
  env$samples               <- samples

  env
}

default.alpha <- function(alpha) {
  K <- dim(alpha)[1]
  L <- dim(alpha)[2]
  result <- array(dim=c(L,L,K))
  for (k in 1:K) {
    result[,,k] <- generate.alpha(alpha[k,])
  }
  result
}

generate.alpha <- function(alpha) {
  result <- outer(rep(1,length(alpha)),alpha)
  result[lower.tri(result)] <- 0
  for (i in 1:length(alpha)) {
    result[i,] <- cumsum(result[i,])
    for (j in i:length(alpha)) {
      result[i,j] <- result[i,j]/(j-i+1)
    }
  }
  result
}

default.beta <- function(n) {
  return(rep(1, n))
}

default.gamma <- function(n) {
  result <- outer(rep(1,n),rep(1,n))
  result[lower.tri(result)] <- 0
  result
}

count.statistic <- function(counts) {
  K <- dim(counts)[1]
  L <- dim(counts)[2]
  result <- array(dim=c(L,L,K))
  for (k in 1:K) {
    result[,,k] <- generate.statistic(counts[k,])
  }
  result
}

generate.statistic <- function(counts) {
  result <- outer(rep(1,length(counts)),counts)
  result[lower.tri(result)] <- 0
  for (i in 1:length(counts)) {
    result[i,] <- cumsum(result[i,])
  }
  result
}

adaptive.sampling <- function(counts, alpha, beta, gamma, ...) {
  L <- dim(counts)[1]
  K <- dim(counts)[3]
  storage.mode(counts) <- "double"
  storage.mode(alpha)  <- "double"
  storage.mode(beta)   <- "double"
  storage.mode(gamma)  <- "double"
  
  options <- make.options(...)

  result <- as.list(.Call("adaptive_sampling",
                counts, alpha, beta, gamma, options))
  attr(result, 'class') <- 'sampling'
  result
}

plotResult <- function(result) {
  expectation <- result$moments[1,]
  utility <- result$utility
  par(mfrow=c(3,1))
  plot(expectation , type="l", ylab=expression(p["x,s"]), xaxs="i", ylim=c(0,1), main="Expectation")
  plot(utility, type="l", ylab=expression(U(x)), xaxs="i", main="Utility")
  barplot(result$mpost, space=0.0, xlab="m", ylab=expression(P(M==m ~~ "|" ~~ Y^X)), xaxs="i", ylim=c(0,1), main="Model posterior")
}

plot.sampling <- function(samplingresult, ...)
	plotResult(unclass(samplingresult), ...)

sampling.demo <- function() {
  L = 6 # number of stimuli
  K = 2 # number of responses
  counts_success <- c(2,3,2,4,7,7)
  counts_failure <- c(8,7,7,6,3,2)
  counts <- count.statistic(t(matrix(c(counts_success, counts_failure), L)))
  alpha_success  <- c(1,1,1,1,1,1)
  alpha_failure  <- c(1,1,1,1,1,1)
  alpha  <- default.alpha(t(matrix(c(alpha_success, alpha_failure), L)))
  beta   <- default.beta(L)
  gamma  <- default.gamma(L)

  adaptive.sampling(counts, alpha, beta, gamma)
}
