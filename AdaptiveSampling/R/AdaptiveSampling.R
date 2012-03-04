default_options <- function() {
  env <- environment()
  env$n_moments             <- 2
  env$model_posterior       <- 1
  env$bprob                 <- 1 
  env$utility               <- 1
  env$kl_component          <- 1
  env$kl_multibin           <- 0
  env$effective_counts      <- 0
  env$marginal              <- 1
  env$marginal_step         <- 0.01
  env$marginal_range        <- c(0,1)
  env$epsilon               <- 0.00001
  env$threads               <- 1
  env$stacksize             <- 256*1024
  env$algorithm             <- 0
  env$which                 <- 0
  env$samples               <- c(100, 2000)

  env
}

default_alpha <- function(alpha) {
  K <- dim(alpha)[1]
  L <- dim(alpha)[2]
  result <- array(dim=c(L,L,K))
  for (k in 1:K) {
    result[,,k] <- generate_alpha(alpha[k,])
  }
  result
}

generate_alpha <- function(alpha) {
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

default_beta <- function(n) {
  return(rep(1, n))
}

default_gamma <- function(n) {
  result <- outer(rep(1,n),rep(1,n))
  result[lower.tri(result)] <- 0
  result
}

count_statistic <- function(counts) {
  K <- dim(counts)[1]
  L <- dim(counts)[2]
  result <- array(dim=c(L,L,K))
  for (k in 1:K) {
    result[,,k] <- generate_statistic(counts[k,])
  }
  result
}

generate_statistic <- function(counts) {
  result <- outer(rep(1,length(counts)),counts)
  result[lower.tri(result)] <- 0
  for (i in 1:length(counts)) {
    result[i,] <- cumsum(result[i,])
  }
  result
}

adaptive_sampling <- function(counts, alpha, beta, gamma, options) {
  L <- dim(counts)[1]
  K <- dim(counts)[3]
  storage.mode(counts) <- "double"
  storage.mode(alpha)  <- "double"
  storage.mode(beta)   <- "double"
  storage.mode(gamma)  <- "double"

  .Call("adaptive_sampling",
     counts, alpha, beta, gamma, options)
}

plot_result <- function(result) {
  expectation <- result$moments[1,]
  utility <- result$utility
  par(mfrow=c(3,1))
  plot(expectation , type="l", ylab=expression(p["x,s"]), xaxs="i", ylim=c(0,1), main="Expectation")
  plot(utility, type="l", ylab=expression(U(x)), xaxs="i", main="Utility")
  barplot(result$mpost, space=0.0, xlab="m", ylab=expression(P(M==m ~~ "|" ~~ Y^X)), xaxs="i", ylim=c(0,1), main="Model posterior")
}

sampling_demo <- function() {
  L = 6 # number of stimuli
  K = 2 # number of responses
  counts_success <- c(2,3,2,4,7,7)
  counts_failure <- c(8,7,7,6,3,2)
  counts <- count_statistic(t(matrix(c(counts_success, counts_failure), L)))
  alpha_success  <- c(1,1,1,1,1,1)
  alpha_failure  <- c(1,1,1,1,1,1)
  alpha  <- default_alpha(t(matrix(c(alpha_success, alpha_failure), L)))
  beta   <- default_beta(L)
  gamma  <- default_gamma(L)

  options <- default_options()

  adaptive_sampling(counts, alpha, beta, gamma, options)
}
