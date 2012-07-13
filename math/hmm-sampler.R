
# forward recursion (marginal posterior)
################################################################################

phi <- function(parameter, phi.result, j, jp, h) {
  a1 <- parameter$alpha.success[1]
  a2 <- parameter$alpha.failure[1]
  n1 <- parameter$alpha.success[1] + sum(parameter$success[1:jp])
  n2 <- parameter$alpha.failure[1] + sum(parameter$failure[1:jp])

  result <- parameter$rho^(jp-1)*beta(n1, n2)/beta(a1, a2)*h(1, jp, n1, n2)
  if (j > 1) {
    for (k in 1:(j-1)) {
      a1 <- parameter$alpha.success[k+1]
      a2 <- parameter$alpha.failure[k+1]
      n1 <- parameter$alpha.success[k+1] + sum(parameter$success[(k+1):jp])
      n2 <- parameter$alpha.failure[k+1] + sum(parameter$failure[(k+1):jp])
      result <- result + (1-parameter$rho)*parameter$rho^(jp-k-1)*phi.result[k]*beta(n1, n2)/beta(a1, a2)*h(k+1, jp, n1, n2)
    }
  }
  return(result)
}
phi.init <- function(parameter) {
  result <- rep(0, parameter$L)
  h <- function(i, j, n1, n2) 1
  for (j in 1:parameter$L) {
    result[j] <- phi(parameter, result, j, j, h)
  }
  return(result)
}
f.forward <- function(parameter, j, h) {
  return (phi(parameter, parameter$phi, j, j, h)/parameter$phi[j])
}
f.forward.expectation <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1/(n1+n2), n2/(n1+n2))
  return(f.forward(parameter, j, h))
}
f.forward.stddev <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1*(n1+n2 - n1)/((n1+n2)^2*(n1+n2+1)), n2*(n1+n2 - n2)/((n1+n2)^2*(n1+n2+1)))
  return(sqrt(f.forward(parameter, j, h)))
}

# backward recursion (marginal posterior)
################################################################################

psi <- function(parameter, psi.result, j, h) {
  a1 <- parameter$alpha.success[j]
  a2 <- parameter$alpha.failure[j]
  n1 <- parameter$alpha.success[j] + sum(parameter$success[j:parameter$L])
  n2 <- parameter$alpha.failure[j] + sum(parameter$failure[j:parameter$L])

  result <- parameter$rho^(parameter$L-j)*beta(n1, n2)/beta(a1, a2)*h(j, parameter$L, n1, n2)
  if (j+1 <= parameter$L) {
    for (k in (j+1):parameter$L) {
      a1 <- parameter$alpha.success[j]
      a2 <- parameter$alpha.failure[j]
      n1 <- parameter$alpha.success[j] + sum(parameter$success[j:(k-1)])
      n2 <- parameter$alpha.failure[j] + sum(parameter$failure[j:(k-1)])
      result <- result + (1-parameter$rho)*parameter$rho^(k-1-j)*psi.result[k]*beta(n1, n2)/beta(a1, a2)*h(j, k-1, n1, n2)
    }
  }
  return(result)
}
psi.init <- function(parameter) {
  result <- rep(0, parameter$L)
  h <- function(i, j, n1, n2) 1
  for (j in parameter$L:1) {
    result[j] <- psi(parameter, result, j, h)
  }
  return(result)
}

# forward and backward recursion (marginal posterior)
################################################################################

g <- function(parameter, j, h) {
  result <- phi(parameter, parameter$phi, j, parameter$L, h)
  if (j+1 <= parameter$L) {
    for (k in (j+1):parameter$L) {
      result <- result + (1-parameter$rho)*phi(parameter, parameter$phi, j, k-1, h)*parameter$psi[k]
    }
  }
  return(result/parameter$phi[j])
}
g.norm <- function(parameter, j) {
  return(parameter$phi[parameter$L]/parameter$phi[j])
}
f <- function(parameter, j, h) {
  return (g(parameter, j, h)/g.norm(parameter, j))
}
f.expectation <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1/(n1+n2), n2/(n1+n2))
  return(f(parameter, j, h))
}
f.stddev <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1*(n1+n2 - n1)/((n1+n2)^2*(n1+n2+1)), n2*(n1+n2 - n2)/((n1+n2)^2*(n1+n2+1)))
  return(sqrt(f(parameter, j, h)))
}

# mixture chain (reverse parameters)
################################################################################

parameter.reverse <- function(parameter.1) {
  parameter.2 <- list()
  parameter.2$success <- rev(parameter.1$success)
  parameter.2$failure <- rev(parameter.1$failure)
  if (is.element('alpha.success', names(parameter.1))) {
    parameter.2$alpha.success <- rev(parameter.1$alpha.success)
  }
  if (is.element('alpha.failure', names(parameter.1))) {
    parameter.2$alpha.failure <- rev(parameter.1$alpha.failure)
  }
  parameter.2$rho <- parameter.1$rho
  
  return(parameter.2)
}

# init parameters and compute posterior
################################################################################

parameter.init <- function(parameter.1, extended=TRUE) {
  parameter.1$L <- length(parameter.1$success)
  if (!is.element('alpha.success', names(parameter.1))) {
    parameter.1$alpha.success <- rep(1, parameter.1$L)
  }
  if (!is.element('alpha.failure', names(parameter.1))) {
    parameter.1$alpha.failure <- rep(1, parameter.1$L)
  }
  parameter.1$phi <- phi.init(parameter.1)
  if (extended) {
    # first chain
    parameter.1$psi     <- psi.init(parameter.1)
    parameter.1$expect  <- sapply(1:parameter.1$L, function(j) f.expectation(parameter.1, j))
    parameter.1$stddev  <- sapply(1:parameter.1$L, function(j) f.stddev(parameter.1, j))
    parameter.1$utility <- sapply(1:parameter.1$L, function(j) f.utility(parameter.1, j))

    # second chain
    parameter.2         <- parameter.reverse(parameter.1)
    parameter.2         <- parameter.init(parameter.2, FALSE)
    parameter.2$psi     <- psi.init(parameter.2)
    parameter.2$expect  <- sapply(1:parameter.2$L, function(j) f.expectation(parameter.2, j))
    parameter.2$stddev  <- sapply(1:parameter.2$L, function(j) f.stddev(parameter.2, j))
    parameter.2$utility <- sapply(1:parameter.2$L, function(j) f.utility(parameter.2, j))

    # mix results
    parameter.1$expect  <- sapply(1:parameter.1$L, function(j) (parameter.1$expect[,j]+parameter.2$expect [,parameter.1$L-j+1])/2)
    parameter.1$stddev  <- sapply(1:parameter.1$L, function(j) (parameter.1$stddev[,j]+parameter.2$stddev [,parameter.1$L-j+1])/2)
    parameter.1$utility <- sapply(1:parameter.1$L, function(j) (parameter.1$utility[j]+parameter.2$utility[ parameter.1$L-j+1])/2)
  }
  return(parameter.1)
}

# utility (on the full joint posterior)
################################################################################

f.utility.prime <- function(parameter, j, l) {
  L <- parameter$L
  parameter.new <- parameter
  if (l == 1) {
    h <- function(p, q, n1, n2) {
      return(digamma(n1) - digamma(n1 + n2))
    }
    parameter.new$success[j] <- parameter.new$success[j] + 1
  }
  else {
    h <- function(p, q, n1, n2) {
      return(digamma(n2) - digamma(n1 + n2))
    }
    parameter.new$failure[j] <- parameter.new$failure[j] + 1
  }
  parameter.new <- parameter.init(parameter.new, FALSE)

  result <- log(parameter$phi[L])-log(parameter.new$phi[L])
  result <- result + f(parameter.new, j, h)
  return(result)
}

f.utility <- function(parameter, j, l) {
  return(parameter$expect[1,j]*f.utility.prime(parameter, j, 1) + parameter$expect[2,j]*f.utility.prime(parameter, j, 2))
}

# plotting
################################################################################

plot.sampling <- function(parameter, groundtruth=c(), density=FALSE) {

  par(mfrow=c(3,1), mar=c(4.3, 4.3, 2, 1))

  plot(c(1,parameter$L),c(0,1.0), type='n', xaxp=c(1,parameter$L,parameter$L-1), xlab="x", ylab=expression(hat(theta)[1]))
  grid(parameter$L,parameter$L)

  if (length(groundtruth) != 0) {
    lines(1:parameter$L, groundtruth, type='l', col='blue', lwd=1.5)
  }
  lines(1:parameter$L, parameter$expect[1,], type='l', col='red', lwd=2.5)
  lines(1:parameter$L, parameter$expect[1,]-parameter$stddev[1,], type='l', col='black', lwd=1)
  lines(1:parameter$L, parameter$expect[1,]+parameter$stddev[1,], type='l', col='black', lwd=1)

  plot(c(1,parameter$L),c(min(parameter$utility),max(parameter$utility)), type='n', xaxp=c(1,parameter$L,parameter$L-1), xlab="x", ylab="Utility")
  grid(parameter$L,parameter$L)
  lines(1:parameter$L, parameter$utility, type='l', col='blue', lwd=2.5)

  counts <- mapply(sum, parameter$success, parameter$failure)
  barplot(counts, ylim=c(0,30), xlab="x", ylab="Counts")
}

# sampling and ground truth
################################################################################

library(nnet) # which.is.max

get.sample <- function(j) {
  if (runif(1, 0, 1) <= groundtruth[j]) {
    return(1)
  }
  else {
    return(2)
  }
}

sample.prime <- function(parameter) {
  j <- which.is.max(parameter$utility)
  l <- get.sample(j)
  print(j)
  if (l == 1) {
    parameter$success[j] <- parameter$success[j] + 1
  }
  else {
    parameter$failure[j] <- parameter$failure[j] + 1
  }
  parameter <- parameter.init(parameter)

  return(parameter)
}

sample <- function(parameter, n, video=FALSE) {
  parameter <- parameter.init(parameter)
  if (video) {
    png(filename=sprintf("plot_%03d.png", 0), width = 800, height = 600)
    plot.sampling(parameter, groundtruth, FALSE)
    dev.off()
  }
  else {
    plot.sampling(parameter, groundtruth)
  }
  if (n > 0) {
    for (i in 1:n) {
      parameter <- sample.prime(parameter)
      print(parameter$utility)
      if (video) {
        png(filename=sprintf("plot_%03d.png", i), width = 800, height = 600)
        plot.sampling(parameter, groundtruth, FALSE)
        dev.off()
      }
      else {
        plot.sampling(parameter, groundtruth, TRUE)
      }
    }
  }
  return(parameter)
}

# data and parameters
################################################################################

groundtruth <- c(0.977895, 0.959606, 0.927331, 0.872769, 0.786948, 0.666468, 0.523201, 0.388603, 0.307012, 0.327954, 0.481978, 0.705735, 0.924494, 0.956063, 0.986731, 0.996743, 0.999621)
#groundtruth <- c(0.983707, 0.977895, 0.970075, 0.959606, 0.945683, 0.927331, 0.903428, 0.872769, 0.834219, 0.786948, 0.730763, 0.666468, 0.59614, 0.523201, 0.452233, 0.388603, 0.338142, 0.307012, 0.301628, 0.327954, 0.38205, 0.481978, 0.593685, 0.705735, 0.901878, 0.924494, 0.944262, 0.956063, 0.975398, 0.986731, 0.993179, 0.996743, 0.998648, 0.999621, 1.00000)

parameter <- list(rho = 0.4,
                  success = rep(0, length(groundtruth)),
                  failure = rep(0, length(groundtruth))
                  )
parameter <- parameter.init(parameter, FALSE)

parameter$alpha.success[ 1] = 10
parameter$alpha.success[17] = 10

parameter <- sample(parameter, 100, video=FALSE)

#mencoder mf://plot_*.png -mf type=png:fps=4 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o plot.avi
