
# forward recursion
################################################################################

phi <- function(parameter, result, i, j) {
  a1 <- parameter$alpha.success[i]
  a2 <- parameter$alpha.failure[i]
  n1 <- parameter$alpha.success[i] + sum(parameter$success[i:j])
  n2 <- parameter$alpha.failure[i] + sum(parameter$failure[i:j])

  if (i == 1) {
    return(beta(n1, n2)/beta(a1, a2))
  }
  else {
#    return(parameter$rho * result[i-1,j] + (1-parameter$rho) * beta(n1, n2)/beta(a1, a2) * result[i-1,i-1])
    return(parameter$rho/result[i-1,i-1] * result[i-1,j] + (1-parameter$rho) * beta(n1, n2)/beta(a1, a2))
  }
}
phi.init <- function(parameter) {
  result <- matrix(0, parameter$L, parameter$L)
  for (i in 1:parameter$L) {
    for (j in i:parameter$L) {
      result[i,j] <- phi(parameter, result, i, j)
    }
  }
  return(result)
}

f.forward.prime <- function(parameter, i,j, h) {
  a1 <- parameter$alpha.success[i]
  a2 <- parameter$alpha.failure[i]
  n1 <- parameter$alpha.success[i] + sum(parameter$success[i:j])
  n2 <- parameter$alpha.failure[i] + sum(parameter$failure[i:j])
  if(i == 1) {
    return(beta(n1, n2)/beta(a1, a2)*h(i, j, n1,n2))
  }
  else {
#    return(parameter$rho * f.forward.prime(parameter, i-1,j, h) + (1-parameter$rho) * beta(n1, n2)/beta(a1, a2)*parameter$phi[i-1,i-1]*h(i, j, n1,n2))
    return(parameter$rho/parameter$phi[i-1,i-1] * f.forward.prime(parameter, i-1,j, h) + (1-parameter$rho) * beta(n1, n2)/beta(a1, a2)*h(i, j, n1,n2))
  }
}
f.forward <- function(parameter, j, h) {
  return (f.forward.prime(parameter, j,j, h)/parameter$phi[j,j])
}
f.forward.expectation <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1/(n1+n2), n2/(n1+n2))
  return(f.forward(parameter, j, h))
}
f.forward.stddev <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1*(n1+n2 - n1)/((n1+n2)^2*(n1+n2+1)), n2*(n1+n2 - n2)/((n1+n2)^2*(n1+n2+1)))
  return(sqrt(f.forward(parameter, j, h)))
}

# backward recursion
################################################################################

psi <- function(parameter, result, j, k) {
  a1 <- parameter$alpha.success[j]
  a2 <- parameter$alpha.failure[j]
  n1 <- parameter$alpha.success[j] + sum(parameter$success[j:k])
  n2 <- parameter$alpha.failure[j] + sum(parameter$failure[j:k])

  if (k == parameter$L) {
    return (beta(n1, n2)/beta(a1, a2))
  }
  else {
    return (parameter$rho * result[j,k+1] + (1-parameter$rho) * beta(n1, n2)/beta(a1, a2) * result[j+1,k+1])
  }
}
psi.init <- function(parameter, result) {
  result <- matrix(0, parameter$L, parameter$L)
  for (j in parameter$L:1) {
    for (k in parameter$L:j) {
      result[j,k] <- psi(parameter, result, j, k)
    }
  }
  return(result)
}

# forward and backward recursion
################################################################################

g <- function(parameter, j, k, h) {
  if (k == parameter$L) {
    return (f.forward.prime(parameter, j,k, h))
  }
  else {
    return (parameter$rho * g(parameter, j, k+1, h) + (1-parameter$rho)*f.forward.prime(parameter, j,k, h)*parameter$psi[j+1,k+1])
  }
}
g.norm <- function(parameter, j, k) {
  if (k == parameter$L) {
    return (parameter$phi[j,k])
  }
  else {
    return (parameter$rho * g.norm(parameter, j, k+1) + (1-parameter$rho)*parameter$phi[j,k]*parameter$psi[j+1,k+1])
  }
}
f <- function(parameter, j,h) {
  return (g(parameter, j,j,h)/g.norm(parameter, j,j))
}
f.expectation <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1/(n1+n2), n2/(n1+n2))
  return(f(parameter, j, h))
}
f.stddev <- function(parameter, j) {
  h <- function(i, j, n1, n2) c(n1*(n1+n2 - n1)/((n1+n2)^2*(n1+n2+1)), n2*(n1+n2 - n2)/((n1+n2)^2*(n1+n2+1)))
  return(sqrt(f(parameter, j, h)))
}
f.density <- function(parameter, j, theta) {
  h <- function(i, j, n1, n2) 1/beta(n1, n2)*theta^n1*(1-theta)^n2
  return(f(parameter, j, h))
}

# init
################################################################################

parameter.init <- function(parameter, extended=TRUE) {
  parameter$L <- length(parameter$success)
  parameter$alpha.success <- rep(1, parameter$L)
  parameter$alpha.failure <- rep(1, parameter$L)
  parameter$phi <- phi.init(parameter)
  if (extended) {
    parameter$psi     <- psi.init(parameter)
    parameter$expect  <- sapply(1:parameter$L, function(j) f.expectation(parameter, j))
    parameter$stddev  <- sapply(1:parameter$L, function(j) f.stddev(parameter, j))
    parameter$utility <- sapply(1:parameter$L, function(j) f.utility(parameter, j))
  }
  return(parameter)
}

# utility
################################################################################

f.utility.rec.prime <- function(parameter, result, i, j, h) {
  a1 <- parameter$alpha.success[i]
  a2 <- parameter$alpha.failure[i]
  n1 <- parameter$alpha.success[i] + sum(parameter$success[i:j])
  n2 <- parameter$alpha.failure[i] + sum(parameter$failure[i:j])
  
  if (i == 1) {
    return(beta(n1, n2)/beta(a1, a2)*h(i, j, n1,n2))
  }
  else {
    return(parameter$rho * result[i-1,j] + (1-parameter$rho) * beta(n1, n2)/beta(a1, a2)*result[i-1,i-1]*h(i, j, n1,n2))
  }
}
f.utility.rec <- function(parameter, i, j, h) {
  result <- matrix(0, parameter$L, parameter$L)
  for (i in 1:parameter$L) {
    for (j in i:parameter$L) {
      result[i,j] <- f.utility.rec.prime(parameter, result, i, j, h)
    }
  }
  return(result[parameter$L,parameter$L])
}

f.utility.prime <- function(parameter, j, l) {
  L <- parameter$L
  parameter.new <- parameter
  if (l == 1) {
    h <- function(p, q, n1, n2) {
      if (p <= j && j <= q) {
        return(digamma(n1) - digamma(n1 + n2))
      }
      else {
        return(1)
      }
    }
    parameter.new$success[j] <- parameter.new$success[j] + 1
  }
  else {
    h <- function(p, q, n1, n2) {
      if (p <= j && j <= q) {
        return(digamma(n2) - digamma(n1 + n2))
      }
      else {
        return(1)
      }
    }
    parameter.new$failure[j] <- parameter.new$failure[j] + 1
  }
  parameter.new <- parameter.init(parameter.new, FALSE)

  result <- log(parameter$phi[1,1])-log(parameter.new$phi[1, 1])
  ml2    <- parameter.new$phi[1, 1]
  for (i in 2:L) {
    result <- result + log(parameter$phi[i,i])-log(parameter.new$phi[i, i])
    ml2    <- ml2*parameter.new$phi[i, i]
  }
  result <- result + f.utility.rec(parameter.new, L, L, h)/ml2
  
  return(result)
}

f.utility <- function(parameter, j, l) {
  return(parameter$expect[1,j]*f.utility.prime(parameter, j, 1) + parameter$expect[2,j]*f.utility.prime(parameter, j, 2))
}

# plotting
################################################################################

plot.density <- function(parameter) {
  K <- 50
  result <- matrix(0, parameter$L, K)
  for (i in 1:parameter$L) {
    for (j in 1:K) {
      result[i, j] <- f.density(parameter, i, j/K)
    }
  }
  image(1:parameter$L,1:K/K,result, col=gray(1-1:K/K/2), xlab="x", ylab=expression(hat(theta)[1]))
}

plot.sampling <- function(parameter, groundtruth=c(), density=FALSE) {

  par(mfrow=c(3,1), mar=c(4.3, 4.3, 2, 1))

  if (density) {
    plot.density(parameter)
  }
  else {
    plot(c(1,parameter$L),c(0,1.0), type='n', xaxp=c(1,parameter$L,parameter$L-1))
    grid(parameter$L,parameter$L)
  }

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
  if (n > 0) {
    parameter <- parameter.init(parameter)
    for (i in 1:n) {
      parameter <- sample.prime(parameter)
      print(parameter$utility)
      if (video) {
        png(filename=sprintf("plot_%03d.png", i), width = 800, height = 600)
        plot.sampling(parameter, groundtruth, TRUE)
        dev.off()
      }
      else {
        plot.sampling(parameter, groundtruth)
      }
    }
  }
  return(parameter)
}

# data and parameters
################################################################################

groundtruth <- c(0.977895, 0.959606, 0.927331, 0.872769, 0.786948, 0.666468, 0.523201, 0.388603, 0.307012, 0.327954, 0.481978, 0.705735, 0.924494, 0.956063, 0.986731, 0.996743, 0.999621)
#groundtruth <- c(0.983707, 0.977895, 0.970075, 0.959606, 0.945683, 0.927331, 0.903428, 0.872769, 0.834219, 0.786948, 0.730763, 0.666468, 0.59614, 0.523201, 0.452233, 0.388603, 0.338142, 0.307012, 0.301628, 0.327954, 0.38205, 0.481978, 0.593685, 0.705735, 0.901878, 0.924494, 0.944262, 0.956063, 0.975398, 0.986731, 0.993179, 0.996743, 0.998648, 0.999621, 1.00000)

parameter   <- list(rho = 0.4,
                    success = rep(0, length(groundtruth)),
                    failure = rep(0, length(groundtruth))
                    )

parameter <- sample(parameter, 250, video=FALSE)

#mencoder mf://plot_*.png -mf type=png:fps=4 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o plot.avi
