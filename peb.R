library("Matrix")
library("dplyr")
library("collapse")
library("mvtnorm")
library("Rfast")
library("Rcpp")

sourceCpp("peb.cpp")


tournament_matrix <- function(tournament)
{
  teams <- levels(tournament$home)
  nteams <- length(teams)
  ihome <- match(tournament$home, teams)
  iaway <- match(tournament$away, teams)
  ngames <- nrow(tournament)
  dims <- c(ngames, nteams)
  dimnames <- list(seq_len(ngames), teams)
  ones <- rep(1, ngames)
  sparseMatrix(
    i = rep(seq_len(ngames), 2L),
    j = c(ihome, iaway),
    x = c(ones, -ones),
    dims = dims,
    dimnames = dimnames)
}

lambda2tau <- function(lambda) {
  2 / pi * asin(1 / (lambda + 2))
}

tau2lambda <- function(tau) {
  (1 - 2 * sin(0.5 * pi * tau)) / sin(0.5 * pi * tau)
}

optimR <- function(par, fun, grad, method = "BFGS", ...) {
  op <- optim(par, fn = fun, gr = grad, method = method, ...)
  ans <- list()
  ans$converged <- (op$convergence == 0)
  if (ans$converged)
    ans$par <- op$par
  else
    ans$par <- rep(NA, length(par))
  ans$value <- op$value
  ## bye bye
  ans
}

compute_tau <- function(stats, delta, gamma, a, eps)
{
  ## closed-form solution if no home effect and no ties
  if (delta == 0 && gamma == 0) {
    tau <- (stats[1,1] + stats[3,3] - stats[1,3] - stats[3,1]) /
      (stats[1,1] + stats[3,3] + stats[1,3] + stats[3,1] + a)
    tau <- max(eps, min(tau, 1/3 - eps))
  }
  else {
    ##
    probs <- matrix(0.0, 3, 3)
    ## penalized negative pairwise log-likelihood
    pair_lik <- function(tau) {
      rho <- sin(0.5 * pi * tau)
      ## compute bivariate probs
      cuts <- c(-Inf, -gamma, gamma, Inf) - delta
      R <- matrix(c(1, rho, rho, 1), nrow = 2L, ncol = 2L)
      for (i in 1:3) {
        for (j in i:3) {
          probs[i, j] <- pmvnorm(lower = c(cuts[i], cuts[j]),
                                 upper = c(cuts[i + 1], cuts[j + 1]),
                                 corr = R)
          if (i < j) probs[j, i] <- probs[i, j]
        }
      }
      ## bye bye
      -sum(log(probs[probs > 0]) * stats[probs > 0]) - a * log(1 - tau ^ 2)
    }
    tau <- optimise(pair_lik, interval = c(eps, 1/3 - eps))$minimum
  }
  ## bye bye
  tau
}

peb <- function(tournament, lambda = NULL, home_effect = FALSE, a = 2 * nlevels(tournament$home), eps = 1e-4, ...)
{
  ## extract the data
  teams <- levels(tournament$home)
  nteams <- length(teams)
  home <- tournament$home
  away <- tournament$away
  y <- tournament$outcome ## 1 = visitor wins; 2 = tie; 3 = host wins

  ## answer
  ans <- list()
  ans$ties <- any(y == 2)
  ans$home_effect <- home_effect

  ## 1- estimation of the tuning parameter
  ## 1.1- estimation of delta and gamma
  ngames <- length(y)
  marginals <- c(sum(y == 1), sum(y == 2), sum(y == 3))
  q1 <- qnorm(marginals[1] / (ngames + 1))
  q12 <- qnorm(1 - marginals[3] / (ngames + 1))
  if (ans$home_effect)
    ans$delta <- -0.5 * (q1 + q12)
  else
    ans$delta <- 0.0
  if (ans$ties)
    ans$gamma <- max(0.0, 0.5 * (q12 - q1))
  else
    ans$gamma <- 0.0
  if (is.null(lambda)) {
    ## 1.2- estimation of lambda
    stats <- compute_stats_rcpp(y, match(home, teams), match(away, teams),
                                home_effect)
    ans$tau <- compute_tau(stats, ans$delta, ans$gamma, a, eps)
    ans$lambda <- tau2lambda(ans$tau)
  }
  else {## lambda is given
    ans$lambda <- lambda
    ans$tau <- lambda2tau(lambda)
  }
  ## 2- estimation of the strength parameters
  ans$X <- tournament_matrix(tournament)
  cuts <- c(-Inf, -ans$gamma, ans$gamma, Inf) - ans$delta
  ## penalized negative log-likelihood
  ans$lik <- function(mu) {
    eta <- mu[home] - mu[away]
    probs <- pnorm(cuts[y + 1] - eta) - pnorm(cuts[y] - eta)
    ## bye bye
    -sum(log(probs)) + 0.5 * ans$lambda * sum(mu * mu)
  }
  ## gradient of the penalized negative log-likelihood
  ans$score <- function(mu) {
    eta <- mu[home] - mu[away]
    probs <- pnorm(cuts[y + 1] - eta) - pnorm(cuts[y] - eta)
    dprobs <- dnorm(cuts[y + 1] - eta) - dnorm(cuts[y] - eta)
    ## bye bye
    as.vector(crossprod(ans$X, dprobs / probs)) + ans$lambda * mu
  }
  ## go!
  op <- optimR(rep(0.0, nteams), ans$lik, ans$score, ...)
  ans$mu <- op$par
  names(ans$mu) <- teams
  ## bye bye
  ans
}

compute_logscore <- function(object, newdata, link = c("logit", "probit"),
                             small = 1e-8)
{
  link <- match.arg(link)
  y <- newdata$outcome
  eta <- object$mu[newdata$home] - object$mu[newdata$away]
  cuts <- c(-Inf, -object$gamma, object$gamma, Inf) - object$delta
  F <- function(x)
    if (link == "logit") plogis(x) else pnorm(x)
  probs <- F(cuts[y + 1] - eta) - F(cuts[y] - eta)
  probs[probs < small] <- small
  ## bye bye
  -mean(log(probs))
}

compute_skillscore <- function(object, newdata, link = c("logit", "probit"),
                               ref_probs, small = 1e-8)
{
  logscore <- compute_logscore(object, newdata, link, small)
  ref_logscore <- -sum(ref_probs * log(ref_probs))
  ## bye bye
  1 - logscore / ref_logscore
}

brmle <- function(tournament, br = TRUE, ...) {
  ans <- list()
  outcome <- as.factor(tournament$outcome)
  X <- as.matrix(tournament_matrix(tournament))
  foot_data <- data.frame(outcome, X)
  teams <- colnames(foot_data)[-1L]
  nteams <- length(teams)
  fo <- paste("outcome", paste0(teams[-1L], collapse=" + "), sep = " ~ ")
  ## go
  fit <- try(bpolr(fo, data = foot_data, method = if (br) "BR" else "ML", ...), silent = TRUE)
  if(inherits(brmle, "try-error")) {
    ans$mu <- rep(NA, nteams - 1)
    ans$gamma <- NA
    ans$delta <- NA
  }
  else {
    ans$mu <- c(0.0, fit$beta)
    ans$mu <- ans$mu - mean(ans$mu)
    ans$gamma <- unname(0.5 * (fit$alpha[2] - fit$alpha[1]))
    ans$delta <- unname(-0.5 * (fit$alpha[1] + fit$alpha[2]))
  }
  names(ans$mu) <- teams
  ## bye bye
  ans
}





