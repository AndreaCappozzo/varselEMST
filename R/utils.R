estep_EMST <- function(data, fitm, F_subset){
  # this is necessary because I will use only the sigma ( restricted to F_subset), whereas estep from Mclust uses the decomposition in scale, shape, orientation
  log_numerator <-
    sapply(1:fitm$G, function(g)
      log(fitm$parameters$pro[g]) + mvnfast::dmvn(
        X = as.matrix(data[,F_subset]),
        mu = fitm$parameters$mean[F_subset, g],
        sigma = fitm$parameters$variance$sigma[F_subset,F_subset , g],
        log = T
      ))

  log_denominator <- mvnfast::dmixn(
    X = as.matrix(data[, F_subset, drop = FALSE]),
    mu = t(fitm$parameters$mean[F_subset, , drop = FALSE]),
    sigma = purrr::array_tree(fitm$parameters$variance$sigma[F_subset, F_subset, , drop =
                                                               FALSE], 3),
    w = fitm$parameters$pro, log = T
  )

  z <- exp(log_numerator-log_denominator)
  res <- list()
  res$z <- z # to be coherent with estep from mclust
  res
}

H_F_subset <- function(S_R_w, S_R, fitm, col_subset) {
  # (5.20) Ritter
  log_det_S_R_w_F <-
    apply(S_R_w[col_subset, col_subset, , drop = FALSE], 3, function(Sigma)
      determinant(Sigma, logarithm = TRUE)$modulus)
  log_det_S_R_F <-
    determinant(S_R[col_subset, col_subset, drop = FALSE], logarithm = TRUE)$modulus

  as.vector(log_det_S_R_w_F%*%fitm$parameters$pro) - log_det_S_R_F
}

# Density Singular multinormal distribution
dsmvnorm <- function(x, mean = rep(0, p), sigma = diag(p), log = FALSE, eigen_sigma, eigen_tol=1e-12){
  # eigen_sigma might be considered redundant, in my specific case I insert it as a function arg for avoiding  computing it again
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean)))
      dim(mean) <- NULL
    if (length(mean) != p)
      stop("mean and sigma have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
  }
  rank_sigma <- sum(eigen_sigma>eigen_tol) # number of positive eigenvalues

  # dec <- tryCatch(chol(sigma), error = function(e) e)
  # if (inherits(dec, "error")) {
  #   x.is.mu <- colSums(t(x) != mean) == 0
  #   logretval <- rep.int(-Inf, nrow(x))
  #   logretval[x.is.mu] <- Inf
  # }
  # else {
  rss <- mahalanobis(
    x = x,
    center = mean,
    cov = MASS::ginv(sigma),
    inverted = TRUE
  )

  logretval <-
    -0.5 * sum(log(eigen_sigma[1:rank_sigma])) - 0.5 * rank_sigma * log(2 *
                                                                          pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  if (log)
    logretval
  else exp(logretval)
}

#
# fitness_b_with_penalty <- function(string,
#                                    S_R_w,
#                                    S_R,
#                                    fitm,
#                                    n_relevant_variables) {
#   inc <- which(string == 1)
#   obj <- -H_F_subset(
#     col_subset = inc,
#     S_R_w = S_R_w,
#     S_R = S_R,
#     fitm = fitm
#   )
#   penalty_upper <-
#     ifelse(sum(string) > n_relevant_variables, -Inf, 0)
#   penalty_lower <-
#     ifelse(sum(string) < n_relevant_variables, -Inf, 0)
#   obj + penalty_lower + penalty_upper
# }
