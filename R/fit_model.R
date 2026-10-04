
#taufunc()
#penalty_eta()
#grad_penalty_eta()
#fit_model()

taufunc <- function(y, lambda) {
  1 / (1 + (1 - y)^(1 / lambda))
}


penalty_eta <- function(eta, nr, sigma0 = 2, sigma = 1, eps = 1e-12) {
  p0 <- min(max(1 / nr, eps), 1 - eps)
  mu0 <- stats::qlogis(p0)

  eta0 <- eta[1]
  eta_nonint <- eta[-1]

  -((eta0 - mu0)^2) / (2 * sigma0^2) -
    sum(eta_nonint^2) / (2 * sigma^2)
}


grad_penalty_eta <- function(eta, nr, sigma0 = 2, sigma = 1, eps = 1e-12) {
  p0 <- min(max(1 / nr, eps), 1 - eps)
  mu0 <- stats::qlogis(p0)

  grad <- numeric(length(eta))
  grad[1] <- -(eta[1] - mu0) / sigma0^2
  grad[-1] <- -eta[-1] / sigma^2

  grad
}


#' Fit the penalized composite global likelihood model
#'
#' Estimates the exponential transmission-effect rate parameter and logistic
#' enrichment coefficients using nested profile optimization.
#'
#' @param data Output from [make_data()].
#' @param lambda_interval Search interval for lambda.
#' @param tol_loglambda Optimization tolerance on log lambda.
#' @param eta_start Optional starting values for eta.
#' @param use_H Logical. If TRUE, include age-at-onset unusualness score.
#' @param sigma0 Prior standard deviation for eta0.
#' @param sigma Prior standard deviation for non-intercept eta coefficients.
#' @param TRACE_TIMING Logical. If TRUE, record timing information.
#' @param eps Numerical guard for likelihood/ initial probabilities.
#' @param rel.tol Relative tolerance for numerical integration.
#'
#' @return An object of class `"fbf_fit"`.
#' @export

fit_model <- function(data,
                      lambda_interval = c(1e-4, 50),
                      tol_loglambda = 1e-3,
                      eta_start = NULL,
                      use_H = TRUE,
                      sigma0 = 2,
                      sigma = 1,
                      TRACE_TIMING = TRUE,
                      eps = 1e-10,
                      rel.tol = NULL) {

  stopifnot(inherits(data, "fbf_data"))

  fam_ids <- data$settings$fam_ids
  pc <- data$settings$pc
  if (is.null(rel.tol)) rel.tol <- data$settings$rel.tol

  var_ids_all <- data$var_ids_all
  F_vec <- as.numeric(data$F_vec_global[var_ids_all]) #Functional status Fv
  H_vec <- as.numeric(data$H_vec_global[var_ids_all]) # unusualness score Hv
  if (!use_H) H_vec <- NULL

  names(F_vec) <- var_ids_all
  if (!is.null(H_vec)) names(H_vec) <- var_ids_all

  pieces_mendel_global <- data$pieces_mendel_global # Mendelian integral pieces pG_half and p0_half
  global_var_info <- data$global_var_info # family-specific carrier configurations for each variant
  zero_procPed_by_fam <- data$zero_procPed_by_fam # zero-configuration pedigree objects

  config_key_vec <- data$pool$config_key_vec
  rep_vid_by_config <- data$pool$rep_vid_by_config # one representative variant per repeated global configuration

  #Set the prior center for baseline over-transmission
  nr <- length(var_ids_all) # number of rare variants
  p0_start <- min(max(1 / nr, eps), 1 - eps)

  if (is.null(eta_start)) {
    eta_start <- if (use_H) {
      c(stats::qlogis(p0_start), 0.1, 0.1) # initial value of eta (logit(1/nr),0.1,0.1)
    } else {
      c(stats::qlogis(p0_start), 0.1)
    }
  }

  .now <- function() as.numeric(Sys.time())

  timings <- data.frame(
    iter = integer(),
    log_lambda = double(),
    lambda = double(),
    secs = double(),
    eta0 = double(),
    eta1 = double(),
    eta2 = double(),
    loglik = double()
  )

  # Define global observed-configuration probability: P(Gv=gv|tau=prod P(Cvi|tau))
  # Each family probability is computed by computeSharingProb_recursive
  # For every tau, the function calculates the probability that the observed carrier pattern arose under the one-founder-entry assumption with founder carrier probability pc
  global_config_prob <- function(vid, tau_vec) {
    tau_vec <- as.numeric(tau_vec)
    info <- global_var_info[[vid]]

    out <- rep(1, length(tau_vec))
    for (k in seq_along(fam_ids)) {
      out <- out * computeSharingProb_recursive(info$procPed_by_fam[[k]], tau_vec, pc)
    }
    out
  }

  # Global zero-configuration probability: P(Gv=0|tau)=prod P(Cvi=0|tau)
  # The probability that no sequenced affected individual in any family carries the variant
  global_zero_prob <- function(tau_vec) {
    tau_vec <- as.numeric(tau_vec)

    out <- rep(1, length(tau_vec))
    for (fid in fam_ids) {
      out <- out * computeSharingProb_recursive(zero_procPed_by_fam[[fid]], tau_vec, pc)
    }
    out
  }

  # Numerical Integration for Mg(lambda)=int(0,1) P(Gv=gv|tau(y,\lambda))
  int_pair_one_variant_global <- function(vid, lambda, rel.tol) {
    vG <- function(y) {
      tau <- taufunc(y, lambda)
      global_config_prob(vid, tau)
    }

    intG <- stats::integrate(vG, 0, 1, rel.tol = rel.tol)$value
    c(intG = intG)
  }

  #Numerical Integration for M0(lambda)= int(0,1) P(Gv=0|tau(y,lambda))
  int_zero_global <- function(lambda, rel.tol) {
    v0 <- function(y) {
      tau <- taufunc(y, lambda)
      global_zero_prob(tau)
    }

    stats::integrate(v0, 0, 1, rel.tol = rel.tol)$value
  }

  # Build all integrals for a fixed lambda
  # For each lambda this build intG(MG(lambda)) and int0(M0(lambda)) for every variant
  build_int_list_at_lambda_global_pooled <- function(lambda) {
    int0 <- int_zero_global(lambda, rel.tol)

    int_by_config <- lapply(rep_vid_by_config, function(rep_vid) {
      c(int_pair_one_variant_global(rep_vid, lambda, rel.tol), int0 = int0)
    })

    out <- vector("list", length(var_ids_all))
    names(out) <- var_ids_all

    for (vid in var_ids_all) {
      key <- config_key_vec[vid]
      out[[vid]] <- int_by_config[[key]]
    }

    out
  }

  # Inner log-likelihood over eta
  #DG=(1-pi)*pG_half + pi*intG
  # D0 = (1-pi)*pG_half + pi*int0
  loglik_eta <- function(eta, int_list) {
    if (any(eta[-1] < 0)) return(-Inf) # eta1>=0, eta2>=0

    X_eta <- eta[1] + eta[2] * F_vec
    if (use_H) X_eta <- X_eta + eta[3] * H_vec # pi_v = expit(eta)+eta1*F_v+eta2*H_v
    pi_vec <- stats::plogis(X_eta)

    pG_half <- vapply(var_ids_all, function(vid) pieces_mendel_global[[vid]]$pG_half, numeric(1))
    p0_half <- vapply(var_ids_all, function(vid) pieces_mendel_global[[vid]]$p0_half, numeric(1))
    intG <- vapply(var_ids_all, function(vid) as.numeric(int_list[[vid]]["intG"]), numeric(1))
    int0 <- vapply(var_ids_all, function(vid) as.numeric(int_list[[vid]]["int0"]), numeric(1))

    DG <- (1 - pi_vec) * pG_half + pi_vec * intG
    D0 <- (1 - pi_vec) * p0_half + pi_vec * int0

    # The log composite likelihood: l(eta,lambda)=sum(logDG - log(1-D0))
    ll <- sum(log(pmax(DG, eps)) - log(pmax(1 - D0, eps)))

    # Penalty term is added: l_pen = l(eta,lambda) + logp(eta)
    ll + penalty_eta(eta, nr = nr, sigma0 = sigma0, sigma = sigma, eps = eps)
  }

  grad_eta <- function(eta, int_list) {
    if (any(eta[-1] < 0)) return(rep(NA_real_, length(eta)))

    X_eta <- eta[1] + eta[2] * F_vec
    if (use_H) X_eta <- X_eta + eta[3] * H_vec
    pi_vec <- stats::plogis(X_eta)

    pG_half <- vapply(var_ids_all, function(vid) pieces_mendel_global[[vid]]$pG_half, numeric(1))
    p0_half <- vapply(var_ids_all, function(vid) pieces_mendel_global[[vid]]$p0_half, numeric(1))
    intG <- vapply(var_ids_all, function(vid) as.numeric(int_list[[vid]]["intG"]), numeric(1))
    int0 <- vapply(var_ids_all, function(vid) as.numeric(int_list[[vid]]["int0"]), numeric(1))

    DG <- (1 - pi_vec) * pG_half + pi_vec * intG
    D0 <- (1 - pi_vec) * p0_half + pi_vec * int0

    w <- pi_vec * (1 - pi_vec) *
      ((intG - pG_half) / pmax(DG, eps) +
         (int0 - p0_half) / pmax(1 - D0, eps))

    Xmat <- if (use_H) {
      cbind(1, F_vec, H_vec)
    } else {
      cbind(1, F_vec)
    }

    as.numeric(colSums(Xmat * w)) +
      grad_penalty_eta(eta, nr = nr, sigma0 = sigma0, sigma = sigma, eps = eps)
  }

  maximize_eta <- function(int_list) {
    lower <- c(-Inf, rep(0, length(eta_start) - 1))
    upper <- rep(Inf, length(eta_start))

    opt <- stats::optim(
      par = eta_start,
      fn = function(e) {
        val <- loglik_eta(e, int_list)
        if (!is.finite(val)) return(1e100)
        -val
      },
      gr = function(e) {
        g <- grad_eta(e, int_list)
        if (any(!is.finite(g))) return(rep(0, length(e)))
        -g
      },
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 1000)
    )

    list(
      eta_hat = opt$par,
      loglik = -opt$value,
      convergence = opt$convergence,
      message = opt$message
    )
  }

  profile_ll_outer <- local({
    iter <- 0L
    function(log_lambda) {
      lambda <- exp(log_lambda)
      t0 <- .now()

      int_list_global <- build_int_list_at_lambda_global_pooled(lambda)
      fit_eta <- maximize_eta(int_list_global)
      prof_loglik <- fit_eta$loglik

      iter <<- iter + 1L
      secs <- .now() - t0

      eta_tmp <- fit_eta$eta_hat
      eta2_tmp <- if (length(eta_tmp) >= 3) eta_tmp[3] else NA_real_

      if (TRACE_TIMING) {
        timings[nrow(timings) + 1L, ] <<- data.frame(
          iter = iter,
          log_lambda = log_lambda,
          lambda = lambda,
          secs = secs,
          eta0 = eta_tmp[1],
          eta1 = eta_tmp[2],
          eta2 = eta2_tmp,
          loglik = prof_loglik
        )
      }

      list(
        neg_prof = -prof_loglik,
        eta_hat = fit_eta$eta_hat,
        lambda = lambda,
        int_list = int_list_global,
        inner = fit_eta
      )
    }
  })

  outer_obj <- function(log_lambda) profile_ll_outer(log_lambda)$neg_prof

  t0 <- Sys.time()
  opt_outer <- stats::optimize(
    outer_obj,
    interval = log(lambda_interval),
    tol = tol_loglambda
  )

  lambda_hat <- exp(opt_outer$minimum)
  prof_at_hat <- profile_ll_outer(log(lambda_hat))
  t1 <- Sys.time()

  eta_hat <- prof_at_hat$eta_hat
  names(eta_hat) <- if (use_H) c("eta0", "eta1_F", "eta2_H") else c("eta0", "eta1_F")

  pi00_hat <- stats::plogis(eta_hat[1])
  pi10_hat <- stats::plogis(eta_hat[1] + eta_hat[2])
  pi01_hat <- if (use_H) stats::plogis(eta_hat[1] + eta_hat[3]) else NA_real_
  pi11_hat <- if (use_H) stats::plogis(sum(eta_hat)) else NA_real_

  fit <- list(
    call = match.call(),
    settings = list(
      lambda_interval = lambda_interval,
      tol_loglambda = tol_loglambda,
      rel.tol = rel.tol,
      eps = eps,
      sigma0 = sigma0,
      sigma = sigma,
      use_H = use_H,
      TRACE_TIMING = TRACE_TIMING
    ),
    lambda_hat = lambda_hat,
    eta_hat = eta_hat,
    pi00_hat = pi00_hat,
    pi10_hat = pi10_hat,
    pi01_hat = pi01_hat,
    pi11_hat = pi11_hat,
    loglik_prof = -prof_at_hat$neg_prof,
    int_list_hat = prof_at_hat$int_list,
    outer = opt_outer,
    inner_at_hat = prof_at_hat$inner,
    timings = timings,
    total_elapsed_sec = as.numeric(difftime(t1, t0, units = "secs"))
  )

  class(fit) <- "fbf_fit"
  fit
}
