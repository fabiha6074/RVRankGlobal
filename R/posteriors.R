# transform b to y

.b_from_y <- function(y, lambda) {
  y <- pmin(pmax(y, 0), 1 - 1e-12)
  -log1p(-y) / lambda
}

# compute posterior measures

#' It takes the fitted model and computes variant level ranking quantities
#' For each variant, it calculates posterior probability of over-transmission and posterior mean of random transmission effect
#'
#' @param fit Output from fit_model()
#' @param data Output from make_data()
#' @param pc Founder carrier probability
#' @param rel.tol Relative tolerance for numerical integration
#' @param eps umerical guard for likelihood/ initial probabilities.
#' @param pool_by_config keep record of variants with same configuration
#' @return An object of class `"post"`
#' @export

variant_posteriors <- function(fit,
                               data,
                               pc = NULL,
                               rel.tol = NULL,
                               eps = 1e-12,
                               pool_by_config = TRUE) {

  stopifnot(inherits(data, "fbf_data"))
  stopifnot(inherits(fit, "fbf_fit") || is.list(fit))

  fam_ids <- data$settings$fam_ids
  if (is.null(pc)) pc <- data$settings$pc
  if (is.null(rel.tol)) rel.tol <- data$settings$rel.tol

  var_ids_all <- data$var_ids_all
  F_vec <- as.numeric(data$F_vec_global[var_ids_all])
  H_vec <- as.numeric(data$H_vec_global[var_ids_all])

  pieces_mendel_global <- data$pieces_mendel_global
  global_var_info <- data$global_var_info

  lambda_hat <- fit$lambda_hat
  eta_hat <- fit$eta_hat
  int_list_hat <- fit$int_list_hat

  if (is.null(lambda_hat) || is.null(eta_hat) || is.null(int_list_hat)) {
    stop("fit must contain lambda_hat, eta_hat, and int_list_hat.")
  }

  use_H <- isTRUE(fit$settings$use_H) && length(eta_hat) >= 3

  linpred <- eta_hat[1] + eta_hat[2] * F_vec
  if (use_H) linpred <- linpred + eta_hat[3] * H_vec

  pi_by_var <- stats::plogis(linpred)
  names(pi_by_var) <- var_ids_all

  posterior_O1_global <- function(pi_i, mendel_piece, int_piece, eps = eps) {
    pG_half <- mendel_piece$pG_half
    intG <- as.numeric(int_piece["intG"])

    DG <- (1 - pi_i) * pG_half + pi_i * intG
    pi_i * intG / pmax(DG, eps)
  }

  post_O1_vec <- mapply(
    function(vid, pi_i) {
      posterior_O1_global(
        pi_i = pi_i,
        mendel_piece = pieces_mendel_global[[vid]],
        int_piece = int_list_hat[[vid]],
        eps = eps
      )
    },
    vid = var_ids_all,
    pi_i = pi_by_var,
    SIMPLIFY = TRUE
  )
  names(post_O1_vec) <- var_ids_all

  global_config_prob <- function(vid, tau_vec) {
    tau_vec <- as.numeric(tau_vec)
    info <- global_var_info[[vid]]

    out <- rep(1, length(tau_vec))
    for (k in seq_along(fam_ids)) {
      out <- out * computeSharingProb_recursive(info$procPed_by_fam[[k]], tau_vec, pc)
    }
    out
  }

  intB_one_variant_global <- function(vid, lambda, rel.tol) {
    num_fun <- function(y) {
      tau <- taufunc(y, lambda)
      b <- .b_from_y(y, lambda)
      b * global_config_prob(vid, tau)
    }

    stats::integrate(num_fun, 0, 1, rel.tol = rel.tol)$value
  }

  if (pool_by_config) {
    config_key_vec <- data$pool$config_key_vec
    rep_vid_by_config <- data$pool$rep_vid_by_config

    intB_by_config <- vapply(rep_vid_by_config, function(rep_vid) {
      intB_one_variant_global(rep_vid, lambda_hat, rel.tol)
    }, numeric(1))

    intB_vec <- vapply(var_ids_all, function(vid) {
      key <- config_key_vec[[vid]]
      intB_by_config[[key]]
    }, numeric(1))
    names(intB_vec) <- var_ids_all
  } else {
    intB_vec <- vapply(
      var_ids_all,
      function(vid) intB_one_variant_global(vid, lambda_hat, rel.tol),
      numeric(1)
    )
    names(intB_vec) <- var_ids_all
  }

  intG_vec <- vapply(var_ids_all, function(vid) as.numeric(int_list_hat[[vid]]["intG"]), numeric(1))
  names(intG_vec) <- var_ids_all

  EB_given_O1_vec <- intB_vec / pmax(intG_vec, eps)
  names(EB_given_O1_vec) <- var_ids_all

  # E[B | G,F,H] = pi_v * intB_v / D_Gv
  EB_vec <- vapply(var_ids_all, function(vid) {
    pi_i <- pi_by_var[[vid]]
    pG_half <- pieces_mendel_global[[vid]]$pG_half
    intG <- as.numeric(int_list_hat[[vid]]["intG"])
    DG <- (1 - pi_i) * pG_half + pi_i * intG
    pi_i * intB_vec[[vid]] / pmax(DG, eps)
  }, numeric(1))
  names(EB_vec) <- var_ids_all

  out <- data.frame(
    variant_id = var_ids_all,
    class = as.integer(data$class_vec_global[var_ids_all]),
    F = as.integer(F_vec),
    H_raw = as.numeric(data$H_raw_global[var_ids_all]),
    H = as.numeric(H_vec),
    pi_hat = as.numeric(pi_by_var[var_ids_all]),
    posterior_O1 = as.numeric(post_O1_vec[var_ids_all]),
    EB_given_O1 = as.numeric(EB_given_O1_vec[var_ids_all]),
    EB = as.numeric(EB_vec[var_ids_all]),
    stringsAsFactors = FALSE
  )

  out$pG_half <- vapply(var_ids_all, function(vid) pieces_mendel_global[[vid]]$pG_half, numeric(1))
  out$p0_half <- vapply(var_ids_all, function(vid) pieces_mendel_global[[vid]]$p0_half, numeric(1))
  out$intG <- as.numeric(intG_vec[var_ids_all])
  out$int0 <- vapply(var_ids_all, function(vid) as.numeric(int_list_hat[[vid]]["int0"]), numeric(1))
  out$intB <- as.numeric(intB_vec[var_ids_all])

  if (!is.null(data$meta_all) && "variant_id" %in% names(data$meta_all)) {
    out <- merge(data$meta_all, out, by = "variant_id", all.y = TRUE, sort = FALSE)
  }

  out[order(-out$posterior_O1), , drop = FALSE]
}
