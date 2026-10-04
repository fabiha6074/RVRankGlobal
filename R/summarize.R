
#' Summarize rare variant posterior ranking results
#'
#' @param fit Output from [fit_model()].
#' @param post Output from [variant_posteriors()].
#' @param data Output from [make_data()].
#' @param chr_keep Optional vector of chromosomes to keep. If NULL, all chromosomes are included.
#' @param rank_by Metric used to rank variants.
#' @param metrics Metrics to include in plots.
#' @param top_n Number of top variants to return.
#' @param facet_by_chr Logical. If TRUE, facet plots by chromosome.
#' @param plot Logical. If TRUE, create plots.
#' @param base_size Base font size for ggplot objects.
#'
#' @return A list with ranked results, top variants, causal variant summaries, plots, and settings.
#' @export
#'
# summarize and rank the rare variants

summarize_rv <- function(fit,
                          post,
                          data,
                          chr_keep = NULL,
                          rank_by = c("posterior_O1", "EB", "EB_given_O1", "pi_hat", "H"),
                          metrics = c("posterior_O1", "EB"),
                          top_n = 50,
                          facet_by_chr = TRUE,
                          plot = TRUE,
                          base_size = 12) {

  stopifnot(inherits(data, "fbf_data"))
  stopifnot(is.data.frame(post))

  chr_col <- "chr"
  pos_col <- "pos"
  flag_col <- "causal_flag"
  pos_unit <- 1e6

  rank_by <- match.arg(rank_by)

  if (is.null(data$meta_all)) stop("data$meta_all is NULL.")
  meta_all <- data$meta_all

  if (!(chr_col %in% names(meta_all))) stop("meta_all missing chr_col: ", chr_col)
  if (!(pos_col %in% names(meta_all))) stop("meta_all missing pos_col: ", pos_col)

  need_post <- c("variant_id", "posterior_O1", "EB", "EB_given_O1", "F", "H", "pi_hat")
  miss_post <- setdiff(need_post, names(post))
  if (length(miss_post) > 0) stop("post missing columns: ", paste(miss_post, collapse = ", "))

  variant_family_label <- function(vid) {
    vid <- as.character(vid)

    if (is.null(data$global_var_info[[vid]])) {
      return("Unknown")
    }

    info <- data$global_var_info[[vid]]

    fam_has_carrier <- vapply(seq_along(data$settings$fam_ids), function(k) {
      length(info$carriers_by_fam[[k]]) > 0
    }, logical(1))

    fams <- data$settings$fam_ids[fam_has_carrier]

    if (length(fams) == 0) return("No carrier")

    paste(fams, collapse = "+")
  }

  # Avoid duplicate chr/pos columns if variant_posteriors() already merged meta.
  if (all(c(chr_col, pos_col) %in% names(post))) {
    results0 <- tibble::as_tibble(post)
  } else {
    results0 <- tibble::as_tibble(post) |>
      dplyr::left_join(meta_all, by = "variant_id")
  }

  results <- results0 |>
    dplyr::filter(!is.na(.data[[chr_col]]), !is.na(.data[[pos_col]])) |>
    dplyr::mutate(
      chr_plot = as.character(.data[[chr_col]]),
      pos_num = suppressWarnings(as.numeric(.data[[pos_col]])),
      pos_plot = .data$pos_num / pos_unit,
      fam_shape = vapply(.data$variant_id, variant_family_label, character(1))
    ) |>
    dplyr::filter(is.finite(.data$pos_num))

  if (!is.null(chr_keep)) {
    results <- results |>
      dplyr::filter(.data$chr_plot %in% as.character(chr_keep))
  }

  ranked <- results |>
    dplyr::arrange(dplyr::desc(.data[[rank_by]])) |>
    dplyr::mutate(rank = dplyr::row_number())

  top_tbl <- ranked |>
    dplyr::slice_head(n = min(top_n, nrow(ranked)))

  fam_txt <- if (!is.null(data$settings$fam_ids)) paste(data$settings$fam_ids, collapse = ",") else "?"

  chr_vals <- unique(ranked$chr_plot)
  chr_vals_num <- suppressWarnings(as.numeric(chr_vals))

  if (all(is.finite(chr_vals_num))) {
    chr_vals <- chr_vals[order(chr_vals_num)]
  } else {
    chr_vals <- sort(chr_vals)
  }

  chr_txt <- paste(chr_vals, collapse = ",")

  lam_txt <- if (!is.null(fit$lambda_hat)) sprintf("lambda_hat=%.3g", fit$lambda_hat) else "lambda_hat=?"
  eta_txt <- if (!is.null(fit$eta_hat)) {
    paste(sprintf("%s=%.2f", names(fit$eta_hat), fit$eta_hat), collapse = ", ")
  } else {
    "eta_hat=?"
  }

  has_flag <- !is.null(flag_col) && flag_col %in% names(ranked)

  make_metric_plot <- function(df, metric) {
    ylab <- switch(metric,
                   posterior_O1 = "Posterior P(O=1 | G,F,H)",
                   EB = "Posterior mean E[B | G,F,H]",
                   EB_given_O1 = "Posterior mean E[B | O=1,G,F,H]",
                   pi_hat = "Prior probability pi_hat(F,H)",
                   H = "Age-onset unusualness score H",
                   metric)

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$pos_plot,
        y = .data[[metric]],
        shape = .data$fam_shape
      )
    ) +
      ggplot2::geom_point(color = "grey70", size = 1.8)

    if (has_flag) {
      p <- p +
        ggplot2::geom_point(
          data = dplyr::filter(df, .data[[flag_col]] %in% c(TRUE, 1, "1")),
          color = "red",
          size = 2.6
        ) +
        ggplot2::geom_vline(
          data = dplyr::filter(df, .data[[flag_col]] %in% c(TRUE, 1, "1")),
          ggplot2::aes(xintercept = .data$pos_plot),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "orange",
          show.legend = FALSE
        )
    }

    if (facet_by_chr) {
      p <- p + ggplot2::facet_wrap(~ chr_plot, scales = "free_x")
    }

    p +
      ggplot2::labs(
        x = sprintf("Genomic position (%s)", if (pos_unit == 1e6) "Mb" else "bp"),
        y = ylab,
        shape = "Carrier family",
        title = sprintf("chr=%s | fam=%s | %s | %s", chr_txt, fam_txt, metric, lam_txt),
        subtitle = eta_txt
      ) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(legend.position = "bottom")
  }

  plots <- list()
  if (plot) {
    for (m in metrics) {
      if (!(m %in% names(ranked))) stop("Requested metric not in results/post: ", m)
      plots[[m]] <- make_metric_plot(ranked, m)
    }
  }

  list(
    results = results,
    ranked = ranked,
    top = top_tbl,
    plots = plots,
    settings = list(
      chr_keep = chr_keep,
      rank_by = rank_by,
      metrics = metrics,
      facet_by_chr = facet_by_chr,
      top_n = top_n
    )
  )
}
