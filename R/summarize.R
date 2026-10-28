#' Summarize rare variant posterior ranking results
#'
#' @param fit Output from [fit_model()].
#' @param post Output from [variant_posteriors()].
#' @param data Output from [make_data()].
#' @param chr_keep Optional vector of chromosomes to keep. If `NULL`, all
#'   chromosomes are included.
#' @param rank_by Metric used to rank variants.
#' @param metrics Metrics to plot.
#' @param top_n Number of top variants to return.
#' @param plot Logical. If `TRUE`, draw the plots immediately.
#' @param show_estimates Logical. If `TRUE`, show fitted parameter estimates
#'   in mathematical notation below the title.
#' @param position_unit Unit used for genomic-position labels on the x-axis.
#'   Use `"Mb"` for megabases or `"bp"` for base pairs.
#' @param point_cex Point-size multiplier.
#' @param main_cex Title-size multiplier.
#'
#' @return A list containing ranked results, top variants, reusable base-R
#'   plotting functions, and settings.
#' @export
summarize_rv <- function(fit,
                         post,
                         data,
                         chr_keep = NULL,
                         rank_by = c("posterior_O1", "EB", "EB_given_O1", "pi_hat", "H"),
                         metrics = c("posterior_O1", "EB"),
                         top_n = 50,
                         plot = TRUE,
                         show_estimates = TRUE,
                         position_unit = c("Mb", "bp"),
                         point_cex = 0.9,
                         main_cex = 1.05) {

  if (!inherits(data, "fbf_data")) {
    stop("data must be an object returned by make_data().")
  }
  if (!(inherits(fit, "fbf_fit") || is.list(fit))) {
    stop("fit must be an object returned by fit_model().")
  }
  if (!is.data.frame(post)) {
    stop("post must be a data frame returned by variant_posteriors().")
  }

  rank_by <- match.arg(rank_by)
  position_unit <- match.arg(position_unit)

  required_post <- c(
    "variant_id", "posterior_O1", "EB", "EB_given_O1", "F", "H", "pi_hat"
  )
  missing_post <- setdiff(required_post, names(post))
  if (length(missing_post) > 0L) {
    stop("post missing columns: ", paste(missing_post, collapse = ", "))
  }

  if (is.null(data$meta_all)) {
    stop("data$meta_all is NULL.")
  }

  meta_all <- as.data.frame(data$meta_all)
  if (!("variant_id" %in% names(meta_all))) {
    stop("data$meta_all is missing variant_id.")
  }
  if (!("chr" %in% names(meta_all))) {
    stop("data$meta_all is missing chr.")
  }
  if (!("pos" %in% names(meta_all))) {
    stop("data$meta_all is missing pos.")
  }

  # Add missing metadata without allowing merge() to reorder variants.
  results <- as.data.frame(post)
  meta_match <- match(results$variant_id, meta_all$variant_id)

  metadata_columns <- intersect(
    c("chr", "pos", "Type", "causal_flag"),
    names(meta_all)
  )

  for (column_name in metadata_columns) {
    if (!(column_name %in% names(results))) {
      results[[column_name]] <- meta_all[[column_name]][meta_match]
    }
  }

  if (!("chr" %in% names(results)) || !("pos" %in% names(results))) {
    stop("Chromosome or position metadata could not be added to post.")
  }

  variant_family_label <- function(variant_id) {
    variant_id <- as.character(variant_id)
    info <- data$global_var_info[[variant_id]]

    if (is.null(info)) return("Unknown")

    has_carrier <- vapply(
      seq_along(data$settings$fam_ids),
      function(k) length(info$carriers_by_fam[[k]]) > 0L,
      logical(1)
    )

    families <- data$settings$fam_ids[has_carrier]
    if (length(families) == 0L) return("No carrier")

    paste(families, collapse = "+")
  }

  results$chr_plot <- as.character(results$chr)
  results$pos_num <- suppressWarnings(as.numeric(results$pos))
  results$fam_shape <- vapply(
    results$variant_id,
    variant_family_label,
    character(1)
  )

  keep <- !is.na(results$chr_plot) & is.finite(results$pos_num)
  results <- results[keep, , drop = FALSE]

  if (!is.null(chr_keep)) {
    results <- results[
      results$chr_plot %in% as.character(chr_keep),
      ,
      drop = FALSE
    ]
  }

  if (nrow(results) == 0L) {
    stop("No variants remain after applying chromosome and position filters.")
  }

  # Determine chromosome order: numeric chromosomes first, then X/Y/MT.
  chromosome_names <- unique(results$chr_plot)
  chromosome_numeric <- suppressWarnings(as.numeric(chromosome_names))
  chromosome_order <- order(
    is.na(chromosome_numeric),
    chromosome_numeric,
    chromosome_names
  )
  chromosome_names <- chromosome_names[chromosome_order]

  # Maximum observed position on each chromosome.
  chromosome_max <- vapply(
    chromosome_names,
    function(chromosome) {
      max(results$pos_num[results$chr_plot == chromosome], na.rm = TRUE)
    },
    numeric(1)
  )

  chromosome_offset <- c(0, utils::head(cumsum(chromosome_max), -1))
  names(chromosome_offset) <- chromosome_names

  chromosome_centre <- chromosome_offset + chromosome_max / 2
  chromosome_band <- (seq_along(chromosome_names) - 1L) %% 2L
  names(chromosome_band) <- chromosome_names

  results$pos_manhattan <- results$pos_num +
    unname(chromosome_offset[results$chr_plot])
  results$chr_band <- unname(chromosome_band[results$chr_plot])

  ranked <- results[
    order(results[[rank_by]], decreasing = TRUE, na.last = TRUE),
    ,
    drop = FALSE
  ]
  ranked$rank <- seq_len(nrow(ranked))

  top_count <- min(top_n, nrow(ranked))
  top_tbl <- ranked[seq_len(top_count), , drop = FALSE]

  # Assign one plotting symbol to each carrier family.
  family_names <- sort(unique(ranked$fam_shape))
  available_shapes <- c(16, 17, 15, 18, 8, 3, 7, 4)
  family_shapes <- rep(available_shapes, length.out = length(family_names))
  names(family_shapes) <- family_names

  chromosome_colors <- c("#27408B", "#87A7D9")

  # Build within-chromosome genomic-position ticks, translated to the
  # cumulative Manhattan coordinates used for plotting.
  position_scale <- if (position_unit == "Mb") 1e6 else 1
  position_axis_title <- paste0("Genomic position (", position_unit, ")")
  position_tick_at <- numeric(0)
  position_tick_labels <- character(0)

  for (i in seq_along(chromosome_names)) {
    chromosome <- chromosome_names[i]
    local_ticks <- pretty(c(0, chromosome_max[i]), n = 3)
    local_ticks <- local_ticks[
      local_ticks >= 0 & local_ticks <= chromosome_max[i]
    ]

    # Do not print two zero labels at a boundary between chromosomes.
    if (i > 1L) local_ticks <- local_ticks[local_ticks > 0]

    position_tick_at <- c(
      position_tick_at,
      chromosome_offset[chromosome] + local_ticks
    )
    position_tick_labels <- c(
      position_tick_labels,
      format(
        local_ticks / position_scale,
        trim = TRUE,
        scientific = FALSE,
        digits = 4
      )
    )
  }

  metric_axis_label <- function(metric) {
    switch(
      metric,
      posterior_O1 = "Posterior P(O=1 | G,F,H)",
      EB = "Posterior mean E[B | G,F,H]",
      EB_given_O1 = "Posterior mean E[B | O=1,G,F,H]",
      pi_hat = "Prior probability pi_hat(F,H)",
      H = "Age-onset unusualness score H",
      metric
    )
  }

  metric_title <- function(metric) {
    switch(
      metric,
      posterior_O1 = "Posterior Evidence of Variant Over-Transmission",
      EB = "Estimated Transmission Effects Across the Genome",
      EB_given_O1 = "Transmission Effects Among Over-Transmitted Variants",
      pi_hat = "Estimated Prior Probability of Over-Transmission",
      H = "Age-at-Onset Unusualness Across the Genome",
      "Rare-Variant Prioritization"
    )
  }

  parameter_expression <- function() {
    if (!show_estimates) return(NULL)

    format_estimate <- function(x) {
      trimws(formatC(x, digits = 3, format = "g"))
    }

    expression_parts <- character(0)

    if (!is.null(fit$lambda_hat) && is.finite(fit$lambda_hat)) {
      expression_parts <- c(
        expression_parts,
        paste0("hat(lambda) == ", format_estimate(fit$lambda_hat))
      )
    }

    if (!is.null(fit$eta_hat)) {
      eta_values <- as.numeric(fit$eta_hat)

      for (j in seq_along(eta_values)) {
        if (is.finite(eta_values[j])) {
          expression_parts <- c(
            expression_parts,
            paste0(
              "hat(eta)[", j - 1L, "] == ",
              format_estimate(eta_values[j])
            )
          )
        }
      }
    }

    if (length(expression_parts) == 0L) return(NULL)

    parse(
      text = paste(expression_parts, collapse = ' ~~ "," ~~ ')
    )[[1]]
  }

  causal_indicator <- rep(FALSE, nrow(ranked))
  if ("causal_flag" %in% names(ranked)) {
    causal_indicator <- tolower(as.character(ranked$causal_flag)) %in%
      c("true", "1", "yes", "y", "causal")
  }

  draw_metric <- function(metric, show_family_legend = TRUE) {
    y <- ranked[[metric]]
    finite_y <- is.finite(y)

    if (!any(finite_y)) {
      stop("Metric contains no finite values: ", metric)
    }

    point_colors <- chromosome_colors[ranked$chr_band + 1L]
    point_shapes <- unname(family_shapes[ranked$fam_shape])

    graphics::plot(
      x = ranked$pos_manhattan[finite_y],
      y = y[finite_y],
      type = "n",
      xaxt = "n",
      xlab = "",
      ylab = metric_axis_label(metric),
      main = metric_title(metric),
      cex.main = main_cex,
      bty = "l"
    )

    graphics::axis(
      side = 1,
      at = position_tick_at,
      labels = position_tick_labels,
      cex.axis = 0.7,
      tck = -0.012
    )

    # Keep chromosome context on a second, unobtrusive row.
    graphics::axis(
      side = 1,
      at = chromosome_centre,
      labels = paste0("Chr ", chromosome_names),
      tick = FALSE,
      line = 1.15,
      cex.axis = 0.75
    )
    graphics::mtext(
      position_axis_title,
      side = 1,
      line = 2.75,
      cex = 0.85
    )

    # Light vertical separators between chromosomes.
    if (length(chromosome_offset) > 1L) {
      graphics::abline(
        v = chromosome_offset[-1],
        col = "grey90",
        lty = 3
      )
    }

    graphics::points(
      x = ranked$pos_manhattan[finite_y],
      y = y[finite_y],
      pch = point_shapes[finite_y],
      col = point_colors[finite_y],
      cex = point_cex
    )

    causal_to_plot <- causal_indicator & finite_y
    if (any(causal_to_plot)) {
      graphics::abline(
        v = ranked$pos_manhattan[causal_to_plot],
        col = "orange",
        lty = 2,
        lwd = 1
      )

      graphics::points(
        x = ranked$pos_manhattan[causal_to_plot],
        y = y[causal_to_plot],
        pch = point_shapes[causal_to_plot],
        col = "red",
        cex = point_cex * 1.5,
        lwd = 1.5
      )
    }

    if (show_family_legend) {
      plot_limits <- graphics::par("usr")
      legend_x <- plot_limits[2] + 0.04 * diff(plot_limits[1:2])
      legend_y <- plot_limits[4]

      graphics::legend(
        x = legend_x,
        y = legend_y,
        legend = family_names,
        pch = family_shapes,
        col = "grey20",
        title = "Carrier family",
        bty = "n",
        cex = 0.72,
        ncol = 1,
        xpd = NA,
        xjust = 0,
        yjust = 1,
        x.intersp = 0.7,
        y.intersp = 0.85
      )
    }

    estimate_label <- parameter_expression()
    if (!is.null(estimate_label)) {
      graphics::mtext(
        estimate_label,
        side = 3,
        line = 0.15,
        cex = 0.75
      )
    }

    invisible(NULL)
  }

  make_plotter <- function(metric) {
    force(metric)
    function() {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mar = c(5.4, 4.5, 4.2, 7.2))
      draw_metric(metric, show_family_legend = TRUE)
    }
  }

  plotters <- setNames(
    lapply(metrics, make_plotter),
    metrics
  )

  if (plot && length(metrics) > 0L) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    number_of_plots <- length(metrics)
    number_of_columns <- if (number_of_plots <= 2L) number_of_plots else 2L
    number_of_rows <- ceiling(number_of_plots / number_of_columns)

    graphics::par(
      mfrow = c(number_of_rows, number_of_columns),
      mar = c(5.4, 4.5, 4.2, 1.2),
      oma = c(0, 0, 0, 7.2)
    )

    for (metric_index in seq_along(metrics)) {
      metric <- metrics[metric_index]
      if (!(metric %in% names(ranked))) {
        stop("Requested metric not in results/post: ", metric)
      }
      draw_metric(
        metric,
        show_family_legend = metric_index == length(metrics)
      )
    }
  }

  list(
    results = results,
    ranked = ranked,
    top = top_tbl,
    plots = plotters,
    settings = list(
      chr_keep = chr_keep,
      rank_by = rank_by,
      metrics = metrics,
      show_estimates = show_estimates,
      position_unit = position_unit,
      top_n = top_n
    )
  )
}
