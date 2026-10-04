

#' Build data for the age-at-onset global likelihood model
#'
#' Converts pedigree, genotype, variant characteristics, and age-specific hazard data
#' into the data object required by the penalized composite global likelihood model.
#'
#' @param fam_ids Character or numeric vector of family IDs to include.
#' @param chr_ids Character or numeric vector of chromosome IDs to include.
#' @param ordered_study_peds Pedigree data frame.
#' @param geno_list List of genotype matrices by chromosome.
#' @param variant_list List of variant metadata by chromosome.
#' @param hazards Optional age-specific hazard table.
#' @param h_col Name of individual-level unusualness score column. Used if `hazards` is NULL and the pedigree already contains unusualness scores.
#' @param birth_col Name of birth year column.
#' @param onset_col Name of onset year column.
#' @param hazard_col Name of population onset hazard column.
#' @param pc Founder carrier probability.
#' @param rel.tol Relative tolerance for numerical integration.
#'
#' @return An object of class `"fbf_data"`.
#' @export


make_data <- function(fam_ids,
                      chr_ids,
                      ordered_study_peds,
                      geno_list,
                      variant_list,
                      hazards = NULL,
                      birth_col = "birthYr",
                      onset_col = "onsetYr",
                      hazard_col = "pop_onset_hazard",
                      h_col = "h_indiv",
                      pc = 0.00032,
                      rel.tol = 1e-8) {

  stopifnot(length(fam_ids) > 0, length(chr_ids) > 0)
  fam_ids <- as.character(fam_ids)

  # -----------------------
  # 0) Add individual-level unusualness if hazards are supplied
  # -----------------------
  if (!is.null(hazards)) {
    ordered_study_peds <- add_age_unusualness(
      ped = ordered_study_peds,
      hazards = hazards,
      birth_col = birth_col,
      onset_col = onset_col,
      hazard_col = hazard_col
    )
    h_col <- "h_indiv"
  }

  has_h <- h_col %in% names(ordered_study_peds)

  # -----------------------
  # 1) Assemble geno_data
  # -----------------------
  geno_data <- do.call(cbind, geno_list[chr_ids])
  if ("unique_ID" %in% colnames(geno_data)) {
    geno_data <- geno_data[, setdiff(colnames(geno_data), "unique_ID"), drop = FALSE]
  }
  stopifnot(!is.null(rownames(geno_data)))

  # -----------------------
  # 2) Build meta_all and F map
  # -----------------------
  meta_all <- build_meta_all(chr_ids, variant_list)
  class_map <- setNames(map_class(meta_all$Type), meta_all$variant_id)

  # -----------------------
  # 3) Build long genotype data
  # -----------------------
  long_geno <- geno_data |>
    as.data.frame() |>
    tibble::rownames_to_column("unique_ID") |>
    tidyr::pivot_longer(
      cols = -unique_ID,
      names_to = "variant_id",
      values_to = "copies"
    ) |>
    dplyr::mutate(copies_num = suppressWarnings(as.numeric(copies)))

  # -----------------------
  # 4) Family objects
  # -----------------------
  fam_data <- setNames(vector("list", length(fam_ids)), fam_ids)

  for (fid in fam_ids) {
    ped_onefam <- ordered_study_peds |>
      dplyr::mutate(FamID_chr = as.character(FamID)) |>
      dplyr::filter(FamID_chr == fid)

    aff_tbl <- ped_onefam |>
      dplyr::mutate(affected = affected %in% c(TRUE, 1, "1")) |>
      dplyr::filter(affected) |>
      dplyr::arrange(ID) |>
      dplyr::select(ID, unique_ID, dplyr::any_of(h_col))

    aff_ids_order  <- aff_tbl$ID
    aff_uids_order <- aff_tbl$unique_ID

    join_cols <- c("unique_ID", "FamID_chr", "affected", "ID")
    if (has_h) join_cols <- c(join_cols, h_col)

    carriers_long <- long_geno |>
      dplyr::left_join(
        ped_onefam |>
          dplyr::mutate(affected = affected %in% c(TRUE, 1, "1")) |>
          dplyr::select(dplyr::all_of(join_cols)),
        by = "unique_ID"
      ) |>
      dplyr::filter(FamID_chr == fid, affected, .data$copies_num > 0)

    carriers_list_all <- carriers_long |>
      dplyr::group_by(variant_id) |>
      dplyr::summarise(
        carriers_ID = {
          ids <- aff_ids_order[aff_uids_order %in% unique(unique_ID)]
          paste(ids, collapse = ",")
        },
        H_family = if (has_h) {
          suppressWarnings(max(as.numeric(.data[[h_col]]), na.rm = TRUE))
        } else {
          NA_real_
        },
        .groups = "drop"
      ) |>
      dplyr::filter(.data$carriers_ID != "") |>
      dplyr::mutate(
        carriers_vec = lapply(.data$carriers_ID, parse_ids),
        H_family = dplyr::if_else(is.infinite(.data$H_family), NA_real_, .data$H_family)
      )

    cl <- carriers_list_all$carriers_vec
    names(cl) <- carriers_list_all$variant_id

    H_family_map <- carriers_list_all$H_family
    names(H_family_map) <- carriers_list_all$variant_id

    fam_data[[fid]] <- list(
      ped = ped_onefam,
      carriers_list_all = cl,
      H_family_map = H_family_map,
      aff_ids_order = aff_ids_order,
      aff_uids_order = aff_uids_order
    )
  }

  ped_by_fam <- lapply(fam_ids, function(fid) fam_data[[fid]]$ped)
  names(ped_by_fam) <- fam_ids

  zero_procPed_by_fam <- lapply(fam_ids, function(fid) {
    processPed(ped_by_fam[[fid]], integer(0))
  })
  names(zero_procPed_by_fam) <- fam_ids

  # -----------------------
  # 5) Global variant set and F_v
  # -----------------------
  var_ids_all <- sort(unique(unlist(lapply(fam_data, function(f) names(f$carriers_list_all)))))

  class_vec_global <- class_map[var_ids_all]
  keep <- !is.na(class_vec_global)
  var_ids_all <- var_ids_all[keep]
  class_vec_global <- setNames(as.integer(class_vec_global[keep]), var_ids_all)

  F_vec_global <- class_vec_global
  names(F_vec_global) <- var_ids_all

  # -----------------------
  # 6) Global variant info and raw H*_v
  # -----------------------
  global_var_info <- lapply(var_ids_all, function(vid) {
    carriers_by_fam <- lapply(fam_ids, function(fid) {
      cl_fam <- fam_data[[fid]]$carriers_list_all

      if (length(cl_fam) == 0L || !(vid %in% names(cl_fam))) {
        integer(0)
      } else {
        cl_fam[[vid]]
      }
    })
    names(carriers_by_fam) <- fam_ids

    I_index <- which(vapply(carriers_by_fam, length, integer(1)) > 0L)

    H_by_fam <- vapply(fam_ids, function(fid) {
      hmap_fam <- fam_data[[fid]]$H_family_map

      if (length(hmap_fam) == 0L || !(vid %in% names(hmap_fam))) {
        NA_real_
      } else {
        as.numeric(hmap_fam[vid])
      }
    }, numeric(1))
    names(H_by_fam) <- fam_ids

    H_star <- if (length(I_index) > 0L && any(is.finite(H_by_fam[I_index]))) {
      mean(H_by_fam[I_index], na.rm = TRUE)
    } else {
      NA_real_
    }

    procPed_by_fam <- lapply(seq_along(fam_ids), function(k) {
      processPed(ped_by_fam[[fam_ids[k]]], carriers_by_fam[[k]])
    })
    names(procPed_by_fam) <- fam_ids

    list(
      variant_id = vid,
      carriers_by_fam = carriers_by_fam,
      I_index = I_index,
      H_by_fam = H_by_fam,
      H_star = H_star,
      procPed_by_fam = procPed_by_fam
    )
  })
  names(global_var_info) <- var_ids_all

  H_raw_global <- vapply(global_var_info, function(x) x$H_star, numeric(1))
  names(H_raw_global) <- var_ids_all

  # If no usable age/hazard information is available, use H_v=0 for everyone.
  if (any(is.finite(H_raw_global))) {
    med_H <- stats::median(H_raw_global[is.finite(H_raw_global)], na.rm = TRUE)
    H_raw_global[!is.finite(H_raw_global)] <- med_H

    H_min <- min(H_raw_global, na.rm = TRUE)
    H_max <- max(H_raw_global, na.rm = TRUE)

    H_vec_global <- if (H_max > H_min) {
      (H_raw_global - H_min) / (H_max - H_min)
    } else {
      rep(0, length(H_raw_global))
    }
  } else {
    H_raw_global <- stats::setNames(rep(NA_real_, length(var_ids_all)), var_ids_all)
    H_vec_global <- stats::setNames(rep(0, length(var_ids_all)), var_ids_all)
  }
  names(H_vec_global) <- var_ids_all

  # -----------------------
  # 7) Pooling keys for G-specific integrals
  # -----------------------
  config_key_vec <- vapply(var_ids_all, function(vid) {
    info <- global_var_info[[vid]]
    fam_strings <- vapply(info$carriers_by_fam, encode_family_config, character(1))
    paste(fam_strings, collapse = "|")
  }, character(1))
  names(config_key_vec) <- var_ids_all

  unique_configs <- unique(config_key_vec)
  rep_vid_by_config <- tapply(names(config_key_vec), config_key_vec, `[`, 1)

  # -----------------------
  # 8) Mendelian pieces at tau=1/2
  # -----------------------
  p0_half_f <- vapply(fam_ids, function(fid) {
    as.numeric(computeSharingProb_recursive(zero_procPed_by_fam[[fid]], 0.5, pc))
  }, numeric(1))

  p0_half_global <- prod(p0_half_f)

  pieces_mendel_global <- lapply(var_ids_all, function(vid) {
    info <- global_var_info[[vid]]

    pG_half_f <- vapply(seq_along(fam_ids), function(k) {
      as.numeric(computeSharingProb_recursive(info$procPed_by_fam[[k]], 0.5, pc))
    }, numeric(1))

    pG_half_i <- prod(pG_half_f)

    list(
      pG_half = pG_half_i,
      p0_half = p0_half_global,
      pG_half_f = pG_half_f,
      p0_half_f = p0_half_f,
      I_index = info$I_index
    )
  })
  names(pieces_mendel_global) <- var_ids_all

  flt <- filter_invalid_by_mendel(pieces_mendel_global, eps = 0)
  dropped_by_mendel <- flt$drop_vids

  out <- list(
    settings = list(
      fam_ids = fam_ids,
      chr_ids = chr_ids,
      pc = pc,
      rel.tol = rel.tol,
      birth_col = birth_col,
      onset_col = onset_col,
      hazard_col = hazard_col,
      h_col = h_col
    ),
    geno_data = geno_data,
    meta_all = meta_all,
    class_map = class_map,
    fam_data = fam_data,
    ped_by_fam = ped_by_fam,
    zero_procPed_by_fam = zero_procPed_by_fam,
    var_ids_all = var_ids_all,
    class_vec_global = class_vec_global,
    F_vec_global = F_vec_global,
    H_raw_global = H_raw_global,
    H_vec_global = H_vec_global,
    global_var_info = global_var_info,
    pool = list(
      config_key_vec = config_key_vec,
      unique_configs = unique_configs,
      rep_vid_by_config = rep_vid_by_config
    ),
    pieces_mendel_global = pieces_mendel_global,
    dropped_by_mendel = dropped_by_mendel,
    deps = list(
      computeSharingProb_recursive = computeSharingProb_recursive
    )
  )

  out <- prune_to_vids(out, flt$keep_vids)

  class(out) <- "fbf_data"
  out
}
