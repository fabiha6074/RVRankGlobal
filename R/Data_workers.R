
#add_age_unusualness()
#build_meta_all()
#map_class()
#parse_ids()
#encode_family_config()
#filter_invalid_by_mendel()
#prune_to_vids()

add_age_unusualness <- function(ped,
                                hazards,
                                birth_col = "birthYr",
                                onset_col = "onsetYr",
                                hazard_col = "pop_onset_hazard",
                                eps = 1e-6) {

  # If hazards has no age column, use row number as age: 1, ..., 100
  if (!"age" %in% names(hazards)) {
    hazards <- hazards |>
      tibble::rownames_to_column("age")
  }

  hazard_age <- hazards |>
    dplyr::mutate(
      age = as.integer(age),
      hazard = as.numeric(.data[[hazard_col]]),
      hazard = pmax(hazard, eps),
      h_indiv = -log(hazard)
    ) |>
    dplyr::select(age, h_indiv)

  ped |>
    dplyr::mutate(
      age_onset = as.numeric(.data[[onset_col]]) - as.numeric(.data[[birth_col]]),
      age_lookup = round(age_onset),
      age_lookup = dplyr::if_else(
        is.na(age_lookup),
        NA_integer_,
        pmin(pmax(as.integer(age_lookup), 1L), 100L)
      )
    ) |>
    dplyr::left_join(
      hazard_age,
      by = c("age_lookup" = "age")
    ) |>
    dplyr::select(-age_lookup)

}


build_meta_all <- function(chr_ids, variant_list) {
  pos_candidates    <- c("pos", "Pos", "POS", "Position", "position", "bp")
  causal_candidates <- c("C_SNV", "causal", "is_causal", "CAUSAL")

  out <- list()

  for (chr in chr_ids) {
    chr_str  <- as.character(chr)
    meta_chr <- variant_list[[chr_str]]
    if (is.null(meta_chr)) meta_chr <- variant_list[[chr]]
    if (is.null(meta_chr)) next

    stopifnot("variant_id" %in% names(meta_chr),
              "Type"       %in% names(meta_chr))

    pos_col    <- intersect(pos_candidates,    names(meta_chr))[1]
    causal_col <- intersect(causal_candidates, names(meta_chr))[1]

    meta_chr_norm <- meta_chr |>
      dplyr::mutate(
        pos = if (!is.na(pos_col)) {
          as.numeric(.data[[pos_col]])
        } else {
          suppressWarnings(as.numeric(sub(".*_", "", variant_id)))
        },
        causal_flag = if (!is.na(causal_col)) {
          x <- .data[[causal_col]]
          dplyr::coalesce(
            as.logical(x),
            tolower(as.character(x)) %in% c("true", "t", "yes", "y", "causal"),
            suppressWarnings(as.integer(as.character(x)) == 1),
            FALSE
          )
        } else FALSE
      ) |>
      dplyr::select(variant_id, Type, pos, causal_flag) |>
      dplyr::mutate(chr = as.integer(chr))

    out[[chr_str]] <- meta_chr_norm
  }

  dplyr::bind_rows(out)
}

# Functional status for the new model:
# F=1 for predicted functional, F=0 otherwise.

map_class <- function(x) {
  x2 <- toupper(trimws(as.character(x)))
  ifelse(x2 == "NS", 1L, ifelse(x2 == "S", 0L, NA_integer_))
}

## parse_ids
parse_ids <- function(s) {
  if (is.na(s) || !nzchar(s)) return(integer(0))
  as.integer(strsplit(s, ",\\s*")[[1]])
}

### encode configurations

encode_family_config <- function(v) {
  if (length(v) == 0L) "0" else paste(sort(v), collapse = "_")
}

### FIlter Invalid Variants after Mendelian computation

filter_invalid_by_mendel <- function(pieces_mendel_global, eps = 0) {
  vids <- names(pieces_mendel_global)

  bad <- vapply(vids, function(vid) {
    pm <- pieces_mendel_global[[vid]]
    any(!is.finite(c(pm$pG_half, pm$p0_half))) ||
      (pm$pG_half <= eps) ||
      (pm$p0_half < 0) ||
      (pm$p0_half >= 1)
  }, logical(1))

  list(
    keep_vids   = vids[!bad],
    drop_vids   = vids[ bad],
    kept_pieces = pieces_mendel_global[vids[!bad]]
  )
}

### Prune everything after we remove invalid variants

prune_to_vids <- function(out, keep_vids) {
  keep_vids <- intersect(keep_vids, out$var_ids_all)

  out$var_ids_all <- keep_vids

  out$meta_all <- out$meta_all[out$meta_all$variant_id %in% keep_vids, , drop = FALSE]
  out$class_vec_global <- out$class_vec_global[keep_vids]
  out$F_vec_global <- out$F_vec_global[keep_vids]
  out$H_raw_global <- out$H_raw_global[keep_vids]
  out$H_vec_global <- out$H_vec_global[keep_vids]

  out$global_var_info <- out$global_var_info[keep_vids]
  out$pieces_mendel_global <- out$pieces_mendel_global[keep_vids]

  out$pool$config_key_vec <- out$pool$config_key_vec[keep_vids]
  out$pool$unique_configs <- unique(out$pool$config_key_vec)
  out$pool$rep_vid_by_config <- tapply(
    names(out$pool$config_key_vec),
    out$pool$config_key_vec,
    `[`, 1
  )

  out
}

