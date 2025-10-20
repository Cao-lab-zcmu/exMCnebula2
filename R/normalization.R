# ===========================================
# normalization.R
# ===========================================

# ---- 1. Sample normalization ----
normalize_sample <- function(df, method = "PQN", ref_group = NULL) {
  message("## [1] Sample normalization (method = ", method, ")")

  if (method == "sum") {
    norm_factors <- colSums(df, na.rm = TRUE)
    df <- sweep(df, 2, norm_factors / median(norm_factors), "/")
    message("\t- Normalized by total ion current (sum normalization)")

  } else if (method == "median") {
    norm_factors <- apply(df, 2, median, na.rm = TRUE)
    df <- sweep(df, 2, norm_factors / median(norm_factors), "/")
    message("\t- Normalized by sample median intensity")

  } else if (method == "PQN") {
    if (is.null(ref_group)) {
      ref_spectrum <- apply(df, 1, median, na.rm = TRUE)
      message("\t- Reference spectrum: median of all samples")
    } else {
      ref_idx <- grep(ref_group, colnames(df))
      ref_spectrum <- apply(df[, ref_idx, drop = FALSE], 1, median, na.rm = TRUE)
      message("\t- Reference spectrum: group = ", ref_group)
    }

    ratios <- apply(df, 2, function(x) x / ref_spectrum)
    pqn_factors <- apply(ratios, 2, median, na.rm = TRUE)
    df <- sweep(df, 2, pqn_factors, "/")
    message("\t- PQN normalization complete")
  } else {
    message("\t! Unknown method: ", method, ". Skipping normalization.")
  }

  return(df)
}

# ---- 2. Data transformation ----
transform_data <- function(df, method = "log2") {
  message("## [2] Data transformation (method = ", method, ")")

  if (method == "log10") {
    df <- log10(df + 1)
    message("\t- Applied log10(x + 1) transformation")

  } else if (method == "log2") {
    df <- log2(df + 1)
    message("\t- Applied natural log(x + 1) transformation")

  } else if (method == "cubeRoot") {
    df <- (df)^(1/3)
    message("\t- Applied cube root transformation")

  } else if (method == "none") {
    message("\t- No transformation applied")

  } else {
    message("\t! Unknown method: ", method, ". Skipping transformation.")
  }

  return(df)
}

# ---- 3. Scaling ----
scale_data <- function(df, method = "pareto") {
  message("## [3] Data scaling (method = ", method, ")")

  if (method == "pareto") {
    scaled_df <- scale(df, center = TRUE, scale = sqrt(apply(df, 2, sd, na.rm = TRUE)))
    message("\t- Applied Pareto scaling (mean-centered, divided by sqrt(SD))")

  } else if (method == "auto") {
    scaled_df <- scale(df, center = TRUE, scale = TRUE)
    message("\t- Applied auto scaling (mean-centered, divided by SD)")

  } else if (method == "range") {
    mins <- apply(df, 2, min, na.rm = TRUE)
    maxs <- apply(df, 2, max, na.rm = TRUE)
    scaled_df <- sweep(df, 2, mins, "-")
    scaled_df <- sweep(scaled_df, 2, maxs - mins, "/")
    message("\t- Applied range scaling (0-1 normalization)")

  } else if (method == "none") {
    scaled_df <- df
    message("\t- No scaling applied")

  } else {
    message("\t! Unknown method: ", method, ". Skipping scaling.")
    scaled_df <- df
  }

  return(scaled_df)
}

# ---- Main normalization pipeline ----
normalize_metabo_data <- function(df,
                                  sample_norm = "PQN",
                                  trans_method = "log2",
                                  scale_method = "pareto",
                                  ref_group = NULL) {
  message("===== Starting normalization pipeline =====")

  n_before <- nrow(df)

  df <- normalize_sample(df, method = sample_norm, ref_group = ref_group)
  df <- transform_data(df, method = trans_method)
  df <- scale_data(df, method = scale_method)

  n_after <- nrow(df)
  message("Normalization summary:")
  message("  Features processed: ", n_before)
  message("  Features after normalization: ", n_after)
  message("===== Normalization complete =====")

  return(df)
}


test <- normalize_metabo_data(test, ref_group = "^Sham")
