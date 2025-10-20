library(MCnebula2)
library(exMCnebula2)

files <- list.files("data_mzmine/neg", pattern = "MSMS\\.csv$", full.names = TRUE)
origin <- data.table::fread(files)
origin <- tibble::as_tibble(origin)

quant <- dplyr::select(
  origin, id = 1, dplyr::contains("Peak area")
)
colnames(quant) <- gsub("\\.mzML Peak area", "", colnames(quant))
quant <- dplyr::mutate(quant, .features_id = as.character(id))

quant <- dplyr::select(quant, -id, -".features_id")
gp <- c(Blank = "^Sham", Model = "^M")
metadata <- MCnebula2:::group_strings(colnames(quant), gp, "sample")

# ============================================
# Modular Filtering Function — Print the number of features filtered at each step
# Assumption: rows = features, columns = samples (including QC / BLANK columns)
# ============================================

# ---- Low-Quality Filtering (blank contamination + missing rate) ----
filter_low_quality <- function(df,
                               blank_cols,
                               sample_cols,
                               contam_cutoff = 0.5,
                               na_cutoff = 0.5) {
  n_init <- nrow(df)

  # blank 污染过滤（若提供 blank 列）
  if (length(blank_cols) > 0) {
    message("## for blank subtraction: remove if blank > ", contam_cutoff * 100, "% of sample signal")
    blank_mean <- rowMeans(df[, blank_cols, drop = FALSE], na.rm = TRUE)
    sample_mean <- rowMeans(df[, sample_cols, drop = FALSE], na.rm = TRUE)
    contam_ratio <- blank_mean / (sample_mean + 1e-9)

    n_before_blank <- nrow(df)
    keep_blank <- contam_ratio < contam_cutoff
    df <- df[keep_blank, , drop = FALSE]
    n_after_blank <- nrow(df)
    removed_blank <- n_before_blank - n_after_blank
    message("  -> blank filter removed ", removed_blank, " features (", 
            round(removed_blank / n_before_blank * 100, 1), "% of step) ",
            " | remaining: ", n_after_blank)
  } else {
    message("## no BLANK columns provided — skip blank contamination filtering")
  }

  # 缺失率过滤
  message("## remove features with missing rate >= ", na_cutoff)
  n_before_na <- nrow(df)
  na_rate <- rowMeans(is.na(df))
  keep_na <- na_rate < na_cutoff
  df <- df[keep_na, , drop = FALSE]
  n_after_na <- nrow(df)
  removed_na <- n_before_na - n_after_na
  message("  -> missing-rate filter removed ", removed_na, " features (",
          round(removed_na / n_before_na * 100, 1), "% of step) ",
          " | remaining: ", n_after_na)

  message("## low-quality filtering done: ", nrow(df), " features remain (from ", n_init, ")")
  return(df)
}


# ---- QC重复性过滤（RSD-based） ----
filter_rsd <- function(df, qc_cols, rsd_cutoff = 30) {
  n_before <- nrow(df)
  if (length(qc_cols) > 1) {
    message("## removing high-RSD features (> ", rsd_cutoff, "%) based on QC columns")
    # 计算 RSD（%），对可能的 NA 做兜底处理（NA -> Inf，表示去除）
    rsd <- apply(df[, qc_cols, drop = FALSE], 1, function(x) {
      m <- mean(x, na.rm = TRUE)
      if (is.na(m) || m == 0) return(Inf)
      sd(x, na.rm = TRUE) / m * 100
    })
    rsd[is.na(rsd)] <- Inf
    keep <- rsd < rsd_cutoff
    df <- df[keep, , drop = FALSE]
    n_after <- nrow(df)
    removed <- n_before - n_after
    message("  -> RSD filter removed ", removed, " features (", round(removed / n_before * 100, 1), "%) | remaining: ", n_after)
  } else {
    message("## no (or only one) QC column provided — skip RSD filtering")
  }
  return(df)
}


# ---- 方差过滤（按 sample_cols 计算 feature 行方差，去掉底部 var_cutoff 百分位） ----
filter_low_variance <- function(df, sample_cols, var_cutoff = 0.1) {
  n_before <- nrow(df)
  message("## removing low-variance features (bottom ", var_cutoff * 100, "% of variance)")
  var_value <- apply(df[, sample_cols, drop = FALSE], 1, var, na.rm = TRUE)
  var_threshold <- quantile(var_value, probs = var_cutoff, na.rm = TRUE)
  keep <- var_value > var_threshold
  df <- df[keep, , drop = FALSE]
  n_after <- nrow(df)
  removed <- n_before - n_after
  message("  -> variance filter removed ", removed, " features (threshold = ", signif(var_threshold, 4), ") | remaining: ", n_after)
  return(df)
}


# ---- 低丰度过滤（按行均值，去掉底部 abun_cutoff 百分位） ----
filter_low_abundance <- function(df, sample_cols, abun_cutoff = 0.1) {
  n_before <- nrow(df)
  message("## removing low-abundance features (bottom ", abun_cutoff * 100, "% by mean intensity)")
  mean_intensity <- rowMeans(df[, sample_cols, drop = FALSE], na.rm = TRUE)
  abun_threshold <- quantile(mean_intensity, probs = abun_cutoff, na.rm = TRUE)
  keep <- mean_intensity > abun_threshold
  df <- df[keep, , drop = FALSE]
  n_after <- nrow(df)
  removed <- n_before - n_after
  message("  -> abundance filter removed ", removed, " features (threshold = ", signif(abun_threshold, 4), ") | remaining: ", n_after)
  return(df)
}

clean_metabo_features <- function(df,
                                  qc_prefix = "QC",
                                  blank_prefix = "BLANK",
                                  contam_cutoff = 0.5,
                                  rsd_cutoff = 30,
                                  na_cutoff = 0.5,
                                  var_cutoff = 0.1,
                                  abun_cutoff = 0.1) {
  message("===== Starting feature filtering =====")

  # ---- 0. Preparation ----
  qc_cols <- grep(paste0("^", qc_prefix), colnames(df))
  blank_cols <- grep(paste0("^", blank_prefix), colnames(df))
  sample_cols <- setdiff(seq_len(ncol(df)), c(qc_cols, blank_cols))
  n_before <- nrow(df)

  # ---- 1. 低质量过滤 ----
  df <- filter_low_quality(df, blank_cols, sample_cols,
                           contam_cutoff = contam_cutoff,
                           na_cutoff = na_cutoff)

  # ---- 2. QC重复性过滤 ----
  df <- filter_rsd(df, qc_cols, rsd_cutoff = rsd_cutoff)

  # ---- 3. 方差过滤 ----
  df <- filter_low_variance(df, sample_cols, var_cutoff = var_cutoff)

  # ---- 4. 低丰度过滤 ----
  df <- filter_low_abundance(df, sample_cols, abun_cutoff = abun_cutoff)

  # ---- 总结 ----
  n_after <- nrow(df)
  message("Filtering summary:")
  message("  Features before filtering: ", n_before)
  message("  Features after filtering : ", n_after)
  message("  Retained proportion       : ", round(n_after / n_before * 100, 1), "%")
  message("===== Filtering complete =====")

  return(df)
}

test <- clean_metabo_features(quant, blank_prefix = "^Sham", contam_cutoff = 0.9)

# 假设 df 是你传入 filter 的原始数据（features x samples）
# blank_cols 与 sample_cols 与你之前一致
df=quant
blank_prefix="^Sham"
blank_mean <- rowMeans(quant[, blank_cols, drop = FALSE], na.rm = TRUE)
sample_mean <- rowMeans(quant[, sample_cols, drop = FALSE], na.rm = TRUE)
contam_ratio <- blank_mean / (sample_mean + 1e-9)

summary(contam_ratio)
summary(blank_mean)
summary(sample_mean)

# 比例直方图
hist(contam_ratio, breaks = 50, main = "blank/sample ratio", xlab = "blank_mean / sample_mean")

# blank 与 sample 的散点（若需要看绝对强度关系）
p <- plot(sample_mean, blank_mean, pch = 20,
     xlab = "sample_mean", ylab = "blank_mean",
     main = "sample_mean vs blank_mean")
abline(0,1, col="red", lty=2)
