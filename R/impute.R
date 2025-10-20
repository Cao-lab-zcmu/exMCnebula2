impute_lcms <- function(df, method = c("knn", "median", "min")) {
  method <- match.arg(method)

  # 1️⃣ 缺失率检测
  na_rate_before <- colMeans(is.na(df))
  if (all(na_rate_before == 0)) {
    message("✅ No missing values detected. Skipping imputation.")
    return(df)
  }

  # 绘制插补前缺失率分布
  plot_before <- ggplot(data.frame(feature = names(na_rate_before), na_rate = na_rate_before)) +
    geom_histogram(aes(x = na_rate), bins = 30, fill = "#69b3a2", color = "white") +
    labs(title = "Missing rate before imputation", x = "Missing rate", y = "Count") +
    theme_minimal(base_size = 14)

  # 2️⃣ 根据方法进行插补
  df_imputed <- df
  if (method == "median") {
    df_imputed <- apply(df, 2, function(x) {
      ifelse(is.na(x), median(x, na.rm = TRUE), x)
    })
  } else if (method == "min") {
    df_imputed <- apply(df, 2, function(x) {
      ifelse(is.na(x), min(x, na.rm = TRUE), x)
    })
  } else if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE)) {
      install.packages("impute")
    }
    df_imputed <- impute::impute.knn(as.matrix(df))$data
  }
  df_imputed <- as.data.frame(df_imputed)

  # 3️⃣ 插补后缺失率
  na_rate_after <- colMeans(is.na(df_imputed))

  plot_after <- ggplot(data.frame(feature = names(na_rate_after), na_rate = na_rate_after)) +
    geom_histogram(aes(x = na_rate), bins = 30, fill = "#404080", color = "white") +
    labs(title = "Missing rate after imputation", x = "Missing rate", y = "Count") +
    theme_minimal(base_size = 14)

  # 4️⃣ 并排显示前后图
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(patchwork)
  print(plot_before + plot_after)

  # 返回插补后的数据框
  return(df_imputed)
}
