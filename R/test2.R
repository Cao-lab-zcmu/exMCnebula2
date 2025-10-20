devtools::load_all("~/git/MCnebula2/")
files <- list.files("data_mzmine", pattern = "msms_neg\\.csv$", full.names = TRUE)
origin <- data.table::fread(files)
origin <- tibble::as_tibble(origin)

quant <- dplyr::select(
  origin, id = 1, dplyr::contains("Peak area")
)
colnames(quant) <- gsub("\\.mzML Peak area", "", colnames(quant))
gp <- c(Sham = "^Sham", Model = "^M", QC = "^QC")
metadata <- MCnebula2:::group_strings(colnames(quant), gp, "sample")

df <- impute_mv(quant, "knn")

library(ggplot2)

df <- quant
  qc_cols <- grep(paste0("^", "QC"), colnames(df))
# 假设你的数据 df 和 qc_cols 已经定义
# 计算 RSD
rsd <- apply(df[, qc_cols, drop = FALSE], 1, function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(Inf)
  sd(x, na.rm = TRUE) / m * 100
})
rsd[is.na(rsd)] <- Inf

# 1. 散点图
rsd_df <- data.frame(
  feature = rownames(df),
  RSD = rsd
)

ggplot(rsd_df, aes(x = feature, y = RSD)) +
  geom_point(color = "steelblue") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "RSD Scatter Plot", y = "RSD (%)", x = "Features")

# 2. 累积分布图（cumulative plot）
rsd_sorted <- sort(rsd)
cum_df <- data.frame(
  RSD = rsd_sorted,
  CumFraction = seq_along(rsd_sorted) / length(rsd_sorted)
)

ggplot(cum_df, aes(x = RSD, y = CumFraction)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Cumulative RSD Distribution", x = "RSD (%)", y = "Cumulative Fraction")

scale_quant <- scale(log2(dplyr::select(quant, -id, -M5_NEG, -dplyr::starts_with("QC")) + 1), center = TRUE, scale = TRUE)
#scale_quant <- scale(log2(test) + 1, center = TRUE, scale = TRUE)
colnames(scale_quant)
scale_quant_t <- t(scale_quant)
#scale_quant_t <- t(test2)

cG_fecPosAll.pca <- FactoMineR::PCA(scale_quant_t, graph = FALSE)
metadata <- metadata[match(rownames(scale_quant_t), metadata$sample), ]

tiff('cG_fec_PCA_score_nor_scale.tiff', units = "in", width = 4.7, height = 4, res = 300)
factoextra::fviz_pca_ind(cG_fecPosAll.pca,
             geom.ind = c("point"),
             #geom.ind = c("point", "text"),
             #repel = TRUE,
             fill.ind = metadata$group,
             col.ind = "black",
             pointshape = 21, pointsize = 1.5,
             palette = "jco",
             addEllipses = TRUE,
             #ellipse.level = 0.5,
             #ellipse.type = "euclid",
             #ellipse.alpha = 0.25,
             #ellipse.linetype = 2,   # 虚线
             #ellipse.color = "grey40", # 灰色边框
             ggtheme = theme_gray(),
             legend.title = "Groups",
             title = "Serum metabolome ESI-",
             legend="right",
             font.legend = 9)
dev.off()

scores <- as.data.frame(cG_fecPosAll.pca$ind$coord[, 1:3])
colnames(scores) <- c("PC1", "PC2", "PC3")
scores$Group <- metadata$group

eig <- cG_fecPosAll.pca$eig[, 2]
p1 <- round(eig[1], 1)
p2 <- round(eig[2], 1)
p3 <- round(eig[3], 1)

plotly::plot_ly(
  data = scores,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Group,
  colors = RColorBrewer::brewer.pal(length(unique(scores$Group)), "Set2"),
  marker = list(size = 6, line = list(color = "black", width = 0.5))
) |>
  plotly::layout(
    scene = list(
      xaxis = list(title = paste0("PC1 (", p1, "%)")),
      yaxis = list(title = paste0("PC2 (", p2, "%)")),
      zaxis = list(title = paste0("PC3 (", p3, "%)"))
    ),
    title = "3D PCA Score Plot"
  )


# PCA 得分
scores <- as.data.frame(FactoMineR::PCA(scale_quant_t, graph = FALSE)$ind$coord[, 1:3])
colnames(scores) <- c("PC1", "PC2", "PC3")
scores$Group <- metadata$group
scores$Sample <- rownames(scores)

# 解释方差
eig <- FactoMineR::PCA(scale_quant_t, graph = FALSE)$eig[, 2]
p1 <- round(eig[1], 1)
p2 <- round(eig[2], 1)
p3 <- round(eig[3], 1)

ellipses <- data.frame()

for(grp in unique(scores$Group)){
  grp_data <- scores[scores$Group == grp, c("PC1", "PC2")]
  el <- ellipse::ellipse(cov(grp_data), centre=colMeans(grp_data), level=0.95, npoints=50)
  el_df <- data.frame(el)
  el_df$Group <- grp
  ellipses <- rbind(ellipses, el_df)
}

ellipses$PC3 <- 0  # 固定在 z=0 平面

colors <- RColorBrewer::brewer.pal(length(unique(scores$Group)), "Set2")

fig <- plotly::plot_ly()

# 样本点
fig <- fig %>% plotly::add_trace(
  data = scores,
  x = ~PC1, y = ~PC2, z = ~PC3,
  type = "scatter3d",
  mode = "markers",
  color = ~Group,
  colors = colors,
  marker = list(size = 6, line = list(color = "black", width = 0.5)),
  text = ~Sample,
  hoverinfo = "text+x+y+z"
)

# 椭圆
for(grp in unique(scores$Group)){
  el_df <- ellipses[ellipses$Group == grp, ]
  fig <- fig %>% plotly::add_trace(
    x = el_df$PC1,
    y = el_df$PC2,
    z = el_df$PC3,
    type = "scatter3d",
    mode = "lines",
    line = list(color = colors[which(unique(scores$Group)==grp)], width = 2),
    showlegend = FALSE
  )
}

# 坐标轴和标题
fig <- fig %>% plotly::layout(
  scene = list(
    xaxis = list(title = paste0("PC1 (", p1, "%)")),
    yaxis = list(title = paste0("PC2 (", p2, "%)")),
    zaxis = list(title = paste0("PC3 (", p3, "%)"))
  ),
  title = "3D PCA Score Plot with Ellipses"
)

fig

htmlwidgets::saveWidget(fig, "3D_PCA_Ellipse.html")

quant_long <- dplyr::select(log10(quant), -id) |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Sample", values_to = "Intensity")

ggplot(quant_long, aes(x = Sample, y = Intensity)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Boxplot per Sample", x = "Sample", y = "Intensity")

test <- processing_metabo_data(quant)
test2 <- normalize_metabo_data(dplyr::select(test, -id), ref_group = "^QC")

X <- scale_quant_t
Y <- as.factor(metadata$group)
plsda_res <- mixOmics::plsda(X, Y, ncomp = 3)

tiff("plsda_res.tiff", units = "in", width = 4.7, height = 4, res = 300)
mixOmics::plotIndiv(plsda_res,
          comp = c(1, 2),
          group = Y,
          ind.names = TRUE,
          ellipse = TRUE,
          legend = TRUE,
          title = "PLS-DA Score Plot")
dev.off()

df_pls <- data.frame(plsda_res$variates$X[, 1:2], Group = Y)
colnames(df_pls)[1:2] <- c("PLS1", "PLS2")
tiff("plsda_res_ggplot.tiff", units = "in", width = 4.7, height = 4, res = 300)
ggplot(df_pls, aes(x = PLS1, y = PLS2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  labs(title = "PLS-DA Score Plot", x = "PLS-DA Component 1", y = "PLS-DA Component 2") +
  scale_color_brewer(palette = "Set1")
dev.off()
