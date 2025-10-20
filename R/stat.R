#' export
library(pheatmap)
library(dendextend)
library(RColorBrewer)
library(MCnebula2)
library(exMCnebula2)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

mcn <- readRDS("data_mcnebula_mzmine/neg/mcn.rds")
files <- list.files("data_mzmine/neg", pattern = "MSMS\\.csv$", full.names = TRUE)
origin <- data.table::fread(files)
origin <- tibble::as_tibble(origin)
head(origin)

quant <- dplyr::select(
  origin, id = 1, dplyr::contains("Peak area")
)
colnames(quant) <- gsub("\\.mzML Peak area", "", colnames(quant))
quant <- dplyr::mutate(quant, .features_id = as.character(id))
head(quant)
gp <- c(Blank = "^Sham", Model = "^M")
metadata <- MCnebula2:::group_strings(colnames(quant), gp, "sample")

classes_heatmap <- c("Indoles and derivatives")
ids_heatmap <- select_features(
  mcn, classes_heatmap, q.value = 0.05, logfc = 0.3
)
quant <- dplyr::filter(quant, .features_id %in% ids_heatmap[[classes_heatmap]])
rownames(quant) <- quant$.features_id #quant行名无法修改就修改矩阵
sample_col <- data.frame(sample = metadata$group)
rownames(sample_col) <- metadata$sample
class_row <- data.frame(class = rep(classes_heatmap, length(quant[[1]])))#数字对应quant
rownames(class_row) <- rownames(quant)
quant <- dplyr::select(quant, -(.features_id | id))
quant_matrix <- t(apply(quant, 1, cal_z_score))
rownames(quant_matrix) <- rownames(class_row)

ann_colors <- list(sample = c(Blank = "#0073C299", Model = "#EFC00099"),
                   class = c("Indoles and derivatives" = "#D2B4DE"))

class_pheatmap <- pheatmap(
  quant_matrix,
  color = colorRampPalette(c("navy", "#FEF9E7", "firebrick3"))(500),
  annotation_col = sample_col, 
  annotation_row = class_row,
  cutree_cols = 2,
  cutree_rows = 2,
  border_color = "black",
  fontsize_row = 8,
  fontsize_col = 7,
  annotation_colors = ann_colors
)
save_pheatmap_png <- function(x, filename, width = 4000, height = 6000, res = 600) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(class_pheatmap, paste0(classes_heatmap, ".png"))

