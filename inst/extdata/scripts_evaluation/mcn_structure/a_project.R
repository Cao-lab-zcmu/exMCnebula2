# ==========================================================================
# external grob
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## external element
grob.lcms <- ex_grob("lcms")
grob.sample <- ex_grob("sample")
grob.convert <- ex_grob("convert")
grob.sirius <- ex_grob("sirius")
grob.collate <- ex_grob("collate_data")

# ==========================================================================
# sub.frame
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## melody: palette
type <- c("class", "slot", "sub.slot", "function", "custom")
pal <- MCnebula2:::.as_dic(palette_set(mcn_5features)[-3], type, na.rm = T)
pal[[ "report" ]] <- "black"

grobs.project <- lst_grecti("dataset", pal, , grecti2)

# ==========================================================================
# content
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## special arrow
spArr <- polygonGrob(
  c(0, 1, 1, .5, 0), c(1, 1, .4, 0, .4),
  gp = gpar(lwd = u(2, line), col = "grey30")
)

# alignment
shift <- rnorm(3, 2, 1)
all_range <- list(1:30, 31:60, 61:100)
set.seed(1)
lst <- mapply(shift, 1:length(shift), SIMPLIFY = F, FUN = function(shift, id){
  peak <- mapply(all_range, SIMPLIFY = F,
    FUN = function(range){
      peak <- dnorm(range, median(range) + shift, rnorm(1, 5, 1.2)) *
        rnorm(1, 0.7, 0.15)
    })
  feature <- mapply(1:length(all_range), lengths(all_range),
    FUN = function(seq, rep){
      rep(paste0("peak", seq), rep)
    })
  tibble::tibble(x = unlist(all_range), y = unlist(peak),
    sample = paste0("sample", id),
    peak = unlist(feature)
  )
})
data <- data.table::rbindlist(lst)
data <- dplyr::mutate(data,
  is = ifelse(y >= 0.003, T, F),
  peak = ifelse(is, peak, "non-feature"))
palette <- c(ggsci::pal_npg()(length(all_range)), "transparent")
names(palette) <- c(paste0("peak", 1:length(all_range)), "non-feature")
p <- ggplot(data, aes(x = x, y = y)) +
  geom_area(aes(fill = peak)) +
  geom_line() +
  geom_vline(xintercept = vapply(all_range, median, 1), linetype = "dashed") +
  geom_hline(yintercept = 0.003, linetype = "dashed") +
  facet_grid(sample ~ .) +
  scale_fill_manual(values = palette) +
  labs(x = "retention time", y = "intensity") +
  theme_void() +
  theme(
    text = element_text(family = "Times"), axis.text.y = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )
grob.align <- as_grob(p)
# detection
all_time <- vapply(all_range, median, 1)
anpi <- 0.05
len <- 5
set.seed(1)
ms1_set <- lapply(all_time,
  function(n){
    ms1 <- c(rnorm(3, 3, 2), rnorm(3, 8, 3), rnorm(3, 3, 2))
    ms1 <- ms1[ ms1 > 0 ]
    ms1 <- data.frame(x = 1:length(ms1), xend = 1:length(ms1),
      y = 0, yend = ms1)
    ms1 <- dplyr::mutate(ms1,
      x = x * len, xend = xend * len,
      y = sinpi(anpi) * x + y, yend = sinpi(anpi) * xend + yend,
      x = cospi(anpi) * x + n, xend = x)
    ms1 <- dplyr::bind_rows(ms1,
      c(x = n, xend = max(ms1$xend),
        y = 0, yend = max(ms1$y)))
    dplyr::mutate(ms1, time = n)
  })
ms1_set <- data.table::rbindlist(ms1_set)
p <- ggplot(dplyr::filter(data, sample == "sample1"), aes(x = x, y = y * 50)) +
  geom_segment(data = ms1_set,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "grey30", size = 0.8) +
  geom_area(aes(fill = peak)) +
  geom_line() +
  labs(x = "retention time", y = "intensity") +
  scale_fill_manual(values = palette) +
  theme_void() +
  theme(text = element_text(family = "Times"),
    axis.text.y = element_blank(),
    legend.position = "none"
  )
grob.detect <- as_grob(p)

## project_dataset
## into
grobs.project$dataset %<>% into(grob.collate)

## omit
omit <- circleGrob(seq(.2, .8, , 3), .5, .07, gp = gpar(fill = "grey20"),
                   vp = viewport(, , u(1.5, line), u(1.5, line)))

# ==========================================================================
# gather
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

grob.sampleCir <- ggather(circleGrob(), grob.sample)
spArrScale <- ggather(spArr, vp = viewport(, , .2, .7))
grob.rdetect <- ggather(
  into(grectn(bgp_args = list(lty = "solid")), grob.detect),
  vp = viewport(, , .7, .9)
)

frame.project <- frame_row(
  c(grob.sampleCir = .5, samples = .1, spArrScale = .1,
    grob.lcms = 1, `LC-MS/MS` = .1, spArrScale = .1,
    grob.convert = .6, convert_raw_data = .1, spArrScale = .1,
    grob.rdetect = .7, feature_detection = .1, spArrScale = .1,
    grob.sirius = .4, run_SIRIUS = .1, spArrScale = .1,
    MCnebula2 = .2, omit = .1, dataset = 1),
  c(namel(grob.sampleCir, grob.lcms, grob.convert, grob.rdetect,
      grob.sirius, spArrScale, omit),
    dataset = list(grobs.project[[1]]),
    sapply(c("samples", "LC-MS/MS", "convert_raw_data",
        "feature_detection", "run_SIRIUS", "MCnebula2"),
      simplify = F, function(name) { gtext(name) })
  )
)
.gene.vp <- viewport(, , width = unit(8 * 1.5, "line"),
                     height = unit(35 * 1.5, "line"), clip = "off")
.project <- gTree(children = gList(frame.project), vp = .gene.vp)

# ==========================================================================
# line and arrow
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

baf <- function(x, y, width = u(1, line), height = u(3, line)) {
  rect <- grectn(bgp_args = gpar(lty = "solid"))@grob
  clip <- clipGrob(, , .7)
  ggather(clip, rect, vp = viewport(x, y, width, height))
}
# a1. <- setnullvp("dataset", list(x = 1), .project)
# baf. <- baf(grobX(a1., 0), grobY(a1., 0))

