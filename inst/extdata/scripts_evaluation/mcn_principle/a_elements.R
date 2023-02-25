# ==========================================================================
# grobs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## object
load(paste0(.expath, "/toActiv30.rdata"))
test1 <- toActiv30
ids <- features_annotation(test1)$.features_id
export_path(test1) <- paste0(tempdir(), "/mcnebula_results")
export_path(test1)
the_id <- "2095"
test1 <- draw_structures(test1, .features_id = the_id)

## get classes dataset
class <- latest(filter_ppcp(test1, pp.threshold = 0))
class <- dplyr::mutate(class, class.name = reorder(class.name, dplyr::desc(pp.value)))
class_2095 <- dplyr::filter(class, .features_id == !!the_id)
class_2095 <- dplyr::arrange(class_2095, class.name)
## modified the data
class <- dplyr::select(class, .features_id, class.name, pp.value)
class <- lapply(split(data.frame(class), ~.features_id),
  function(df){
    df <- dplyr::mutate(df, class.name = reorder(class.name, dplyr::desc(pp.value)))
    head(dplyr::arrange(df, class.name), n = 30)
  })
class <- data.table::rbindlist(class)

## for inner filter
class_in <- dplyr::filter(class, .features_id %in% as.character(c(2093, 2095, 2097)),
  pp.value <= 0.9)
class_in <- dplyr::bind_rows(class_in,
  data.frame(class.name = paste0("class", 1:3),
    pp.value = 0.5,
    .features_id = "Others"
    ))

## draw horizontal bar plot for inner filter
y.arrow <- seq(-0.4, 0.4, length.out = 3)
check_pal <- c(true = "#91D1C2FF", false = "#E64B35FF")
df <- dplyr::filter(class_in, .features_id != "Others")
df <- dplyr::mutate(df, .features_id = paste0("ID: ", .features_id))
p1 <- ggplot(df) +
  geom_col(aes(x = class.name, y = pp.value,
      fill = ifelse(pp.value >= 0.5, "true", "false"),
      )) +
  geom_segment(data = dplyr::filter(df, !grepl("[0-9]", class.name), pp.value >= 0.5),
    aes(x = class.name, xend = class.name, y = 1, yend = 1.15),
    color = check_pal[["true"]],
    arrow = arrow(8, unit(0.1, "inches"), "last", "closed")) +
  geom_point(data = dplyr::filter(df, grepl("[0-9]", class.name) | pp.value < 0.5),
    aes(x = class.name, y = 1),
    shape = 4, color = check_pal[["false"]]) +
  geom_text(data = dplyr::filter(df, grepl("[0-9]", class.name)),
    aes(x = class.name, y = -0.05, label = stringr::str_wrap(class.name, 30)),
    color = "red",
    family = "Times", hjust = 1, size = 1.2) +
  geom_text(data = dplyr::filter(df, !grepl("[0-9]", class.name)),
    aes(x = class.name, y = -0.05, label = stringr::str_wrap(class.name, 30)),
    color = "black",
    family = "Times", hjust = 1, size = 1.2) +
  geom_hline(aes(yintercept = 0.5),
    linetype = "dashed", color = "red") +
  coord_flip() +
  ylim(c(-1.3,1.3)) +
  labs(x = "Classification") +
  scale_fill_manual(values = check_pal) +
  theme_minimal() +
  facet_grid(.features_id ~ ., scales = "free") +
  theme(text = element_text(family = "Times", face = "bold"),
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    legend.position = "none") +
  geom_blank()
## arrow
p1.ex <- ggplot() +
  geom_segment(data = dplyr::filter(class_in, .features_id == "Others"),
    aes(x = 5, xend = 1, y = y.arrow, yend = y.arrow),
    linetype = "dashed",
    arrow = arrow(10, unit(0.15, "inches"), "last", "closed")) +
  coord_flip() +
  ylim(c(-1.3,1.3)) +
  labs(y = "Posterior probability") +
  scale_fill_manual(values = check_pal) +
  theme_minimal() +
  facet_grid(.features_id ~ ., scales = "free") +
  theme(text = element_text(family = "Times", face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    legend.position = "none") +
  geom_blank()
## as grob
p1 <- as_grob(p1)
p1.ex <- as_grob(p1.ex)
p1 <- frame_row(c(p1 = 3, p1.ex = .4), namel(p1, p1.ex))

## for stardust
stardust <- dplyr::filter(class, pp.value >= 0.5)
stardust <- dplyr::filter(class, !grepl("[0-9]", class.name))
stardust <- split(data.frame(stardust), ~ as.character(class.name))
hierarchy <- hierarchy(test1) %>% 
  dplyr::select(hierarchy, class.name)
index <- as.list(hierarchy$hierarchy)
names(index) <- hierarchy$class.name
label_color <- palette_label(test1)
## draw network plots
library(ggraph)
n <- 8
lst <- mapply(head(stardust, n = n), names(stardust)[1:n],
  SIMPLIFY = F,
  FUN = function(df, name){
    edges <- data.frame(.features_id1 = df[1, ]$.features_id,
      .features_id2 = df[1. ]$.features_id)
    graph <- igraph::graph_from_data_frame(edges, directed = T, vertices = df)
    graph <- tidygraph::as_tbl_graph(graph)
    layout_n <- create_layout(graph, layout = "kk")
    ggraph(layout_n) +
      geom_node_point(shape = 21, stroke = 0.2, fill = "#D9D9D9", size = 5) +
      theme_void() +
      ggtitle(stringr::str_wrap(name, width = 30)) +
      theme(text = element_text(family = "Times", face = "bold"),
        plot.title = ggtext::element_textbox(
          size = 10, color = "white",
          fill = label_color[ index[[name]] ],
          box.color = "white",
          halign = 0.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3))
        ) +
      geom_blank()
  })
## as grob
lst <- lapply(lst, as_grob)

## ring diagrame
text_gp <- gpar(color = "black", cex = 1, fontfamily = "Times", fontface = "bold",
  fontsize = 20)
p.ring <- ggplot(data.frame(x = 1, y = c(0.1, 0.9), group = c("true", "false"))) +
  geom_col(aes(x = x, y = y, fill = group)) +
  annotate("segment", x = 1.7, xend = 1.7, y = 0.01, yend = 0.08,
    arrow = arrow(angle = 10, type = "closed", length = unit(0.05, "npc"))) +
  annotate("segment", x = 1.6, xend = 1.8, y = 0, yend = 0) +
  annotate("segment", x = 1.6, xend = 1.8, y = 0.1, yend = 0.1) +
  scale_fill_manual(values = check_pal) +
  xlim(c(-1, 2)) +
  coord_polar(theta = "y") +
  ggtitle("Ratio") +
  theme_void() +
  theme(legend.position = "none",
    plot.title = element_text(family = .font, size = 6, hjust = .5)) +
  geom_blank()
## as grob
p.ring <- as_grob(p.ring)
## another ring...
p.ring2 <- ggplot(data.frame(x = 1, y = c(0.7, 0.3), group = c("true", "false"))) +
  geom_col(aes(x = x, y = y, fill = group)) +
  annotate("segment", x = 1.7, xend = 1.7, y = 0.01, yend = 0.68,
    arrow = arrow(angle = 10, type = "closed", length = unit(0.05, "npc"))) +
  annotate("segment", x = 1.6, xend = 1.8, y = 0, yend = 0) +
  annotate("segment", x = 1.6, xend = 1.8, y = 0.7, yend = 0.7) +
  scale_fill_manual(values = check_pal) +
  xlim(c(-1, 2)) +
  coord_polar(theta = "y") +
  ggtitle("Score distribution") +
  theme_void() +
  theme(legend.position = "none",
    plot.title = element_text(family = .font, size = 6, hjust = .5)) +
  geom_blank()
## as grob
p.ring2 <- as_grob(p.ring2)

## bar diagrame
p.bar <- ggplot(data.frame(x = 1, y = c(0.1, 0.9), group = c("false", "true"))) +
  geom_col(aes(x = x, y = y, fill = reorder(group, desc(y))), width = 0.5) +
  annotate("segment", x = 1.7, xend = 1.7, y = 0.12, yend = 0.98,
    arrow = arrow(angle = 10, type = "closed", length = unit(0.05, "npc"))) +
  annotate("segment", x = 1.6, xend = 1.8, y = 0.1, yend = 0.1) +
  scale_fill_manual(values = check_pal) +
  xlim(c(0, 2)) +
  ggtitle("Number") +
  theme_void() +
  theme(legend.position = "none",
    plot.title = element_text(family = .font, size = 6, hjust = .5)) +
  geom_blank()
## as grob
p.bar <- as_grob(p.bar)

## peak simbolized for feature 
p.peak <- ggplot(data.frame(x = 1:20, y = dnorm(1:20, 10, 3)), aes(x = x, y = y)) +
  geom_line() +
  theme_minimal() +
  labs(x = "RT", y = "Intensity") +
  theme(text = element_blank())
p.peak.yellow <- p.peak + 
  geom_area(fill = "lightyellow") 
p.peak.grey <- p.peak +
  geom_area(fill = "#D9D9D9") 
## as grob
p.peak <- as_grob(p.peak)
p.peak.yellow <- as_grob(p.peak.yellow)
p.peak.grey <- as_grob(p.peak.grey)

## boxplot
struc_2095 <- structures_grob(child_nebulae(test1))[[ the_id ]]
set.seed(100)
simu_score <- data.frame(.features_id = 1:50, score = rnorm(50, 0.4, 0.15)) %>% 
  dplyr::mutate(group = ifelse(score >= 0.3, "true", "false"))
p.box <- ggplot(simu_score, aes(x = "", y = score)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
  geom_jitter(width = 0.1, aes(color = group)) +
  labs(x = "", y = "") +
  # ggtitle("Identified score") +
  scale_fill_manual(values = check_pal) +
  theme_void() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = .5),
    text = element_text(family = .font, size = 5),
    # panel.grid = element_blank(),
    # plot.background = element_blank(),
    # panel.background = element_blank(),
    axis.ticks = element_blank()) +
  geom_blank()
## as grob
p.box <- as_grob(p.box)

## ratio bar plot
ratio <- table(simu_score$group) %>% 
  prop.table() %>%
  as.data.frame()
p.ratio <- ggplot(ratio, aes(x = 1, y = Freq, fill = reorder(Var1, desc(Freq)))) +
  ggchicklet::geom_chicklet(radius = unit(0.1, "npc")) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "red") +
  scale_fill_manual(values = check_pal) +
  xlim(c(0, 2)) +
  theme_minimal() +
  theme(text = element_blank(),
    legend.position = "none") +
  geom_blank()
  ## as_grob
p.ratio <- as_grob(p.ratio)

## assessment
df <- dplyr::filter(stardust_classes(test1), hierarchy >= 3)
parents <- get_parent_classes(unique(df$class.name), test1)
set <- parents[!vapply(parents, is.null, logical(1))]
set <- set[vapply(set, function(cl) if (any(cl %in% unique(df$class.name) )) T else F, T)]
## draw
ec <- c("Carboxylic acids")
ec. <- "Child..."
ep <- c("Carboxylic acids and derivatives")
ep. <- "Parent..."
labels <- lapply(c(ec., ep., "Class..."),
  function(lab) {
    zo(into(grectn(bgp_args = list(lty = "solid")),
        gtextp(stringr::str_wrap(lab, 15), list(col = "white", cex = .6))))
  })
names(labels) <- c("ec", "ep", "omit")
labels$null <- nullGrob()
## frame
scolor <- function(grobs, hier, color = palette_label(test1)[hier]) {
  lapply(grobs,
    function(grob) {
      if (is(grob, "null")) return(grob)
      else {
        grob$children[[1]]$children[[1]]$gp$fill <- color
        grob$children[[1]]$children[[1]]$gp$col <- "transparent"
        return(grob)
      }
    })
}
frame1 <- frame_row(c(null = .3, omit = .4, null = .3, omit = .4, null = .3),
  scolor(labels, 2))
frame2 <- frame_row(c(omit = .5, null = .3, ep = 1, null = .3, omit = .5),
  scolor(labels, 3))
frame3 <- frame_row(c(null = .3, ec = 1, null = .3, omit = .5, null = .3),
  scolor(labels, 4))
frame <- frame_col(c(frame1 = 1, null = .6, frame2 = 1, null = .6, frame3 = 1),
  namel(frame1, frame2, frame3, null = nullGrob()))
## arrows
path <- c("frame1::.*::omit.rep.2", "frame2::.*::omit",
  "frame2::.*::ep", "frame2::.*::omit.rep.2", "frame3::.*::ec", "frame3::.*omit")
nulls <- lapply(path, function(p) {
  list(l = setnullvp(p, list(x = 0, y = .5), frame),
    r = setnullvp(p, list(x = 1, y = .5), frame))
  })
names(nulls) <- c("n1.2", "n2.1", "ep", "n2.3", "ec", "n3.2")
arr1_2 <- c("n2.1", "ep", "n2.3")
cur <- c(-.1, .2, .3)
arr_gpar <- gpar(fill = palette_label(test1)[2], col = palette_label(test1)[2])
arr1_2 <- lapply(1:length(arr1_2),
  function(n) {
    no <- arr1_2[n]
    garrow(nulls$n1.2$r, nulls[[no]]$l,
      list(curvature = cur[n], gp = arr_gpar))
  })
arr2_3 <- c("ec", "n3.2")
cur <- c(-.1, .2)
arr2_3 <- lapply(1:length(arr2_3),
  function(n) {
    no <- arr2_3[n]
    garrow(nulls$ep$r, nulls[[no]]$l, 
      list(curvature = cur[n], gp = arr_gpar))
  })
## gather
bg <- rectGrob(.3, , .7, 1.1, just = c("left", "centre"),
  gp = gpar(col = "transparent", fill = "lightyellow"))
hier_tree <- do.call(ggather, c(list(bg), list(frame), arr1_2, arr2_3))

## venn plot
conts <- list(dplyr::filter(df, class.name == ep)$.features_id,
  dplyr::filter(df, class.name == ec)$.features_id)
feas <- unique(c(conts[[1]], conts[[2]]))
data <- data.frame(a = feas %in% conts[[1]], b = feas %in% conts[[2]])
names <- stringr::str_trunc(c(ep, ec), 15, "cen")
cir <- ggather(circleGrob(.3, gp = gpar(alpha = .5, fill = palette_label(test1)[3])),
  circleGrob(.7, gp = gpar(alpha = .5, fill = palette_label(test1)[4])))
venn.text <- frame_col(c(x = 1, null = 1),
  list(x = gtextp("Overlap Ratio", list(cex = .5), x = .8, y = -1),
    null = nullGrob()))
venn <- frame_row(c(null = .05, tit = .1, cont = 2, null = .05),
  list(tit = venn.text, cont = cir, null = nullGrob()))

## upset plot
set.seed(101)
rand <- sample(unique(df$class.name), 3)
comc <- c(ep, ec, rand)
data <- dplyr::filter(df, class.name %in% comc)
mutate <- .as_dic(c(ep., ec., paste0(n(Class, 2), "..."), "..."), c(ep, ec, rand), , F)
data <- dplyr::mutate(data, class.name =
  vapply(class.name, FUN.VALUE = character(1), USE.NAMES = F,
    function(cl) {
      do.call(switch, c(list(cl), mutate))
    }))
data <- lapply(split(data, ~.features_id), function(df) df$class.name)
data <- lapply(data,
  function(cl) {
    unlist(lapply(cl,
        function(cl) {
          mutate <- mutate[mutate != cl]
          paste0(mutate, "_", cl)
        }))
  })
data <- tibble::tibble(Class = unlist(data))
p.upset <- ggplot(data, aes(x = Class)) +
  geom_bar() +
  ggupset::axis_combmatrix(sep = "_", levels = unlist(mutate, use.names = F)) +
  theme(text = element_text(family = .font),
    axis.text.y = element_text(size = 5),
    plot.title = element_text(hjust = .5, size = 6),
    axis.title = element_text(size = 6)) +
  labs(y = "Number") +
  ggtitle("Overlap Number") +
  ggupset::theme_combmatrix(
    combmatrix.label.text = element_text(family = .font, size = 5)) +
  geom_blank()
p.upset <- as_grob(p.upset)


