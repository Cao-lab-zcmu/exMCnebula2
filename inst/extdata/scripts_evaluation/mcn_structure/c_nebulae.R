# ==========================================================================
# sub.frame
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

weight.nebulae <- c(parent_nebula = .7, child_nebulae = 1.5)
grobs.nebulae <- lst_grecti(names(weight.nebulae), pal, "slot", grecti2)

# ==========================================================================
# content
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## parent_nebula
keep <- c("grid_layout", "viewports", "ggset", "...")
name.par <- c("...", slotNames(parent_nebula(toBinary5)))
fun <- function(names) {
  names <- names[ names %in% keep ]
  weight <- rep(10, length(names))
  weight[which(names == "...")] <- 2
  names(weight) <- names
  weight
}
weight.par <- fun(name.par)
## draw
grobs.par <- lst_grecti(names(weight.par), pal, "sub.slot", grecti2)
omit <- circleGrob(seq(.2, .8, , 3), .5, .07, gp = gpar(fill = "grey20"),
                   vp = viewport(, , u(1.5, line), u(1.5, line)))
grobs.par$... <- omit
## ggset
ids <- n(n, 5)
set.seed(99)
graph <- random_graph(ids, layout = "kk")
theme <- theme_void() + theme(legend.position = "none")
layer0 <- ggraph::ggraph(graph)
scale_size <- scale_size(range = c(1, 3))
layer1 <- layer0 + ggraph::geom_node_point(aes(size = size), shape = 21) +
  scale_size + theme
layer1 <- as_grob(layer1)
edge_color <- "lightblue"
layer2 <- layer0 +
  ggraph::geom_edge_fan(aes(edge_width = width), color = edge_color) + theme
layer2 <- as_grob(layer2)
layer4 <- ggplot(data.frame(x = 1:2, y = 1:2, fill = n(g, 2))) + 
  geom_tile(aes(x = x, y = y, fill = fill)) + 
  ggsci::scale_fill_npg() + theme
layer4 <- as_grob(layer4)
layer5 <- layer0 + 
  ggraph::geom_edge_fan(aes(edge_width = width), color = edge_color) +
  ggraph::geom_node_point(aes(size = size, fill = name), shape = 21) +
  scale_size +
  ggsci::scale_fill_npg() +
  theme
layer5 <- layer5_ <- ggather(as_grob(layer5))
## ggset gather
layers <- lapply(namel(layer1, layer2, layer4, layer5),
                 function(lay){
                   into(grectn(), lay)
                 })
wei_ggset <- c(layer1 = 1, plus = .5, layer2 = 1, plus = .5, 
               layer4 = 1, equals = 1, layer5 = 1)
grob_ggset <- frame_col(wei_ggset, c(layers, list(plus = gtext("+"),
                                                  equals = gtext("=>"))))
grob_ggset <- ggather(grob_ggset, vp = viewport(, , .9, .5))
grobs.par$ggset %<>% into(grob_ggset)
## parent gather
.par <- frame_row(weight.par, grobs.par)
.par <- ggather(.par, vp = viewport(, , .size_subslot, .size_subslot))
## into
grobs.nebulae$parent_nebula %<>% into(.par)

## child_nebulae
name.chi <- c("...", slotNames(child_nebulae(toBinary5)), "...")
weight.chi <- fun(name.chi)
## draw
grobs.chi <- lst_grecti(names(weight.chi), pal, "sub.slot", grecti2)
grobs.chi$... <- omit
## signal grid layout
n <- 3
seq <- seq(.2, .8, length.out = n)
grid <- segmentsGrob(c(rep(0, n), seq),
                     c(seq, rep(1, n)),
                     c(rep(1, n), seq),
                     c(seq, rep(0, n)), gp = gpar(lty = "dashed"))
grid <- ggather(clipGrob(height = .8), rectGrob(), grid, vp = viewport(, , .8, .8))
grobs.chi$grid_layout %<>% into(grid)
## signal viewport
n <- 6
seq <- seq(135, 0, length.out = n)
seq2 <- seq(45, 0, length.out = n)
vps <- lapply(1:length(seq),
              function(n, fx = .5, fy = .5, x = .5) {
                ang <- seq[n]
                ang2 <- seq2[n]
                vp <- viewport(cospi(ang / 180) * fx + x, sinpi(ang / 180) * fy + x,
                               u(2, line), u(2, line), angle = ang2,
                               just = c("centre", "bottom"))
                ggather(rectGrob(gp = gpar(fill = "lightyellow", col = pal[["sub.slot"]])),
                        gtext(paste0("n", rev(1:length(seq))[n]),
                              gpar(fontface = "plain"), form = F), vp = vp)
              })
vps <- do.call(ggather, c(vps, list(vp = viewport(, .1, .5, .5))))
grobs.chi$viewports %<>% into(vps)
## signal ggset
ggsets <- frame_col(c(n = 1, x = 1, mn = 1),
                    list(n = gtext("n", list(cex = 2.5), form = F),
                         x = gtext("×", list(cex = 2, font = c(plain = 1))),
                         mn = into(glayer(3), layer5_)))
ggsets <- ggather(ggsets, vp = viewport(, , .7, .7))
grobs.chi$ggset %<>% into(ggsets)
## child... gather
.chi <- frame_row(weight.chi, grobs.chi)
.chi <- ggather(.chi, vp = viewport(, , .size_subslot, .size_subslot))
## ino
grobs.nebulae$child_nebulae %<>% into(.chi)

## visualization
load(paste0(.expath, "/toActiv30.rdata"))
test1 <- toActiv30
set.seed(10)
## filter
reference(mcn_dataset(test1))$nebula_index %<>%
  dplyr::filter(class.name %in% sample(unique(class.name), 9))
test1 <- create_parent_nebula(test1, 0.01, 5, T)
test1 <- create_child_nebulae(test1, 0.01, 5)
test1 <- create_parent_layout(test1)
test1 <- create_child_layouts(test1)
test1 <- activate_nebulae(test1)
## gather child grobs
chAsGrob <- function(ch, x) {
  ggset <- modify_default_child(ch)
  ggset@layers$ggtitle@command_args$label %<>%
    gsub("[a-zA-Z]", ".", .) %>% 
    gsub("^....", "Class", .)
  as_grob(call_command(ggset))
}
sets <- lapply(ggset(child_nebulae(test1)), chAsGrob, x = test1)
sets <- lapply(names(sets),
               function(name) {
                 ggather(sets[[name]],
                         vp = viewports(child_nebulae(test1))[[name]])
               })
sets_vp <- viewport(layout = grid_layout(child_nebulae(test1)))
sets <- do.call(ggather, c(sets, list(vp = sets_vp)))
legendH <- MCnebula2:::.legend_hierarchy(child_nebulae(test1), test1)
export_name(test1)[["tani.score"]] <- "TS"
export_name(test1)[["similarity"]] <- "SS"
legendG <- ggset(child_nebulae(test1))[[1]] %>%
  modify_default_child(x = test1) %>%
  call_command() %>%
  MCnebula2:::.get_legend()
## integrate
vis <- frame_row(list(sets = 5, legendH = .5), namel(sets, legendH))
vis <- frame_col(list(vis = 4, legendG = 1), namel(vis, legendG))
## final
.grob.vis <- vis

## report
rept <- grecti("Report: Workflow", tfill = pal[["report"]],
               tgp_args = list(col = "transparent"), borderF = 1.5)
head <- gtext("4.4 Create Nebulae", x = 0, hjust = 0)
head2 <- gtext("4.5. Visualize Nebulae", x = 0, hjust = 0)
fig <- into(grectn(), gtext("Fig.", list(cex = 3)))
des <- gtext("Description...", x = .1, hjust = 0)
code <- into(grectn("grey95", , list(col = "transparent")),
             gtext("## code",
                   list(font = c(plain = 1)), x = .1, hjust = 0))
mores <- gltext("More sections")
cont <- frame_row(c(mores = .3, head = .3, des = .2, code = .3,
                    head2 = .3, fig = 1, des = .2, code = .3, mores = .3),
                  namel(head2, head, fig, des, code, mores))
cont <- ggather(cont, vp = viewport(, , .8, .8))
rept <- into(rept, cont)
## final
.grob.report <- rept

# ==========================================================================
# gather
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## nebulae # .size_subslot <- u(1, npc) - u(.5, line)
weight.nebulae[[ "ne_null" ]] <- .1
weight.nebulae <- weight.nebulae[c(1, 3, 2)]
.nebulae <- frame_row(weight.nebulae, c(grobs.nebulae, list(ne_null = nullGrob())))
## more
.vis <- frame_col(c(slotNebula = 1, null = .3, report = 1),
                  list(slotNebula = .nebulae, report = .grob.report, null = nullGrob()))
.vis <- frame_row(c(vis1 = 1, vis_null = .05, vis2 = 1),
                  list(vis1 = .vis, vis2 = .grob.vis, vis_null = nullGrob()))
## panel vp
.gene.vp.2w <- .gene.vp
.gene.vp.2w$width <- .gene.vp$width * 2
.nebulae <- ggather(.vis, vp = .gene.vp.2w)

## show
# x11(, 7, 11)
# draw(.nebulae)

# ==========================================================================
# line and arrow
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## link
link <- c(
  "grid_layout", "vis2", "visualize_all",
  "viewports", "vis2", "visualize_all",
  "ggset", "vis2", "visualize_all"
)
link <- split(link, rep(1:(length(link) / 3), each = 3))
link <- data.frame(t(do.call(cbind, link)))
names(link) <- c("from", "to", "fun")
## attributes
fun_pal <- ggsci::pal_d3("category20")(20)
fun_pal <- fun_pal[!fun_pal %in% unique(link_b$color)]
link <- dplyr::mutate(link,
  group = agroup(to, 1:10, integer(1)),
  color = agroup(group, fun_pal),
  left = F, up = F,
  shift = 2,
  dup = duplicated(group)
)

## tools
tools_c <- maparrow(.nebulae, link, list(ggset = "child_nebulae.*::ggset"))
link <- tools_c$data

## bafs
bafs <- unlist(tools_c$bafs, F)
bafs <- lapply(bafs, function(fun) fun(u(1, line), u(1, line)))
bafs$vis2.r <- NULL
bafs$vis2 <- ({
  function(){
    g <- tools_c$nulls$vis2$t
    rectGrob(grobX(g, 0), grobY(g, 0), u(25, line), u(.05, line),
             gp = gpar(col = "transparent", fill = link$color[1]))
  }})()

## arrow
arr <- apply(link, 1, simplify = F,
             function(v) {
               pos <- ifelse(as.logical(v[[ "left" ]]), "l", "r")
               cur <- if (pos == "l") 1 else -1
               cur <- if (as.logical(v[[ "up" ]])) -cur else cur
               garrow(tools_c$nulls[[ v[["from"]] ]][[ pos ]],
                      tools_c$nulls[[ v[["to"]] ]][[ "t" ]],
                      list(curvature = cur, square = T, ncp = 1,
                           gp = gpar(fill = v[["color"]], col = v[["color"]],
                                     lwd = u(1, line)))
               )
             })

## sagnage
cvp <- garrow_city(tools_c$nulls$ggset$r, tools_c$nulls$ggset$t, F, F,
                   u(link$shift[1], line), gpar())@cvp
sag <- tools_c$sags[[1]]
sag$vp[[1]] <- cvp

.nebulae <- do.call(ggather, c(list(.nebulae), bafs, arr, list(sag)))
# draw(.nebulae)

link_c <- link


