# ==========================================================================
# sub.frame
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
load(paste0(.expath, "/toBinary5.rdata"))
weight.mcn <- weight(mcn_dataset(toBinary5), slotNames(mcn_dataset(toBinary5)))
weight.mcn <- 
  sapply(c("dataset", "reference", "backtrack"), simplify = F,
         function(x) weight.mcn[[x]])

grobs.mcn <- lst_grecti(names(weight.mcn), pal, "sub.slot", grecti2)

# ==========================================================================
# content
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## dataset
## features
mshn <- 5
mshvp <- viewport(, , .7, .7)
fea <- circleGrob(gp = gpar(fill = "grey85"))
fea <- ggather(fea, vp = mshvp)
grob_feas <- frame_col(fill_list(n(fea, mshn), 1), fill_list(n(fea, mshn), list(fea)))
## function to draw
mglayer <- function(from = 1:5, n = 5, seed = 100) {
  set.seed(seed)
  ns <- sample(from, n, T)
  grobs <- lapply(ns, function(n) {
                    ggather(glayer(n), vp = mshvp)
         })
  sig <- paste0(paste0("glayer", seed), 1:n)
  names(grobs) <- sig
  frame_col(fill_list(sig, 1), grobs)
}
## set
grob_fset <- mglayer(, , 100)
grob_sset <- mglayer(2:7, , 110)
grob_cset <- mglayer(5:10, , 120)
## titles
tit_fea <- gltext("features")
tit_cand <- gltext("candidates", flip = T)
tit_formu <- gtext("molecular formula", gpar(fontface = "plain"))
tit_struc <- gtext("chemical structure", gpar(fontface = "plain"))
tit_class <- gtext("chemical classes", gpar(fontface = "plain"))
## combine
names <- c("grob_feas", "tit_formu", "grob_fset",
           "tit_struc", "grob_sset", "tit_class", "grob_cset")
env <- environment()
mdataset <- frame_row(fill_list(names, 1), sapply(names, get, envir = env))
mdataset <- frame_col(list(tit_cand = .1, mdataset = 1),
                      list(tit_cand = tit_cand, mdataset = mdataset))
mdataset <- frame_row(list(tit_fea = .1, mdataset = 1),
                      list(tit_fea = tit_fea, mdataset = mdataset))
.grob.dataset <- ggather(mdataset, vp = viewport(, , .9, .9))
## into
grobs.mcn$dataset %<>% into(.grob.dataset)

## reference
ref <- names(reference(mcn_dataset(toBinary5)))
## spe.table
df <- data.frame(.fea_id = n(id., 3), .cand_id = n(cand., 3))
grob_table <- gridExtra::tableGrob(df, theme = gridExtra::ttheme_default(8))
grobs_ref <- sapply(ref, simplify = F,
                    function(name) {
                      if (name == "nebula_index")
                        res <- into(gshiny(xn = 6), gtext(name))
                      else if (name == "specific_candidate") {
                        res <- frame_row(list(tit = .2, df = 1),
                                         list(tit = gtext(name, y = .2), df = grob_table))
                        res <- into(grectn("transparent"), res)
                      } else
                        res <- into(grectn("transparent"), gtext(name))
                      ggather(res, vp = viewport(, , .9))
                    })
weight.ref <- vapply(ref, FUN.VALUE = numeric(1),
                     function(name) {
                       if (name == "nebula_index") .5
                       else if (name == "specific_candidate") 1.3
                       else .5
                     })
grobs_ref <- frame_row(weight.ref, grobs_ref)
grobs.ref <- ggather(grobs_ref, vp = viewport(, , .9, .9))
## into
grobs.mcn$reference %<>% into(grobs.ref)

## backtrack
grob.trash <- ex_grob("trash")
grob_trash <- ggather(grob.trash, vp = viewport(, , .8, .8))
grobs.mcn$backtrack %<>% into(grob_trash)

# ==========================================================================
# gather
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.size_subslot <- u(1, npc) - u(.5, line)
frame.mcn <- frame_row(weight.mcn, grobs.mcn)
frame.mcn <- ggather(frame.mcn, vp = viewport(, , .size_subslot, .size_subslot))
.mcn <- grecti2(form("mcn_dataset"), tfill = pal[[ "slot" ]])
.mcn <- into(.mcn, frame.mcn)
.mcn <- ggather(.mcn, vp = .gene.vp)

# draw(.mcn)

# ==========================================================================
# line and arrow
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## link
link <- c(
  # "dataset", "stardust_classes", "create_stardust_classes",
  # "dataset", "specific_candidate", "create_reference",
  # "dataset", "features_annotation", "create_features_annotation",
  # "specific_candidate", "features_annotation", "create_features_annotation",
  # "hierarchy", "stardust_classes", "create_stardust_classes",
  # "features_annotation", "nebula_index", "create_nebula_index",
  # "stardust_classes", "nebula_index", "create_nebula_index",
  "stardust_classes", "backtrack", "cross_filter_stardust",
  "backtrack", "stardust_classes", "backtrack_stardust"
)
link <- split(link, rep(1:(length(link) / 3), each = 3))
link <- data.frame(t(do.call(cbind, link)))
names(link) <- c("from", "to", "fun")
## the attributes for arrow
fun_pal <- ggsci::pal_d3("category20")(20)
link <- dplyr::mutate(link,
  group = agroup(to, 1:10, integer(1)),
  group = ifelse(from == "backtrack", max(group) + 1, group),
  color = agroup(group, fun_pal),
  left = ifelse(from %in% c("dataset", "specific_candidate"), F, T),
  left = ifelse(to %in% c("specific_candidate", "stardust_classes"),
    !left, left),
  up = ifelse(from == "backtrack", T, F),
  shift = ifelse(from %in% c("dataset", "backtrack"), 2, 2.5))

## tools
tools <- maparrow(.mcn, link)
link <- tools$data

## complex arr
com_arr <- lapply(1:length(tools$arr_city),
                  function(n) {
                    rect <- c("create_stardust_classes", "cross_filter_stardust",
                              "create_features_annotation", "backtrack_stardust")
                    if (link$fun[[ n ]] %in% rect & !link$dup[[ n ]])
                      tools$arr_city[[n]]
                    else
                      tools$arr_round[[n]]
                  })

## modify bafs
bafs <- sapply(names(tools$bafs), simplify = F,
               function(name){
                 sub_name <- gsub("\\.[a-z]$", "", name)
                 fun <- tools$bafs[[name]]
                 if (sub_name == "backtrack")
                   fun(h = u(2, line))
                 else if (sub_name == "specific_candidate")
                   fun(w = u(3, line))
                 else if (sub_name == "dataset")
                   fun()
                 else
                   fun(u(3, line), u(1, line))
               })

## show
.mcn <- do.call(ggather, c(list(.mcn), bafs, com_arr, tools$sags))
# draw(.mcn)

link_b <- link

