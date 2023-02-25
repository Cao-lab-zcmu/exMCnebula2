# ==========================================================================
# all grobs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.over <- frame_col(list(a = .5, b = .5, c = 1),
                   list(a = .project, b = .mcn, c = .nebulae))

# ==========================================================================
# line and arrow
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## pal
fun_pal <- ggsci::pal_d3("category20")(20)
fun_pal <- fun_pal[!fun_pal %in% unique(c(link_b$color, link_c$color))]
fun_pal <- fun_pal[-1]

## a_dataset to b_dataset
# col <- fun_pal[1]
# a_dataset <- setnullvp("a::.*::dataset", list(x = 1, y = .5), .over)
# b_dataset <- setnullvp("b::.*::dataset", list(x = 0, y = .55), .over)
# arr <- garrow_snake(a_dataset, b_dataset, col)
# sag <- sagnage_shiny("Fil.", F, col)
# sag <- ggather(sag, vp = viewport(grobX(arr, 90), grobY(arr, 90)))

## a_dataset to b_spectral_similarity
# col <- fun_pal[2]
# b_spec <- setnullvp("b::.*::spectral_similarity", list(x = 0, y = .5), .over)
# baf2 <- baf(grobX(b_spec, 0), grobY(b_spec, 0), u(3, line), u(1, line))
# arr2 <- garrow_snake(a_dataset, b_spec, col, cur = -1)
# sag2 <- sagnage_shiny("Com.", F, col)
# sag2 <- ggather(sag2, vp = viewport(grobX(arr2, 270), grobY(arr2, 270)))

## b_nebula_index to c_child_nebula
col <- fun_pal[3]
b_index <- setnullvp("b::.*::nebula_index", list(x = 1, y = .5), .over)
b_index. <- setnullvp("b::.*::nebula_index", list(x = 1.2, y = .5), .over)
c_child <- setnullvp("c::.*::child_nebulae", list(x = 0, y = .95), .over)
baf3 <- baf(grobX(b_index, 0), grobY(b_index, 0), u(3, line), u(1, line))
arr3 <- garrow_snake(b_index., c_child, col)
sag3 <- sagnage_shiny("Cre.", F, col)
sag3 <- ggather(sag3, vp = viewport(grobX(arr3, 90), grobY(arr3, 90)))

## b ... to c_vis2
col <- "grey50"
b_spec <- setnullvp("b::.*::spectral_similarity", list(x = 1.2, y = .5), .over)
c_vis2 <- setnullvp("c::.*::vis2", list(x = 0, y = .1), .over)
c_vis2. <- setnullvp("c::.*::vis2", list(x = 0, y = .6), .over)
seg.x <- ruler(b_spec, c_vis2.)
seg.y <- ruler(c_vis2, c_vis2.)
arr4 <- parrow(, col)
vp <- viewport(grobX(seg.y, 0), grobY(seg.y, 0), grobWidth(seg.x), grobHeight(seg.y),
               just = c("right", "centre"))
arr4$vp <- vp
sag4 <- sagnage_shiny("Cus.", F, col)
sag4 <- ggather(sag4, vp = vp)

# ==========================================================================
# legend
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

le1 <- grecti2("...", tfill = pal[[ "slot" ]])
le2 <- grecti2("...", tfill = pal[[ "sub.slot" ]])
le <- frame_row(c(le1 = 1, leNull = .2, le2 = 1),
                namel(le1, leNull = nullGrob(), le2))
le.text <- frame_row(c(txt1 = 1, txtNull = .2, txt2 = 1),
                     list(txt1 = gtextp("slot", x = 0, hjust = 0),
                          txtNull = nullGrob(),
                          txt2 = gtextp("sub slot", x = 0, hjust = 0)))
le <- frame_col(c(leNull2 = .2, le = .6, leNull2 = .2, le.text = 1),
                namel(le, leNull2 = nullGrob(), le.text))

dsf_fpal <- dplyr::select(link_b, color, fun) %>% 
  rbind(., dplyr::select(link_c, color, fun)) %>% 
  rbind(.,
    # c(fun_pal[1], "filter"), c(fun_pal[2], "compute"),
    c(fun_pal[3], "create"), c("grey50", "custom")) %>% 
  dplyr::mutate(des = stringr::str_extract(fun, "^[a-zA-Z]{1,}")) %>% 
  dplyr::distinct(color, des) %>% 
  split(~ des)
le3 <- sapply(dsf_fpal, simplify = F,
  function(df) {
    n <- if (nrow(df) > 1) 2 else 1
    sags <- lapply(1:n,
      function(n) {
        text <- paste0(substr(df$des[n], 1, 3), ".")
        sag <- sagnage_shiny(text, F, df$color[n])
        arr <- garrow_snake(list(x1 = .5, y1 = 1),
          list(x2 = .5, y2 = 0),
          df$color[n])
        size <- grobHeight(gtext("text")) * 5
        ggather(arr, sag, vp = viewport(, , size, size))
      })
    if (n  == 2)
      sags <- frame_col(c(c1 = 1, c2 = 2),
        list(c1 = sags[[1]], c2 = sags[[2]]))
    else
      sags[[1]]
  })
le3.obj <- sapply(sort(names(le3)), function(name) le3[[ name ]], simplify = F)
le3 <- frame_row(fill_list(names(le3.obj), 1), le3.obj)
dup <- vapply(dsf_fpal, function(df) if (nrow(df) > 1) T else F, T)
texts <- ifelse(dup, paste0("...  ", form(names(le3.obj))), names(le3.obj))
le3.text <- sapply(texts, simplify = F,
  function(lab) gtextp(lab, x = 0, hjust = 0))
le3.text <- frame_row(fill_list(names(le3.text), 1), le3.text)
le3 <- frame_col(c(le3 = 1, le3.text = 1), namel(le3, le3.text))

legend <- frame_row(c(tit1 = .1, le = .3, legendNull = .1,
                      tit2 = .1, le3 = 1),
                    namel(tit1 = gtext("Object", list(cex = 1.2)), le,
                          tit2 = gtext("Function", list(cex = 1.2)), le3,
                          legendNull = nullGrob()))
legend <- ggather(legend, vp = viewport(.2, , , .8))

# ==========================================================================
# Altogether
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

figure <- ggather(.over, baf3, arr3, sag3, arr4, sag4)
figure <- frame_col(c(figure = 1, legend = .15),
                    namel(figure, legend))

# ==========================================================================
# output
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## gather all
pdf(tmp_pdf(), 15, 11)
draw(figure)
dev.off()
# op(tmp_pdf())
# file.copy(tmp_pdf(), "~/Documents/figure.pdf")
