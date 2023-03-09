# ==========================================================================
# gather the grobs as a figure
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

pal <- c(fun = "#BC3C29FF", ex = "#E18727FF")

## inner filter
inner <- grecti3("inner filter", tfill = pal[[ "fun" ]])
inner <- into(inner, p1)

## stardust classes
stardust <- grecti3("stardust classes", tfill = pal[[ "ex" ]])
stardust <- into(stardust, frame_row(fill_list(names(lst), 1), lst))

## quantity
ass <- into(glayer(4), p.peak.yellow)
sum <- into(glayer(6), p.peak.grey)
seg <- segmentsGrob(0, 0, 1, 1)
arr <- parrow(3, "black", "solid")
mp.ring <- ggather(p.ring, vp = viewport(.4, , .8))
c1 <- frame_col(c(ass = 1, seg = .4, sum = 1, arr = .3, mp.ring = 1),
                namel(ass, seg, sum, arr, mp.ring))
mp.bar <- ggather(p.bar, vp = viewport(.2, , .8))
c2 <- frame_col(c(ass = 1, arr = .3, mp.bar = 1), namel(ass, arr, mp.bar))
sep <- segmentsGrob(.5, .1, .5, .9)
## gather
content <- frame_col(c(c1 = 1, c2 = .7), namel(c1, c2))
quantity <- grectn_frame(zo(content, h = .7), gtext0("quantity  (Abundance selection)"), zo = F)

## score
struc <- into(glayer(4), struc_2095)
## gather
content <- frame_col(c(struc = 1, arr = .3, p.box = 1, arr = .3, p.ring2 = 1),
                     namel(struc, arr, p.box, p.ring2))
score <- grectn_frame(zo(content, h = .8), gtext0("score  (Goodness assessment)"), zo = F)

## identical
identical <- frame_col(c(hier_tree = .7, arr = .15, p.upset = 1.4, null = .1, venn = .3),
                       namel(hier_tree, arr, venn, null = nullGrob(), p.upset))
identical <- grectn_frame(identical, gtext0("identical  (Identicality assessment)"))

## cross
cross <- grecti3("cross filter", tfill = pal[[ "fun" ]])
obj <- sapply(namel(quantity, score, identical), zo, simplify = F)
obj <- frame_row(fill_list(names(obj), 1), obj)
cross <- into(cross, obj)
# draw(cross)

# ==========================================================================
# gather all
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

prin <- frame_col(c(inner = 1, null = .1, stardust = .7, null = .1, cross = 3.5),
                  namel(inner, stardust, cross, null = nullGrob()))
vp <- viewport(, , .95, .95)
prin <- ggather(prin, vp = vp)

pdf("figure_mech.pdf", 11, 7)
draw(prin)
dev.off()

