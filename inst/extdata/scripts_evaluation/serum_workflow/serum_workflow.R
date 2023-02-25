# ==========================================================================
# workflow to process data and output report (Serum)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Initialize
s0 <- new_heading("Initialize analysis", 2)

s0.02 <- new_section2(
  c("Set SIRIUS project path and its version to initialize mcnebula object."),
  rblock({
    mcn <- mcnebula()
    mcn <- initialize_mcnebula(mcn, "sirius.v4", ".")
    ion_mode(mcn) <- "pos"
  }, eval = F)
)

s0.1 <- new_section2(
  c("Create a temporary folder to store the output data."),
  rblock({
    tmp <- paste0(tempdir(), "/temp_data")
    dir.create(tmp, F)
    export_path(mcn) <- tmp
  })
)

s0.2 <- new_section2(
  c("In order to demonstrate the process of analyzing data with MCnebula2,",
    "we provide a 'mcnebula' object that was extracted in advance using the",
    "`collate_used` function, which means that all the data",
    "used in the subsequent analysis has already stored in this 'mcnebula'",
    "object, without the need to obtain it from the original Project",
    "directory. This avoids the hassle of downloading and storing a dozen",
    "GB of raw files. The following, we",
    "use the collated dataset containing 6501 features",
    "with chemical formula identification.",
    "This dataset was origin and processed from the research in article:",
    "<https://doi.org/10.1016/j.cell.2020.07.040>"),
  rblock({
    exfiles <- system.file("extdata", package = "exMCnebula2")
  }))

s0.3 <- new_section2(
  c("Load the '.rdata' file."),
  rblock({
    load(paste0(exfiles, "/mcn_serum6501.rdata"))
    mcn <- mcn_serum6501
    rm(mcn_serum6501)
  })
)

## Filter candidates
s1 <- new_heading("Filter candidates", 2)

s1.1 <- new_section2(
  reportDoc$filter,
  rblock({
    mcn <- filter_structure(mcn)
    mcn <- create_reference(mcn)
    mcn <- filter_formula(mcn, by_reference = T)
  })
)

## Filter chemcial classes
s2 <- new_heading("Filter chemical classes", 2)

s2.1 <- new_section2(
  reportDoc$stardust,
  rblock({
    mcn <- create_stardust_classes(mcn)
    mcn <- create_features_annotation(mcn)
    mcn <- cross_filter_stardust(mcn, max_ratio = .05, cutoff = .4, identical_factor = .6)
    classes <- unique(stardust_classes(mcn)$class.name)
    table.filtered.classes <- backtrack_stardust(mcn)
  })
)

s2.5 <- new_section2(
  c("Manually filter some repetitive classes or sub-structural classes.",
    "By means of Regex matching, we obtained a number of recurring",
    "name of chemical classes that would contain manay identical compounds",
    "as their sub-structure."),
  rblock({
    classes
    pattern <- c("stero", "fatty acid", "pyr", "hydroxy", "^orga")
    dis <- unlist(lapply(pattern, grep, x = classes, ignore.case = T))
    dis <- classes[dis]
    dis
    dis <- dis[-1]
  }, args = list(eval = T))
)

## Create Nebulae
s3 <- new_heading("Create Nebulae", 2)

s3.1 <- new_section2(
  c("Create Nebula-Index data. This data created based on 'stardust_classes' data."),
  rblock({
    mcn <- backtrack_stardust(mcn, dis, remove = T)
    mcn <- create_nebula_index(mcn)
  })
)

s3.2 <- new_section2(
  reportDoc$nebulae,
  rblock({
    mcn <- compute_spectral_similarity(mcn)
    mcn <- create_parent_nebula(mcn)
    mcn <- create_child_nebulae(mcn)
  })
)

## Visualize Nebulae
s4 <- new_heading("Visualize Nebulae", 2)

s4.1 <- new_section2(
  c("Create layouts for Parent-Nebula or Child-Nebulae visualizations."),
  rblock({
    mcn <- create_parent_layout(mcn)
    mcn <- create_child_layouts(mcn)
    mcn <- activate_nebulae(mcn)
  })
)

s4.5 <- new_section2(
  c("The available chemical classes for visualization and its",
    "sequence in storage."),
  rblock({
    table.nebulae <- visualize(mcn)
    table.nebulae
  }, args = list(eval = T))
)

s4.6 <- new_section2(
  c("Draw and save as .png or .pdf image files."),
  rblock({
    ## Parent-Nebula
    p <- visualize(mcn, "parent")
    ggsave(f4.61 <- paste0(tmp, "/parent_nebula.png"), p)
    ## Child-Nebulae
    pdf(f4.62 <- paste0(tmp, "/child_nebula.pdf"), 12, 14)
    visualize_all(mcn)
    dev.off()
  })
)

s4.6.fig1 <- include_figure(f4.61, "parent", "Parent-Nebula")
s4.6.fig2 <- include_figure(f4.62, "child", "Child-Nebulae")

ref <- function(x) {paste0("(Fig. ", get_ref(x), ")")}

s4.7 <- c("In general, Parent-Nebulae", ref(s4.6.fig1),
          "is too informative to show, so Child-Nebulae", ref(s4.6.fig2),
          "was used to dipict the abundant classes of features (metabolites)",
          "in a grid panel, intuitively. In a bird's eye view of",
          "Child-Nebulae, we can obtain many characteristics of features,",
          "involving classes distribution, structure identified accuracy, as",
          "well as spectral similarity within classes.")

## Statistic analysis
s5 <- new_heading("Statistic analysis", 2)

s5.1 <- new_section2(
  c("Next we perform a statistical analysis with quantification data of the",
    "features. Note that the SIRIUS project does not contain quantification",
    "data of features, so our object `mcn` naturally does not contain",
    "that either. We need to get it from elsewhere."),
  rblock({
    utils::untar(paste0(exfiles, "/serum.tar.gz"), exdir = tmp)
    origin <- data.table::fread(paste0(tmp, "/serum_origin_mmc3.tsv"), skip = 1)
    origin <- tibble::as_tibble(origin)
  })
)

s5.2 <- new_section2(
  c("Its original data can obtained from:",
    paste0("<https://www.cell.com/cms/10.1016/j.cell.2020.07.040/",
      "attachment/f13178d1-d1ee-4179-9d33-227a02e604f1/",
      "mmc3.xlsx>."),
    "Now, let's check the columns in the table."),
  rblock({
    origin
  }, args = list(eval = T))
)

s5.3 <- new_section2(
  c("Remove the rest of the columns and keep only the columns for ID,",
    "m/z, retention time, and quantification."),
  rblock({
    keep <- grep("^[A-Z]{2}[0-9]{1,3}$", colnames(origin))
    quant <- dplyr::select(origin, oid = 1, mz = 2, rt = 3, dplyr::all_of(keep))
  })
)

s5.4 <- new_section2(
  c("The IDs in the data `quant` are",
    "different from the IDs in the object `mcn`, so we need to align them",
    "first, according to mz and rt (they originate from the same batch of samples).",
    "In the following, we have allowed the two sets of data to be merged as", 
    "accurately as possible in the form of evaluation of score:", "",
    "- Score = (1 - rt.difference / rt.tolerance) * rt.weight +",
    "(1 - mz.defference / mz.tolerance) * mz.weight"),
  rblock({
    meta_col <- dplyr::select(features_annotation(mcn), .features_id, mz, rt.secound)
    meta_col$rt.min <- meta_col$rt.secound / 60
    merged <- align_merge(meta_col, quant, ".features_id", rt.main = "rt.min", rt.sub = "rt")
    merged <- dplyr::select(merged, -mz.main, -mz.sub, -rt.min, -rt, -rt.secound)
  })
)

s5.5 <- new_section2(
  c("Due to the differences in feature detection algorithms,",
    "some of the features inevitably do not get matched."),
  rblock({
    merged
  }, args = list(eval = T))
)

s5.6 <- new_section2(
  c("Create the metadata table and store it in the `mcn` object",
    "along with the quantification data."),
  rblock({
    gp <- c(NN = "^NN", HN = "^HN", HS = "^HS", HM = "^HM")
    metadata <- MCnebula2:::group_strings(colnames(merged), gp, "sample")
    metadata$group_name <-
      vapply(metadata$group, switch, FUN.VALUE = character(1),
        NN = "non-hospital & non-infected",
        HN = "hospital & non-infected",
        HS = "hospital & survival",
        HM = "hospital & mortality")
    metadata$supergroup <- 
      vapply(metadata$group, switch, FUN.VALUE = character(1),
        NN = "control groups",
        HN = "control groups",
        HS = "infection groups",
        HM = "infection groups")
    features_quantification(mcn) <- dplyr::select(merged, -oid)
    sample_metadata(mcn) <- metadata
  })
)

s5.7 <- new_section2(
  c(reportDoc$statistic, "", "In the following we use the",
    "`binary_comparison` function for variance analysis. Note that",
    "the quantification data in `origin` has been normalized.",
    "To accommodate the downstream analysis of gene",
    "expression that the `limma` package was originally used for, we",
    "should log2-transform and centralize this data."),
  rblock({
    mcn <- binary_comparison(mcn, (HS + HM) - (NN + HN), HM - HS,
      fun_norm = function(x) scale(log2(x), scale = F))
    top.list <- top_table(statistic_set(mcn))
  })
)

s5.8 <- new_section2(
  c("To verify the validity of the above variance analysis, the data",
    "columns were merged to obtain the IDs from the original analysis."),
  rblock({
    top.list <- lapply(top.list, merge,
      y = dplyr::select(merged, .features_id, oid),
      by = ".features_id", all.x = T, sort = F)
    top.list <- lapply(top.list, tibble::as_tibble)
  })
)

s5.81 <- new_section2(
  c("Verify with the `EFS_Rank` and `MW_Rank` column in the `origin` data.",
    "(The original authors used the two methods to rank the features.)"),
  rblock({
    origin_top50 <- dplyr::filter(origin, EFS_Rank <= 50 | MW_Rank <= 50)
    inter. <- lapply(top.list,
      function(df) {
        match <- head(df, n = 50)$oid %in% origin_top50$Unique_ID
        oid <- head(df, n = 50)$oid[match]
        list(table.match = table(match), oid = oid)
      })
    lapply(inter., function(x) x$table.match)
  }, args = list(eval = T))
)

s5.82 <- new_section2(
  c("Let's see which compounds were identified that intersected our",
    "ranking and the original ranking of features."),
  rblock({
    inter.compound <- dplyr::filter(origin, Unique_ID %in% inter.[[2]]$oid)
    table(inter.compound$Spectral_Library_Match, useNA = "if")
  }, args = list(eval = T))
)

s5.9 <- c("Interestingly, these two compounds were critical compounds",
            "in the original study.")

## Set tracer in Child-Nebulae
s6 <- new_heading("Set tracer in Child-Nebulae", 2)

s6.1 <- new_section2(
  reportDoc$tracer,
  rblock({
    n <- 50
    tops <- unique(unlist(lapply(top.list, function(df) df$.features_id[1:n])))
    palette_set(melody(mcn)) <- colorRampPalette(palette_set(mcn))(length(tops))
    mcn2 <- set_tracer(mcn, tops)
    mcn2 <- create_child_nebulae(mcn2)
    mcn2 <- create_child_layouts(mcn2)
    mcn2 <- activate_nebulae(mcn2)
    mcn2 <- set_nodes_color(mcn2, use_tracer = T)
  })
)

s6.2 <- new_section2(
  c("Draw and save the image."),
  rblock({
    pdf(f6.2 <- paste0(tmp, "/tracer_child_nebula.pdf"), 12, 14)
    visualize_all(mcn2)
    dev.off()
  })
)

s6.2.fig1 <- include_figure(f6.2, "tracer", "Tracing top features in Child-Nebulae")

s6.3 <- c("A part of the top features are marked with colored nodes in",
          "Child-Nebulae", ref(s6.2.fig1), ".",
          "These features are at least identified with chemical molecular",
          "formula. Those that are not identified, or the Nebula-Index data do",
          "not contain the chemical class to which these features belong, are",
          "absent from the Figure.")

## Quantification in Child-Nebulae
s7 <- new_heading("Quantification in Child-Nebulae", 2)

s7.1 <- new_section2(
  c("Show Fold Change (HM versus HS) in Child-Nebulae."),
  rblock({
    palette_gradient(melody(mcn2)) <- c("blue", "grey90", "red")
    mcn2 <- set_nodes_color(mcn2, "logFC", top.list[[2]])
    pdf(f7.1 <- paste0(tmp, "/logFC_child_nebula.pdf"), 12, 14)
    visualize_all(mcn2, fun_modify = modify_stat_child)
    dev.off()
  })
)

s7.1.fig1 <- include_figure(f7.1, "logFC", "Show log2(FC) in Child-Nebulae")

s7.2 <- c("Each Child-Nebula separately shows the overall content variation of",
          "the chemical class to which it belongs", ref(s7.1.fig1), ".")

## Annotate Nebulae
s8 <- new_heading("Annotate Nebulae", 2)

s8.0 <- new_section2(
  c("Now, the available Nebulae contains:"),
  rblock({
    table.nebulae2 <- visualize(mcn2)
    table.nebulae2
  }, args = list(eval = T))
)

s8.1 <- new_section2(
  c("Next, let us focus on Acyl carnitines, a class that was highlighted",
    "in the original research and also appears in Child-Nebulae, marked by",
    "plural top features (Likewise, we annotated two other chemical",
    "classes of Nebulae)."),
  rblock({
    mcn2 <- set_nodes_color(mcn2, use_tracer = T)
    palette_stat(melody(mcn2)) <- c(
      NN = "#B6DFB6", HN = "#ACDFEE", HS = "#EBA9A7", HM = "grey70"
    )
    ac <- "Acyl carnitines"
    lpc <- "Lysophosphatidylcholines"
    ba <- "Bile acids, alcohols and derivatives"
    mcn2 <- annotate_nebula(mcn2, ac)
    mcn2 <- annotate_nebula(mcn2, lpc)
    mcn2 <- annotate_nebula(mcn2, ba)
  })
)

s8.2 <- new_section2(
  c("Draw and save the image."),
  rblock({
    p <- visualize(mcn2, ac, annotate = T)
    ggsave(f8.2 <- paste0(tmp, "/ac_child.pdf"), p, height = 5)
    p <- visualize(mcn2, lpc, annotate = T)
    ggsave(f8.2.2 <- paste0(tmp, "/lpc_child.pdf"), p, height = 5)
    p <- visualize(mcn2, ba, annotate = T)
    ggsave(f8.2.3 <- paste0(tmp, "/ba_child.pdf"), p, height = 5)
  })
)

s8.2.fig1 <- include_figure(f8.2, "ac", paste0("Annotated Nebulae: ", ac))

s8.3 <- c("See results", ref(s8.2.fig1), ".", reportDoc$annotate)

s8.4 <- new_section2(
  c("Use the `show_node` function to get the annotation details",
    "for a feature. For example:"),
  rblock({
    ef <- "2068"
    pdf(f8.4 <- paste0(tmp, "/features_", ef, ".pdf"), 10, 4)
    show_node(mcn2, ef)
    dev.off()
  })
)

s8.4.fig1 <- include_figure(f8.4, "ef", "The annotated feature of ID: 2068")

s8.5 <- c("See results", ref(s8.4.fig1), ".")

notShow1 <- new_section2(
  c("Combine the images."),
  rblock({
    p.ac <- into(grecta("a"), as_grob(visualize(mcn2, ac, annotate = T)))
    p.lpc <- into(grecta("a"), as_grob(visualize(mcn2, lpc, annotate = T)))
    p.ba <- into(grecta("b"), as_grob(visualize(mcn2, ba, annotate = T)))
    p.node <- into(grecta("b"), grid::grid.grabExpr(show_node(mcn2, ef)))
    dev.off()
    grob.ac_node <- frame_row(c(p.ac = 1, p.node = .5), namel(p.ac, p.node))
    grob.ac_node <- ggather(grob.ac_node, vp = viewport(, , .95, .95))
    grob.lpc_ba <- frame_row(c(p.lpc = 1, p.ba = 1), namel(p.lpc, p.ba))
    grob.lpc_ba <- ggather(grob.lpc_ba, vp = viewport(, , .95, .95))
    pdf(paste0(tmp, "/fig.ac_node.pdf"), 9, 10)
    draw(grob.ac_node)
    dev.off()
    pdf(paste0(tmp, "/fig.lpc_ba.pdf"), 9, 13)
    draw(grob.lpc_ba)
    dev.off()
  })
)

notShow1.5 <- new_section2(
  c("Combine the image..."),
  rblock({
    ## b plot: ac
    mcn2 <- set_nodes_color(mcn2, "logFC", top.list[[2]])
    p.logfc.ac <- visualize(mcn2, ac, modify_stat_child)
    p.logfc.ac <- into(grecta("b"), as_grob(p.logfc.ac))
    ## a plot: ac
    mcn2 <- set_nodes_color(mcn2, use_tracer = T)
    p.tracer.ac <- visualize(mcn2, ac, modify_set_labs_and_unify_scale_limits)
    ## legend for a and b
    child_legend <- MCnebula2:::.get_legend(p.tracer.ac + guides(fill = "none"))
    ## a plot: ac revise
    p.tracer.ac <- visualize(mcn2, ac, modify_default_child)
    p.tracer.ac <- into(grecta("a"), as_grob(p.tracer.ac))
    ## c plot: annotated ac
    p.ac <- into(grecta("c"), as_grob(visualize(mcn2, ac, annotate = T)))
    ## d plot: node
    p.node <- into(grecta("d"), grid::grid.grabExpr(show_node(mcn2, ef)))
    ## gather theme
    vis1 <- frame_col(
      c(p.tracer.ac = 1, p.logfc.ac = 1, child_legend = .5),
      namel(p.tracer.ac, p.logfc.ac, child_legend)
    )
    vis2 <- frame_row(
      c(vis1 = 1, p.ac = 1.5, p.node = 1),
      namel(vis1, p.ac, p.node)
    )
    vis2 <- ggather(vis2, vp = viewport(, , .95, .95))
    pdf(paste0(tmp, "/fig.ac_node2.pdf"), 10, 14)
    draw(vis2)
    dev.off()
  })
)

## Output identification table
s10 <- new_heading("Query compounds", 2)

s10.1 <- c("The `features_annotation(mcn)` contains the main annotation information of all",
           "the features, i.e., the identity of the  compound. Next, we would",
           "query the identified compounds based on the 'inchikey2d' column therein.",
           "Note that the stereoisomerism of the compounds is difficult to be",
           "determined due to the limitations of MS/MS spectra.",
           "Therefore, we used InChIKey 2D (representing the molecular",
           "backbone of the compound) to query",
           "the compound instead of InChI.")

s10.2 <- new_section2(
  c("First we need to format and organize the annotated data of",
    "features to get the non-duplicated 'inchikey2d'.",
    "We provide a function with a pre-defined filtering algorithm to quickly",
    "organize the table.",
    "By default, this function filters the data based on",
    "'tani.score' (Tanimoto similarity),",
    "and then sorts and de-duplicates it."),
  rblock({
    feas <- features_annotation(mcn2)
    feas <- merge(feas, top.list[[2]], by = ".features_id", all.x = T)
    feas <- dplyr::mutate(feas, arrange.rank = adj.P.Val)
    feas <- format_table(feas, export_name = NULL)
    key2d <- feas$inchikey2d
  })
)

s10.3 <- new_section2(
  c("Create a folder to store the acquired data."),
  rblock({
    tmp2 <- paste0(tmp, "/query")
    dir.create(tmp2, F)
  })
)

s10.4 <- new_section2(
  c("Query the compound's InChIKey, chemical class, IUPUA name.",
    "If your system is not Linux, the multithreading below may pose some problems,",
    "please remove the parameters `curl_cl = 4` and `classyfire_cl = 4`."),
  rblock({
    key.rdata <- query_inchikey(key2d, tmp2, curl_cl = 4)
    class.rdata <- query_classification(key2d, tmp2, classyfire_cl = 4)
    iupac.rdata <- query_iupac(key2d, tmp2, curl_cl = 4)
  })
)

s10.5 <- new_section2(
  c("We will also query for synonyms of compounds, but this is done in",
    "'CID' (PubChem's ID), so some transformation is required."),
  rblock({
    key.set <- extract_rdata_list(key.rdata)
    cid <- lapply(key.set, function(data) data$CID)
    cid <- unlist(cid, use.names = F)
    syno.rdata <- query_synonyms(cid, tmp2, curl_cl = 4)
  })
)

s10.6 <- new_section2(
  c("Screen for unique synonyms and chemical classes for all compounds."),
  rblock({
    syno <- pick_synonym(key2d, key.rdata, syno.rdata, iupac.rdata)
    feas$synonym <- syno
    class <- pick_class(key2d, class.rdata)
    feas$class <- class
    feas.table <- rename_table(feas)
    write_tsv(feas.table, paste0(tmp, "/compounds_format.tsv"))
  })
)

s10.7 <- new_section2(
  c("The formatted table as following:"),
  rblock({
    feas.table
  }, args = list(eval = T))
)

s10.8 <- new_section2(
  c( "Filtering our results based on the ranking of 'features' (top 25 of EFS",
      "and MWU) and identification by Wozniak et al."),
  rblock({
    origin_top50 <- dplyr::select(
      origin_top50, oid = 1, EFS_Rank, MW_Rank,
      Spectral_Library_Match
    )
    feas.otop <- dplyr::select(merged, .features_id, oid)
    feas.otop <- merge(origin_top50, feas.otop, all.x = T, by = "oid")
    feas.otop <- merge(
      feas.otop, features_annotation(mcn2),
      all.x = T, by = ".features_id"
    )
    feas.otop.format <- format_table(
      feas.otop, filter = NULL,
      select = c("oid", "EFS_Rank", "MW_Rank", "Spectral_Library_Match", .select_format),
      export_name = c(oid = "# Original ID", EFS_Rank = "# EFS_Rank",
        MW_Rank = "# MW_Rank",
        Spectral_Library_Match = "# Spectral_Library_Match",
        .export_name)
    )
    write_tsv(feas.otop.format, paste0(tmp, "/oTop50_compounds_format.tsv"))
  })
)

## Pathway enrichment 
s11 <- new_heading("Pathway enrichment", 2)

s11.1 <- new_section2(
  c("A plural number of chemical classes of interest were selected and",
    "the IDs of their features were obtained. These",
    "features were filtered with the statistic analysis data",
    "(Q value < 0.05)."),
  rblock({
    focus <- c("Acyl carnitines",
      "Lysophosphatidylcholines",
      "Bile acids, alcohols and derivatives")
    focus <- select_features(mcn, focus, q.value = .05, logfc = .3, coef = 1:2)
  })
)

s11.2 <- new_section2(
  c("Get the InChIkey 2D of these features, then get the possible",
    "InChIkey, then get the CID, and finally get the KEGG ID by the",
    "CID."),
  rblock({
    focus.key2d <- maps(feas, focus, ".features_id", "inchikey2d")
    focus.cid <- lapply(focus.key2d,
      function(key2d) {
        set <- key.set[ names(key.set) %in% key2d ]
        unlist(lapply(set, function(df) df$CID),
          use.names = F)
      })
    keggids <- cid.to.kegg(unlist(focus.cid, use.names = F))
    focus.kegg <- maps(keggids, lapply(focus.cid, as.character), "Query", "KEGG")
  })
)

s11.3 <- new_section2(
  c("Let's see what we get."),
  rblock({
    focus.kegg
  }, args = list(eval = T))
)

s11.4 <- new_section2(
  c("Using package 'FELLA' for pathway enrichment analysis.",
    "The following step create a 'database' for enrichment",
    "(It is quite time consuming and can take up to 30 minutes).",
    "Subsequently, load the data."),
  rblock({
    db.dir <- init_fella(tmp, "hsa")
    db.data <- load_fella(db.dir)
  })
)

s11.5 <- new_section2(
  c("Perform enrichment."),
  rblock({
    focus.enrich <- enrich_fella(focus.kegg, db.data)
    focus.graph <- graph_fella(focus.enrich, db.data, "pagerank")
    names(focus.graph) <- names(focus.kegg)
  })
)

s11.6 <- new_section2(
  c("Some compounds are not present in the KEGG graph and are only background",
    "compounds. Let's check the enrichment results."),
  rblock({
    !vapply(focus.graph, is.null, logical(1))
  }, args = list(eval = T))
)

s11.7 <- new_section2(
  c("Visualization of enrichment results."),
  rblock({
    p2 <- plotGraph_fella(focus.graph[[2]])
    ggsave(f11.71 <- paste0(tmp, "/ba_enrich.pdf"), p2)
    p3 <- plotGraph_fella(focus.graph[[3]])
    ggsave(f11.72 <- paste0(tmp, "/lpc_enrich.pdf"), p3)
  })
)

s11.7.fig1 <- include_figure(f11.71, "ba", "Enrichment of Pagerank of BA compounds")
s11.7.fig2 <- include_figure(f11.72, "lpc", "Enrichment of Pagerank of LPC compounds")

## Heatmap analysis
s12 <- new_heading("Heatmap analysis", 2)

s12.1 <- new_section2(
  c("Since the quantitative data obtained by merging contain many missing values,",
    "they need to be processed before plotting the heat map:",
    "For each subset of data, the missing values will be filled with the average",
    "value; if the set is all missing values, they will be filled with zero."),
  rblock({
    hp.data <- handling_na(
      features_quantification(mcn2),
      metadata = sample_metadata(mcn2)
    )
  })
)

s12.2 <- new_section2(
  c("Convert wide data to long data; log transform the values; if there is a",
    "value 0, replace it with 1/10 of the minimum value of the value column."),
  rblock({
    hp.data <- log_trans(hp.data)
  })
)

s12.3 <- new_section2(
  c("Draw heat maps for three chemical classes."),
  rblock({
    names(focus) <- c("ACs", "BAs", "LPCs")
    hp.lst <- plot_heatmap(focus, hp.data, metadata, pal_group = palette_stat(mcn2))
    ggsave(f12.31 <- paste0(tmp, "/ac_heatmap.pdf"), hp.lst[[1]], width = 13, height = 4)
    ggsave(f12.32 <- paste0(tmp, "/ba_heatmap.pdf"), hp.lst[[2]], width = 13, height = 7)
    ggsave(f12.33 <- paste0(tmp, "/lpc_heatmap.pdf"), hp.lst[[3]], width = 13, height = 4)
  })
)

# test <- plot_heatmap(focus, hp.data, metadata, pal_group = palette_stat(mcn2), clust_col = F)
# test <- lapply(test, function(x) grid.grabExpr(print(x)))
# names(test)
# # [1] "ACs"  "BAs"  "LPCs"
# test.p <- frame_row(c(ACs = 1, BAs = 1.5, LPCs = 1), test)
# pdf(paste0(tmp, "/test_heatmap.pdf"), width = 13, height = 18)
# draw(test.p)
# dev.off()

s12.4.fig1 <- include_figure(f12.31, "acHp", "Heatmap of ACs")
s12.4.fig2 <- include_figure(f12.32, "baHp", "Heatmap of BAs")
s12.4.fig3 <- include_figure(f12.33, "lpcHp", "Heatmap of LPCs")

s12.5 <- c("As shown in figure,",
  paste0("Fig. ", get_ref(s12.4.fig1), ", ", get_ref(s12.4.fig2), ", and ",
    get_ref(s12.4.fig3)),
  "ACs and BAs implied a high correlation with disease,",
  "while LPCs showed a relatively weak correlation."
)

notShow2 <- new_section2(
  c("Combined heatmap."),
  rblock({
    grobs.hp <- lapply(hp.lst, function(p) grid::grid.grabExpr(print(p)))
    grobs.hp <- lapply(1:3,
      function(n, x = c("a", "b", "c")) {
        into(grecta(x[n]), grobs.hp[[ n ]])
      })
    names(grobs.hp) <- names(hp.lst)
    grobs.hp <- frame_row(c(ACs = 1, BAs = 1.5, LPCs = 1), grobs.hp)
    grobs.hp <- ggather(grobs.hp, vp = viewport(, , .95, .95))
    pdf(paste0(tmp, "/hps.pdf"), 13, 15)
    draw(grobs.hp)
    dev.off()
  })
)

s15 <- new_heading("Verify Identification", 1)

s15.1 <- new_section2(
  c("In research of Wozniak et al, a subset of ACs compounds were identified. In",
    "addition, four top rank metabolites were identified. However, Most of the",
    "top 25 metabolites (EFS rank) were not identified."),
  rblock({
    ac_names <- c("Palmitoyl-carnitine", "Octanoyl-carnitine",
      "Acetyl-carnitine", "Hexanoyl-carnitine",
      "Decanoyl-carnitine")
    ac_inchikey2d <- c("XOMRRQXKHMYMOC", "CXTATJFJDMJMIY",
      "RDHQFKQIGNGIED", "VVPRQWTYSNDTEA",
      "LZOSYCMHQXPBFU")
    names(ac_inchikey2d) <- ac_names
    other_4m <- c(
      "Hepc", "D-erythro-Sphingosine-1-phosphate", "L-THYROXINE",
      "Decanoyl-L-carnitine"
    )
    other_4m_inchikey2d <- c(
      "BDPQVGIMLZYZQA", "DUYSYHSSBDVJSM", "XUIIKFGFIJCVMT",
      "LZOSYCMHQXPBFU"
    )
    other_4m.seq <- unlist(lapply(other_4m, grep,
        x = origin$Spectral_Library_Match, ignore.case = T))
    other_4m <- origin[other_4m.seq, ]
    origin.top25 <- c(349, 746, 854, 228, 320, 971, 2532, 670, 92, 1363,
      13, 798, 1379, 1947, 4146, 736, 1656, 464, 731, 289,
      4431, 3865, 476, other_4m$Unique_ID)
    origin.top25.featureID <-
      dplyr::filter(merged, oid %in% origin.top25)$.features_id
  })
)

s15.2 <- new_section2(
  c("Confirm which compounds were identified:"),
  rblock({
    ac_inchikey2d %in% features_annotation(mcn2)$inchikey2d
    other_4m_inchikey2d %in% features_annotation(mcn2)$inchikey2d
    reIdentify.origin.top25 <- dplyr::filter(features_annotation(mcn2),
      .features_id %in% origin.top25.featureID)
    reIdentify.origin.top25 <- dplyr::select(
      reIdentify.origin.top25, .features_id, tani.score, inchikey2d
    )
    print(reIdentify.origin.top25, n = Inf)
  }, args = list(eval = T))
)

## Session infomation
s100 <- new_heading("Session infomation", 1)

s100.1 <- rblock({
  sessionInfo()
}, args = list(eval = T))

# ==========================================================================
# output report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sections <- gather_sections()
prologue <- list(
  new_section("Abstract", 1, reportDoc$abstract, NULL),
  new_section("Introduction", 1, reportDoc$introduction, NULL),
  new_section("Set-up", 1, reportDoc$setup,
    rblock({
      library(MCnebula2)
      library(exMCnebula2)
    }, F)
  )
) 
report <- do.call(new_report, c(prologue, sections))
yaml(report)[1] <- c("title: Analysis on serum dataset")

## post-modify, add heading
h1 <- new_heading("Integrate data and Create Nebulae", 1)
h2 <- new_heading("Nebulae for Downstream analysis", 1)
seqs <- search_heading(report, "^Initialize|^Statistic")
report <- insert_layers(report, seqs[1], h1)
report <- insert_layers(report, seqs[2], h2)

render_report(report, file.report <- paste0(tmp, "/report.rmd"))
rmarkdown::render(file.report)
