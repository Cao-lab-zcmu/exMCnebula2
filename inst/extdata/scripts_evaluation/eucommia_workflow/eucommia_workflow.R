# ==========================================================================
# workflow to process data and output report (Eucomma)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# workflow(mode = "print")

s0.1 <- new_section("Abstract", 1, reportDoc$abstract, NULL)

s0.2 <- new_section("Introduction", 1,
  c(reportDoc$introduction, "",
    "*Eucommia ulmoides Oliv.* (*E.  ulmoides*), as a traditional",
    "Chinese medicine (TCM), after being processed with saline water, was",
    "applied to the treatment of renal diseases for a long time in China.",
    "Due to its complex composition, discovering chemical",
    "changes during processing (such as processed with saline water) is challenging.",
    "We would next demonstrate the addressing of this challenge with MCnebula, which",
    "may be enlightening for the study of phytopharmaceuticals."), NULL)

s0.3 <- new_section("Set-up", 1, reportDoc$setup,
  rblock({
    library(MCnebula2)
    library(exMCnebula2)
  }, F)
)

s0.9 <- new_heading("Integrate data and Create Nebulae", 1)

s1 <- new_heading("Initialize analysis", 2)

s1.1 <- new_section2(
  c("Set SIRIUS project path and its version to initialize mcnebula object."),
  rblock({
    mcn <- mcnebula()
    mcn <- initialize_mcnebula(mcn, "sirius.v4", ".")
    ion_mode(mcn) <- "neg"
  }, eval = F)
)

s1.5 <- new_section2(
  c("Create a temporary folder to store the output data."),
  rblock({
    tmp <- paste0(tempdir(), "/temp_data")
    dir.create(tmp, F)
  })
)

s1.6 <- new_section2(
  c("In order to demonstrate the process of analyzing data with MCnebula2,",
    "we provide a 'mcnebula' object that was extracted in advance using the",
    "`collate_used` function, which means that all the data",
    "used in the subsequent analysis has already stored in this 'mcnebula'",
    "object, without the need to obtain it from the original Project",
    "directory. This avoids the hassle of downloading and storing a dozen",
    "GB of raw files. The following, we",
    "use the collated dataset containing 1612 features",
    "with chemical formula identification."),
  rblock({
    exfiles <- system.file("extdata", package = "exMCnebula2")
  }))

s1.7 <- new_section2(
  c("Load the '.rdata' file."),
  rblock({
    load(paste0(exfiles, "/mcn_herbal1612.rdata"))
    mcn <- mcn_herbal1612
    rm(mcn_herbal1612)
    export_path(mcn) <- tmp
  })
)

s2 <- new_heading("Filter candidates", 2)

s2.1 <- new_section2(
  reportDoc$filter,
  rblock({
    mcn <- filter_structure(mcn)
    mcn <- create_reference(mcn)
    mcn <- filter_formula(mcn, by_reference = T)
  })
)

s3 <- new_heading("Filter chemical classes", 2)

s3.1 <- new_section2(
  reportDoc$stardust,
  rblock({
    mcn <- create_stardust_classes(mcn)
    mcn <- create_features_annotation(mcn)
    mcn <- cross_filter_stardust(mcn, max_ratio = .1, cutoff = .4, identical_factor = .6)
    classes <- unique(stardust_classes(mcn)$class.name)
    table.filtered.classes <- backtrack_stardust(mcn)
  })
)

s3.5 <- new_section2(
  c("Manually filter some repetitive classes or sub-structural classes.",
    "By means of Regex matching, we obtained a number of recurring",
    "name of chemical classes that would contain manay identical compounds",
    "as their sub-structure."),
  rblock({
    classes
    pattern <- c("fatty acid", "hydroxy")
    dis <- unlist(lapply(pattern, grep, x = classes, ignore.case = T))
    dis <- classes[dis]
    dis
  }, args = list(eval = T))
)

s3.6 <- new_section2(
  c("Remove these classes."),
  rblock({
    mcn <- backtrack_stardust(mcn, dis, remove = T)
  })
)

s4 <- new_heading("Create Nebulae", 2)

s4.1 <- new_section2(
  c("Create Nebula-Index data. This data created based on 'stardust_classes' data."),
  rblock({
    mcn <- create_nebula_index(mcn)
  })
)

s4.5 <- new_section2(
  reportDoc$nebulae,
  rblock({
    mcn <- compute_spectral_similarity(mcn)
    mcn <- create_parent_nebula(mcn)
    mcn <- create_child_nebulae(mcn)
  })
)

s5 <- new_heading("Visualize Nebulae", 2)

s5.1 <- new_section2(
  c("Create layouts for Parent-Nebula or Child-Nebulae visualizations."),
  rblock({
    mcn <- create_parent_layout(mcn)
    mcn <- create_child_layouts(mcn)
    mcn <- activate_nebulae(mcn)
  })
)

s5.3 <- new_section2(
  c("The available chemical classes for visualization and its",
    "sequence in storage."),
  rblock({
    table.nebulae <- visualize(mcn)
    table.nebulae
  }, args = list(eval = T))
)

s5.6 <- new_section2(
  c("Draw and save as .png or .pdf image files."),
  rblock({
    p <- visualize(mcn, "parent")
    ggsave(f5.61 <- paste0(tmp, "/parent_nebula.png"), p)
    pdf(f5.62 <- paste0(tmp, "/child_nebula.pdf"), 12, 14)
    visualize_all(mcn)
    dev.off()
  })
)

s5.6.fig1 <- include_figure(f5.61, "parent", "Parent-Nebula")
s5.6.fig2 <- include_figure(f5.62, "child", "Child-Nebulae")

ref <- function(x) {
  paste0("(Fig. ", get_ref(x), ")")
}

s5.8 <- c(
  "In general, Parent-Nebulae", ref(s5.6.fig1),
  "is too informative to show, so Child-Nebulae", ref(s5.6.fig2),
  "was used to dipict the abundant classes of features (metabolites)",
  "in a grid panel, intuitively. In a bird's eye view of",
  "Child-Nebulae, we can obtain many characteristics of features,",
  "involving classes distribution, structure identified accuracy, as",
  "well as spectral similarity within classes."
)

s6 <- new_heading("Nebulae for Downstream analysis", 1)

## Statistic analysis
s7 <- new_heading("Statistic analysis", 2)

s7.1 <- new_section2(
  c("Next we perform a statistical analysis with quantification data of the",
    "features. Note that the SIRIUS project does not contain quantification",
    "data of features, so our object `mcn` naturally does not contain",
    "that either. We need to get it from elsewhere."),
  rblock({
    utils::untar(paste0(exfiles, "/herbal.tar.gz"), exdir = tmp)
    origin <- data.table::fread(paste0(tmp, "/features.csv"))
    origin <- tibble::as_tibble(origin)
  })
)

s7.2 <- new_section2(
  c("Now, let's check the columns in the table."),
  rblock({
    origin
  }, args = list(eval = T))
)

s7.3 <- new_section2(
  c("Remove the rest of the columns and keep only the columns for ID,",
    "m/z, retention time, and quantification."),
  rblock({
    quant <- dplyr::select(
      origin, id = 1, dplyr::contains("Peak area")
    )
    colnames(quant) <- gsub("\\.mzML Peak area", "", colnames(quant))
    quant <- dplyr::mutate(quant, .features_id = as.character(id))
  })
)

s7.6 <- new_section2(
  c("Create the metadata table and store it in the `mcn` object",
    "along with the quantification data."),
  rblock({
    gp <- c(Blank = "EU-BlANK", Raw = "EU-Raw", Pro = "EU-Pro")
    metadata <- MCnebula2:::group_strings(colnames(quant), gp, "sample")
    metadata$annotation <-
      vapply(metadata$group, switch, FUN.VALUE = character(1),
        Blank = "methanol/water (1:1, v/v)",
        Raw = "Raw bark of Eucommia ulmoides Oliv.",
        Pro = "Precessed (with saline water) bark of Eucommia ulmoides Oliv."
      )
    features_quantification(mcn) <- dplyr::select(quant, -id)
    sample_metadata(mcn) <- metadata
  })
)

s7.7 <- new_section2(
  c(reportDoc$statistic, "", "In the following we use the",
    "`binary_comparison` function for variance analysis.",
    "To accommodate the downstream analysis of gene",
    "expression that the `limma` package was originally used for, we",
    "should log2-transform and centralize this data",
    "(use the default parameter 'fun_norm' of `binary_comparison()`)."),
  rblock({
    mcn <- binary_comparison(mcn, Pro - Raw)
    top.list <- top_table(statistic_set(mcn))
  })
)

s7.8 <- new_section2(
  c("Check the results."),
  rblock({
    top.list[[1]]
  }, args = list(eval = T, echo = T))
)

## Set tracer in Child-Nebulae
s8 <- new_heading("Set tracer in Child-Nebulae", 2)

s8.1 <- new_section2(
  reportDoc$tracer,
  rblock({
    n <- 20
    tops <- select_features(
      mcn, tani.score_cutoff = .5, order_by_coef = 1, togather = T
    )
    top20 <- tops[1:n]
    palette_set(melody(mcn)) <- colorRampPalette(palette_set(mcn))(n)
    mcn2 <- set_tracer(mcn, top20)
    mcn2 <- create_child_nebulae(mcn2)
    mcn2 <- create_child_layouts(mcn2)
    mcn2 <- activate_nebulae(mcn2)
    mcn2 <- set_nodes_color(mcn2, use_tracer = T)
  })
)

s8.2 <- new_section2(
  c("Draw and save the image."),
  rblock({
    pdf(f8.2 <- paste0(tmp, "/tracer_child_nebula.pdf"), 12, 14)
    visualize_all(mcn2)
    dev.off()
  })
)

s8.2.fig1 <- include_figure(f8.2, "tracer", "Tracing top features in Child-Nebulae")

s8.3 <- c("A part of the top features are marked with colored nodes in",
          "Child-Nebulae", paste0(ref(s8.2.fig1), "."))

## Quantification in Child-Nebulae
s9 <- new_heading("Quantification in Child-Nebulae", 2)

s9.1 <- new_section2(
  c("Show Fold Change (Pro versus Raw) in Child-Nebulae."),
  rblock({
    palette_gradient(melody(mcn2)) <- c("blue", "grey90", "red")
    mcn2 <- set_nodes_color(mcn2, "logFC", top.list[[1]])
    pdf(f9.1 <- paste0(tmp, "/logFC_child_nebula.pdf"), 12, 14)
    visualize_all(mcn2, fun_modify = modify_stat_child)
    dev.off()
  })
)

s9.1.fig1 <- include_figure(f9.1, "logFC", "Show log2(FC) in Child-Nebulae")

s9.2 <- c("Each Child-Nebula separately shows the overall content variation of",
          "the chemical class to which it belongs", ref(s9.1.fig1), ".")

## Annotate Nebulae
s10 <- new_heading("Annotate Nebulae", 2)

s10.0 <- new_section2(
  c("Now, the available Nebulae contains:"),
  rblock({
    table.nebulae2 <- visualize(mcn2)
    print(table.nebulae2, n = Inf)
  }, args = list(eval = T))
)

s10.1 <- new_section2(
  c("Next, let us focus on 'Lignans, neolignans and related compounds' and",
    "'Iridoids and derivatives'. They were representative chemical classes in",
    "E. ulmoides."),
  rblock({
    mcn2 <- set_nodes_color(mcn2, use_tracer = T)
    palette_stat(melody(mcn2)) <- c(
      Pro = "#EBA9A7", Raw = "#ACDFEE", Blank = "grey80"
    )
    lig <- "Lignans, neolignans and related compounds"
    iri <- "Iridoids and derivatives"
    mcn2 <- annotate_nebula(mcn2, lig)
    mcn2 <- annotate_nebula(mcn2, iri)
  })
)

s10.2 <- new_section2(
  c("Draw and save the image."),
  rblock({
    p <- visualize(mcn2, lig, annotate = T)
    ggsave(f10.2.1 <- paste0(tmp, "/lig_child.pdf"), p, width = 6, height = 4)
    p <- visualize(mcn2, iri, annotate = T)
    ggsave(f10.2.2 <- paste0(tmp, "/iri_child.pdf"), p, width = 6, height = 4)
  })
)

# pngGrob <- zoom_pdf(f10.2.2, dpi = 1000, position = c(.3, .5))
# grid.draw(pngGrob)

s10.2.fig1 <- include_figure(f10.2.1, "lig", paste0("Annotated Nebulae: ", lig))
s10.2.fig2 <- include_figure(f10.2.2, "iri", paste0("Annotated Nebulae: ", iri))

s10.3 <- c(
  "See results (Fig.", paste0(get_ref(s10.2.fig1), " and ", get_ref(s10.2.fig2)),
  ").", "", reportDoc$annotate
)

s10.4 <- new_section2(
  c("Use the `show_node` function to get the annotation details",
    "for a feature. For example:"),
  rblock({
    ef <- "1642"
    pdf(f10.4 <- paste0(tmp, "/features_", ef, ".pdf"), 10, 6)
    show_node(mcn2, ef)
    dev.off()
  })
)

s10.4.fig1 <- include_figure(f10.4, "ef", "The annotated feature of ID: 1642")

s10.5 <- c(
  "See results", ref(s10.4.fig1),
  paste0("(feature in ", get_ref(s10.2.fig2), ")"), "."
)

## Output identification table
s11 <- new_heading("Query compounds", 2)

s11.1 <- c("The `features_annotation(mcn)` contains the main annotation information of all",
           "the features, i.e., the identity of the  compound. Next, we would",
           "query the identified compounds based on the 'inchikey2d' column therein.",
           "Note that the stereoisomerism of the compounds is difficult to be",
           "determined due to the limitations of MS/MS spectra.",
           "Therefore, we used InChIKey 2D (representing the molecular",
           "backbone of the compound) to query",
           "the compound instead of InChI.")

s11.2 <- new_section2(
  c("First we need to format and organize the annotated data of",
    "features to get the non-duplicated 'inchikey2d'.",
    "We provide a function with a pre-defined filtering algorithm to quickly",
    "organize the table.",
    "By default, this function filters the data based on",
    "'tani.score' (Tanimoto similarity),",
    "and then sorts and de-duplicates it."),
  rblock({
    feas <- features_annotation(mcn2)
    feas <- merge(feas, top.list[[1]], by = ".features_id", all.x = T)
    feas <- dplyr::mutate(feas, arrange.rank = adj.P.Val)
    feas <- format_table(feas, export_name = NULL)
    key2d <- feas$inchikey2d
  })
)

s11.3 <- new_section2(
  c("Create a folder to store the acquired data."),
  rblock({
    tmp2 <- paste0(tmp, "/query")
    dir.create(tmp2, F)
  })
)

s11.4 <- new_section2(
  c("Query the compound's InChIKey, chemical class, IUPUA name.",
    "If your system is not Linux, the multithreading below may pose some problems,",
    "please remove the parameters `curl_cl = 4` and `classyfire_cl = 4`."),
  rblock({
    key.rdata <- query_inchikey(key2d, tmp2, curl_cl = 4)
    class.rdata <- query_classification(key2d, tmp2, classyfire_cl = 4)
    iupac.rdata <- query_iupac(key2d, tmp2, curl_cl = 4)
  })
)

s11.5 <- new_section2(
  c("We will also query for synonyms of compounds, but this is done in",
    "'CID' (PubChem's ID), so some transformation is required."),
  rblock({
    key.set <- extract_rdata_list(key.rdata)
    cid <- lapply(key.set, function(data) data$CID)
    cid <- unlist(cid, use.names = F)
    syno.rdata <- query_synonyms(cid, tmp2, curl_cl = 4)
  })
)

s11.6 <- new_section2(
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

s11.7 <- new_section2(
  c("The formatted table as following:"),
  rblock({
    feas.table
  }, args = list(eval = T))
)

## Plot spectra of top 'features'
s12 <- new_heading("Plot spectra of top 'features'", 2)

s12.1 <- new_section2(
  c("Drawing of MS/MS spectra of top 'features'."),
  rblock({
    mcn2 <- draw_structures(mcn2, .features_id = top20)
    pdf(f12.1 <- paste0(tmp, "/msms_tops_identified.pdf"), 12, 10)
    plot_msms_mirrors(mcn2, top20)
    dev.off()
  })
)

s12.2 <- include_figure(f12.1, "msmsTops", "MS/MS spectra of top features (identified)")

s12.3 <- c("See results", paste0(ref(s12.2), "."))

s12.5 <- new_section2(
  c("Plot EIC spectra of top 'features'.",
    "(Since the code below requires .mzML mass spectrometry data, these files are too",
    "large to be stored in a package, so the code below will not be run. But we have",
    "saved the result data.)"),
  rblock({
    metadata$file <- paste0(metadata$sample, ".mzML")
    data <- plot_EIC_stack(top20, metadata,
      quant.path = paste0(tmp, "/features.csv"),
      mzml.path = "/media/echo/back/thermo_mzML_0518/",
      palette = palette_stat(mcn2)
    )
    save(data, file = paste0(tmp, "/eic_data.rdata"))
  }, F)
)

s12.6 <- new_section2(
  c("Load the saved data and draw the figure."),
  rblock({
    load(paste0(tmp, "/eic_data.rdata"))
    pdf(f12.6 <- paste0(tmp, "/eic_tops_identified.pdf"), 12, 10)
    print(data$p)
    dev.off()
  })
)

s12.7 <- new_section2(
  c("Or use following to re-plot the 'ggplot' object."),
  rblock({
    data <- plot_EIC_stack(data = data, palette = palette_stat(mcn2))
  }, F, args = list(eval = F))
)

s12.8 <- include_figure(f12.6, "eic",
  "Extracted Ions Chromatograph (EIC) of top features (identified)")

s12.9 <- c("See results", paste0(ref(s12.8), "."),
  "It can be assumed that 1642, 1785, and 2321 are newly generated compounds after",
  "the processing. According to Fig. \\@ref(fig:tracer), they were belong to chemical",
  "classes of 'Iridoids and derivatives', 'Dialkyl ethers' and",
  "'Phenylpropanoids and polyketides'."
)

s13 <- new_heading("Discover more around top 'features' in Child-Nebulae", 2)

s13.1 <- new_section2(
  c("Plot the Child-Nebulae of the two chemical classes that we have not annotated",
    "yet."),
  rblock({
    dia <- "Dialkyl ethers"
    phe <- "Phenylpropanoids and polyketides"
    mcn2 <- annotate_nebula(mcn2, dia)
    mcn2 <- annotate_nebula(mcn2, phe)
    p <- visualize(mcn2, dia, annotate = T)
    ## c(.3, .4)
    ggsave(f13.1.1 <- paste0(tmp, "/dia_child.pdf"), p, width = 6, height = 4)
    p <- visualize(mcn2, phe, annotate = T)
    ## c(.5, .8)
    ggsave(f13.1.2 <- paste0(tmp, "/phe_child.pdf"), p, width = 6, height = 4)
  })
)

s13.2 <- new_section2(
  c("Draw their partial views and put them together."),
  rblock({
    grob_iri <- zoom_pdf(f10.2.2, c(.28, .515), c(.1, .1), dpi = 3000)
    grob_dia <- zoom_pdf(f13.1.1, c(.25, .3), c(.1, .1), dpi = 3000)
    grob_phe <- zoom_pdf(f13.1.2, c(.465, .825), c(.06, .06), dpi = 3000)
    local_1642 <- into(grecta("a"), grob_iri)
    local_1785 <- into(grecta("b"), grob_dia)
    local_2321 <- into(grecta("c"), grob_phe)
    locals <- frame_row(c(local_1642 = 1, local_1785 = 1, local_2321 = 1),
      namel(local_1642, local_1785, local_2321))
  })
)

s13.3 <- new_section2(
  c("We found interesting adjacent compounds (ID:2110 and ID:854) in `local_1642`,",
    "which has a similar chemical structure to 'feature' 1642."),
  rblock({
    grob_struc2110 <- grid.grabExpr(show_structure(mcn2, "2110"))
    grob_struc2110 <- into(grecta("d"), grob_struc2110)
    grob_struc854 <- grid.grabExpr(show_structure(mcn2, "854"))
    grob_struc854 <- into(grecta("e"), grob_struc854)
    grob_msms2110 <- as_grob(plot_msms_mirrors(mcn2, c("2110", "854"), structure_vp = NULL))
    grob_msms2110 <- into(grecta("f"), grob_msms2110)
  })
)

s13.4 <- new_section2(
  c("Again, we don't run the following code, but we save the results."),
  rblock({
    data <- plot_EIC_stack(c("2110", "854"), metadata,
      quant.path = paste0(tmp, "/features.csv"),
      mzml.path = "/media/echo/back/thermo_mzML_0518/",
      palette = palette_stat(mcn2)
    )
    save(data, file = paste0(tmp, "/eic_data2110.rdata"))
  }, F)
)

s13.5 <- new_section2(
  c("Load the data and convert picture."),
  rblock({
    load(paste0(tmp, "/eic_data2110.rdata"))
    grob_eic2110 <- as_grob(data$p)
    grob_eic2110 <- into(grecta("g"), grob_eic2110)
  })
)

s13.6 <- new_section2(
  c("Draw the final figure."),
  rblock({
    frame <- frame_col(c(grob_struc2110 = 1, grob_struc854 = 1),
      namel(grob_struc2110, grob_struc854))
    frame <- frame_row(c(frame = 1, grob_msms2110 = 1, grob_eic2110 = 1),
      namel(frame, grob_msms2110, grob_eic2110))
    frame <- frame_col(c(locals = 1, frame = 1.5),
      namel(locals, frame))
    frame <- ggather(frame, vp = viewport(, , .95, .95))
    pdf(f13.6 <- paste0(tmp, "/complex.pdf"), 14, 10)
    draw(frame)
    dev.off()
  })
)

s13.7 <- include_figure(f13.6, "complex",
  "Discover chemical changes using MCnebula")

s13.8 <- c("It can be speculated that the changes in the levels of ID 1642 and ID 845 were",
  "caused by structural changes of ID 2110 during the processing, which may have",
  "involved reactions such as dehydration and rearrangement (Fig. \\@ref(fig:complex)).")

s100 <- new_heading("Session infomation", 1)

s100.1 <- rblock({
  sessionInfo()
}, args = list(eval = T))

# ==========================================================================
# gather
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sections <- gather_sections()
report <- do.call(new_report, sections)
yaml(report)[1] <- c("title: Analysis on *E. ulmoides* dataset")
render_report(report, file <- paste0(tmp, "/report.Rmd"))
rmarkdown::render(file)

# ==========================================================================
# as biocStyle
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# library(MCnebula2)

# write_biocStyle(report, file2 <- paste0(tmp, "/report_biocStyle_nofloat.Rmd"),
#   title <- paste0(yaml(report)[1], "\nauthor: 'LiChuang Huang'")
# )

# rmarkdown::render(file2)
