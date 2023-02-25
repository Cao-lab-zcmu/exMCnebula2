# ==========================================================================
# script for evaluation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

s0.1 <- new_heading("Introduction", 1)

s0.2 <- c(
  "This document provides the code to evaluate MCnebula2. Most of the data used",
  "in this is already included in the package 'exMCnebula2'. By downloading and",
  "installing the 'exMCnebula2' package, the following code can be run without",
  "trouble. The code in this document can be divided into two parts:",
  "",
  "- The first part is the code used to generate the predicate dataset (for",
  "evaluation), which is time consuming to run, but we have already run it in",
  "advance and included the outcome data in the package 'exMCnebula2', so users",
  "do not need to rerun them, unless they are tested for feasibility.",
  "",
  "- The second part of the code is the code for evaluating the results, including",
  "the data processing and data visualization modules, which is lightweight and",
  "can be run quickly to get the results."
)

s0.3 <- new_heading("Set-up", 1)

s0.4 <- new_section2(
  c("Load the R package used for analysis.  In the following analysis process, to",
    "illustrate the source of the function, we use the symbol `::` to mark the",
    "functions, e.g., `dplyr::filter`.  The functions that were not marked may",
    "source from MCnebula2 or the packages that R (version 4.2) loaded by default."),
  rblock({
    library(MCnebula2)
    library(exMCnebula2)
  }, F)
)

s0.5 <- new_heading("Initialisation", 1)

s0.9 <- new_section2(
  c("Create a temporary folder to store the output data."),
  rblock({
    tmp <- paste0(tempdir(), "/temp_data")
    dir.create(tmp, F)
    exfiles <- system.file("extdata", package = "exMCnebula2")
    utils::untar(paste0(exfiles, "/evaluation.tar.gz"), exdir = tmp)
    source(paste0(tmp, "/evaluation.R"))
  })
)

s1 <- new_heading("Use simulated dataset", 1)

s1.1 <- new_heading("Convert .msp as .mgf", 2)

s1.2 <- new_section2(
  c("We downloaded a collated spectral data (.msp file) originating from GNPS from a",
    "third-party website.", "(<http://prime.psc.riken.jp/compms/msdial/main.html#MSP>)",
    "The .msp was converted as .mgf file and, at the same time, the simulated",
    "isotope peaks were added to MS1. Then the .mgf file could be used for",
    "SIRIUS computation. The following would not be evaluated, as the converted",
    ".mgf file has been included in .tar file."),
  rblock({
    if (!requireNamespace("rcdk", quietly = T))
      install.packages("rcdk")
    msp_to_mgf(
      name = "origin_gnps_pos.msp",
      id_prefix = "gnps",
      path = tmp,
      mass_level = "all",
      fun = "deal_with_msp_record")
  }, F)
)

s1.3 <- new_heading("Query classification for compounds", 2)

s1.4 <- new_section2(
  c("In order for this data to be used for evaluation, we need to first obtain the",
    "classes of these compounds. (The results has been included in .tar files.",
    "So the following need not to be run again, as its time-consuming.).",
    "(Please note that the following cl (`curl_cl`, `classyfire_cl`) sets the number of threads,",
    "if it is not a Unix like system, there may be an error, please set this parameter to",
    "`NULL`, or refer to <https://github.com/psolymos/pbapply>.)",
    ),
  rblock({
    mgf_metadata <- data.table::fread(
      paste0(tmp, "/", "origin_gnps_pos.meta.tsv")
    )
    # mgf_metadata <- data.table::fread(
    # paste0("/media/echo/back/test_mcnebula/mgf/used", "/", "origin_gnps_pos.meta.tsv")
    # )
    dir.create(tmp1 <- paste0(tmp, "/query"))
    key2d <- unique(stringr::str_extract(mgf_metadata$INCHIKEY, "^[A-Z]*"))
    key.rdata <- query_inchikey(key2d, tmp1, curl_cl = 10)
    class.rdata <- query_classification(key2d, tmp1, classyfire_cl = 20)
  }, F)
)

s1.5 <- new_heading("Add noise peak", 2)

s1.6 <- new_section2(
  c("The 'noise' include mass shift, peak intensity shift, and external",
    "noise peaks for original spectra. (refer to <https://doi.org/10.1038/s41587-021-01045-9>).",
    "First, we collected a subset of the noise peaks used to insert the spectra:"),
  rblock({
    load(files[1])
    non_noise <- dplyr::select(
      latest(mcn, "project_dataset", ".f3_spectra"),
      .features_id, mz, rel.int.
    )
    non_noise <- split(non_noise, ~ .features_id)
    non_noise <- lapply(non_noise, dplyr::select, -.features_id)
    rm(mcn)
    origin_lst <- filter_mgf(NULL, paste0(tmp, "/origin_gnps_pos.mgf"))
    # origin_lst <- filter_mgf(NULL, "/media/echo/back/test_mcnebula/mgf/used/origin_gnps_pos.mgf")
    noise_pool <- collate_as_noise_pool(origin_lst, non_noise)
  }, F)
)

s1.7 <- new_section2(
  c("Two noise models were simulated:"),
  rblock({
    medium_noise_lst <- spectrum_add_noise(
      origin_lst,
      int.sigma = 1,
      global.sigma = 10/3 * 1e-6,
      indivi.sigma = 10/3 * 1e-6,
      sub.factor = 0.03,
      alpha = 0.2,
      .noise_pool = noise_pool
    )
    high_noise_lst <- spectrum_add_noise(
      origin_lst,
      int.sigma = 2^(1/2),
      global.sigma = 15/3 * 1e-6,
      indivi.sigma = 15/3 * 1e-6,
      sub.factor = 0.03,
      alpha = 0.4,
      .noise_pool = noise_pool
    )
  }, F)
)

s1.8 <- new_heading("Output .mgf for SIRIUS", 2)

s1.85 <- new_section2(
  c("Export the above three list data to the .mgf format required by SIRIUS."),
  rblock({
    lst <- list(
      origin = origin_lst,
      medium_noise = medium_noise_lst, 
      high_noise = high_noise_lst
    )
    dir.create(tmp2 <- paste0(tmp, "/simutate"))
    lapply(names(lst),
      function(name) {
        data <- data.table::rbindlist(lst[[ name ]])
        write.table(
          data, paste0(tmp2, "/", name, "_gnps_pos.sirius.mgf"),
          quote = F, col.names = F, row.names = F
        )
      })
  }, F)
)

s1.9 <- new_heading("Output .mgf for GNPS", 2)

s1.95 <- new_section2(
  c("Export the above three list data to the .mgf format required by GNPS (FBMN).",
    "(<https://doi.org/10.1038/s41592-020-0933-6>).",
    "FBMN worklow required two types of data for upload, .mgf file and",
    ".csv (quantification table)"),
  rblock({
    quant_table <- simulate_gnps_quant(
      mgf_metadata, tmp2, return_df = T
    )
    lapply(names(lst), 
      function(name) {
        lst <- discard_level1(lst[[ name ]])
        lst <- pbapply::pblapply(lst, mgf_add_anno.gnps)
        data <- data.table::rbindlist(lst)
        data <- dplyr::mutate(
          data, V1 = gsub("RTINSECONDS=", "RTINSECONDS=1000", V1),
          V1 = gsub("CHARGE=+1", "CHARGE=1+", V1),
          V1 = gsub("FEATURE_ID=gnps", "FEATURE_ID=", V1)
        )
        write.table(
          data, paste0(tmp2, "/", name, "_gnps_pos.gnps.mgf"),
          quote = F, col.names = F, row.names = F
        )
        quant <- dplyr::filter(
          quant_table, `row m/z` <= 800, `row ID` %in% names(lst)
        )
        quant$`row ID` <- stringr::str_extract(quant$`row ID`, "[0-9]{1,}$")
        write.table(
          quant, paste0(tmp2, "/", name, "_gnps_pos.gnps.meta.csv"),
          sep = ",", row.names = F, col.names = T, quote = F)
      })
  }, F)
)

s2 <- new_heading("Evaluate MCnebula2", 1)

s2.05 <- new_heading("Use pre-computed data", 2)

s2.05 <- new_section2(
  c("Extract the required data from the pre-computed SIRIUS projects, using the",
    "following code:"),
  rblock({
    dirs <- c("origin", "medium_noise", "high_noise")
    lst <- lapply(dirs,
      function(dir){
        mcn <- collated_used()
        mcn <- mcnebula()
        mcn <- initialize_mcnebula(mcn, "sirius.v4", dir)
        mcn <- collate_used(mcn)
      })
  }, F)
)

s2.1 <- new_heading("Integrate data", 2)

s2.11 <- new_section2(
  c("Integrate data to classifiy the 'features' via MCnebula2.",
    "(The following .rdata files were not provided in packages of 'exMCnebula2'.",
    "But you can downloaded them via URL such as:",
    "<https://raw.githubusercontent.com/Cao-lab-zcmu/utils_tool/master/inst/extdata/evaluationLarge/origin_gnps_pos.rdata>)"
    ),
  rblock({
    files <- paste0(
      system.file("extdata/evaluationLarge", package = "utils.tool"), "/",
      c("origin_gnps_pos.rdata", "medium_noise_gnps_pos.rdata", "high_noise_gnps_pos.rdata")
    )
    lst <- lapply(files,
      function(file) {
        load(file)
        mcn <- filter_structure(mcn)
        mcn <- create_reference(mcn)
        mcn <- filter_formula(mcn, by_reference = T)
        mcn <- create_stardust_classes(mcn)
        mcn <- create_features_annotation(mcn)
        mcn <- cross_filter_stardust(
          mcn, min_number = 50, max_ratio = .1, cutoff = .4, identical_factor = .6
        )
        mcn <- create_nebula_index(mcn, force = T)
        if (file == files[1]) {
          mcn <- create_hierarchy(mcn)
          mcnHier. <- mcnebula()
          reference(mcn_dataset(mcnHier.))$hierarchy <- hierarchy(mcn)
          save(mcnHier., file = paste0(system.file("extdata/evaluation",
                package = "utils.tool"), "/mcnHier.rdata"))
        }
        list(features_annotation = features_annotation(mcn),
          nebula_index = nebula_index(mcn))
      })
    names(lst) <- dirs
    save(lst, file = paste0(system.file("extdata/evaluation",
          package = "utils.tool"),
        "/integrated.rdata"))
  }, F)
)

s2.2 <- new_heading("Use pre-integrate data", 2)

s2.21 <- new_section2(
  c("Load the integrated data in the 'exMCnebula2' package. This data contains",
    "the results of the three levels, as well as the corresponding classified",
    "results and identification results."),
  rblock({
    load(paste0(tmp, "/integrated.rdata"))
    # load(paste0(exfiles, "/evaluation/integrated.rdata"))
  })
)

s2.3 <- new_heading("Download results of finished jobs from GNPS service", 2)

s2.31 <- new_section2(
  c("In order to compare MCnebula2 with MolNetEnhancer in GNPS, we have",
    "pre-uploaded and completed jobs of FBMN and MolNetEnhancer.",
    "The URL of finished jobs in GNPS service were as following.",
    "(MolNetEnhancer: <https://doi.org/10.1101/654459>)"),
  rblock({
    fbmn_lst <- list(
      origin =
        "https://gnps.ucsd.edu/proteosafe/status.jsp?task=05f492249df5413ba72a1def76ca973d",
      medium_noise = 
        "https://gnps.ucsd.edu/proteosafe/status.jsp?task=c65abe76cd9846c99f1ae47ddbd34927",
      high_noise =
        "https://gnps.ucsd.edu/proteosafe/status.jsp?task=62b25cf2dcf041d3a8b5593fdbf5ac5e"
    )
    molnet_lst <- list(
      origin =
        "https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=9d9c7f83fa2046c2bf615a3dbe35ca62",
      medium_noise = 
        "https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=7cc8b5a2476f4d4e90256ec0a0f94ca7",
      high_noise =
        "https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=f6d08a335e814c5eac7c97598b26fb80"
    )
  }, F)
)

s2.32 <- new_section2(
  c("The data that would be used has been extracted."),
  rblock({
    molnet_files <- paste0(tmp, "/", "gnps_results")
    # molnet_files <- paste0(exfiles, "/evaluation/gnps_results")
    molnet_lst <- sapply(c("origin", "medium_noise", "high_noise"), simplify = F,
      function(dir) {
        data.table::fread(paste0(molnet_files, "/", dir, "/ClassyFireResults_Network.txt"))
      })
  })
)

s2.4 <- new_heading("Evaluate accuracy of classify", 2)

s2.5 <- c("Evaluate the 'features' of each chemical class within the Nebula-Index:",
  "whether the compounds relative to the chemical class.")

s2.51 <- new_heading("Load assessment data and reference data.", 3)

s2.511 <- new_section2(
  c("The evaluation and reference data have already been saved in the previous",
    "steps and only need to be loaded:"),
  rblock({
    mgf_metadata <- data.table::fread(
      paste0(tmp, "/", "origin_gnps_pos.meta.tsv")
    )
    id2key <- stringr::str_extract(mgf_metadata$INCHIKEY, "^[A-Z]*")
    names(id2key) <- mgf_metadata$.id
    id2key <- as.list(id2key)
    load(paste0(tmp, "/mcnHier.rdata"))
    class.db <- extract_rdata_list(paste0(tmp, "/classification.rdata"))
  })
)

s2.52 <- new_heading("Filter assessment data.", 3)

s2.521 <- new_section2(
  c("Filtered according to: the common 'features' (at least identified to the",
    "chemical formula) in the three levels; whether included in reference data (`class.db`)."),
  rblock({
    common <- lapply(lst, function(l)  l$features_annotation$.features_id )
    common <- common[[1]][
      (common[[1]] %in% common[[2]]) &
      (common[[1]] %in% common[[3]]) 
    ]
    keep <- id2key %in% names(class.db)
    common <- common[ common %in% names(id2key)[keep] ]
    ## length(common) # [1] 7524
    lst <- lapply(lst,
      function(l) {
        l$nebula_index <- dplyr::filter(
          l$nebula_index, .features_id %in% common
        )
        return(l)
      })
  })
)

s2.53 <- new_heading("Count the results.", 3)

s2.531 <- new_section2(
  c("Assess the accuracy and derive ratios for each type of result."),
  rblock({
    res <- lapply(lst,
      function(l) {
        l <- lapply(split(l$nebula_index, ~ class.name),
          function(df) df[[ ".features_id" ]])
        stat_list <- sapply(names(l), simplify = F,
          function(class.name) {
            stat_classify(
              l[[ class.name ]], class.name,
              id2key, mcnHier., class.db
            )
          })
        stat_table <- data.table::rbindlist(
          lapply(stat_list, table_app, prop = T), fill = T, idcol = T
        )
        stat_table <- dplyr::rename(stat_table, class.name = .id)
        stat_table <- dplyr::summarise_all(
          stat_table, function(x) ifelse(is.na(x), 0, x)
        )
      })
  })
)

s2.54 <- new_heading("Post-filtering.", 3)

s2.541 <- new_section2(
  c("Filter out chemical classes with insufficient number of features; keep only",
    "those contained in 'origin'."),
  rblock({
    res[[1]] <- dplyr::filter(res[[1]], sum >= 50)
    keep <- res[[1]]$class.name
    res[2:3] <- lapply(res[2:3], dplyr::filter, class.name %in% keep)
  })
)

s2.55 <- new_heading("Distinguish the chemical class of the dominant structure.", 3)

s2.551 <- new_section2(
  c("Our reference data (`class.db`) contains only chemical classes for the",
    "dominant structure of a compound, but not for the sub-structure of a",
    "compound.  Our evaluation data (`lst`) contain both chemical classes of",
    "substructures and classes of dominant structures. This makes the assessment",
    "biased. But they are easily distinguishable:"),
  rblock({
    dominant_res <- lapply(res, dplyr::filter, false < .4)
    sub_res <- lapply(res, dplyr::filter, false >= .4)
    summarise <- function(df) {
      false <- round(mean(df$false) * 100, 1)
      sum <- round(mean(df$sum), 1)
      list(false = false, sum = sum)
    }
    summary <- lapply(dominant_res, summarise)
  })
)

s2.56 <- new_section2(
  c("These chemical classes represent mostly local structures in compounds (they",
    "represent very small structures that are present in many other classes of",
    "compounds):"),
  rblock({
    sub_res[[1]]$class.name
  }, args = list(eval = T))
)

s2.57 <- new_heading("Visualizaion", 3)

s2.571 <- new_section2(
  c("Visualize the results of evaluation:"),
  rblock({
    vis_lst <- lapply(dominant_res, visualize_stat, mcn = mcnHier.)
    pdfs <- list()
    for (i in names(vis_lst)) {
      pdfs[[ i ]] <- paste0(tmp, "/", i, "_classified_accuracy.pdf")
      pdf(pdfs[[ i ]], 13, 11)
      draw(vis_lst[[ i ]])
      dev.off()
    }
  })
)

s2.572 <- include_figure(pdfs[[ 1 ]], "origin",
  "Classified accuracy (MCnebula2) of origin dataset")
s2.573 <- include_figure(pdfs[[ 2 ]], "medium",
  "Classified accuracy (MCnebula2) of medium noise dataset")
s2.574 <- include_figure(pdfs[[ 3 ]], "high",
  "Classified accuracy (MCnebula2) of high noise dataset")

ref <- function(x) {
  paste0("(Fig. ", get_ref(x), ")")
}
s2.575 <- c(
  paste0("Totally ", length(common), " compounds used for evaluation."),
  "For the origin dataset, the average false rate of MCnebula2 classifying is",
  paste0(summary$origin$false, "% ", ref(s2.572), ";"),
  "the average classified number of 'features' is",
  paste0(summary$origin$sum, "."),
  "For the medium noise dataset, the average false rate of MCnebula2 classifying is",
  paste0(summary$medium_noise$false, "% ", ref(s2.573), ";"),
  "the average classified number of 'features' is",
  paste0(summary$medium_noise$sum, "."),
  "For the high noise dataset, the average false rate of MCnebula2 classifying is",
  paste0(summary$high_noise$false, "% ", ref(s2.574), "."),
  "the average classified number of 'features' is",
  paste0(summary$high_noise$sum, ".")
)

s2.576 <- new_section2(
  c("See following:"),
  rblock({
    summary
  }, args = list(eval = T, echo = T))
)

s2.577 <- new_section2(
  c("Gather three levels (origin, medium_noise, high_noise) of results:"),
  rblock({
    vis <- visualize_statComplex(
      dominant_res, mcnHier., weight = c(pl = .7, pm = 1.1, pr = .8)
    )
    pdf(f2.577 <- paste0(tmp, "/gather_classified_accuracy.pdf"), 15.5, 13)
    draw(vis)
    dev.off()
  })
)

s2.578 <- include_figure(f2.577, "gather",
  "Classified accuracy (MCnebula2) of three levels dataset")

s2.579 <- c("See results", paste0(ref(s2.578), "."))

s2.58 <- new_heading("Compare with MolNetEnhancer", 3)

s2.581 <- new_section2(
  c("Pre-process the data."),
  rblock({
    molnet_lst <- lapply(molnet_lst,
      function(df) {
        df <- dplyr::select(df, .id = `cluster index`, ends_with("class"))
        df <- dplyr::mutate(df, .id = paste0("gnps", .id))
        df <- dplyr::filter(df, .id %in% !!common)
        lst <- lapply(2:4,
          function(n) {
            lst <- split(df, df[[ n ]])
            lst <- lst[!names(lst) %in% c("", "no matches")]
            lst <- lapply(lst,
              function(df) {
                if (length(df[[ ".id" ]]) >= 50) df[[ ".id" ]]
                else NULL
              })
          })
        lst <- unlist(lst, recursive = F)
        lst[!vapply(lst, is.null, logical(1))]
      })
  })
)

s2.582 <- new_section2(
  c("Count the results."),
  rblock({
    res_molnet <- lapply(molnet_lst,
      function(l) {
        stat_list <- sapply(names(l), simplify = F,
          function(class.name) {
            stat_classify(
              l[[ class.name ]], class.name,
              id2key, mcnHier., class.db
            )
          })
        stat_table <- data.table::rbindlist(
          lapply(stat_list, table_app, prop = T), fill = T, idcol = T
        )
        stat_table <- dplyr::rename(stat_table, class.name = .id)
        stat_table <- dplyr::summarise_all(
          stat_table, function(x) ifelse(is.na(x), 0, x)
        )
      })
    res_molnet[[1]] <- dplyr::filter(res_molnet[[1]], sum >= 50)
    keep <- res_molnet[[1]]$class.name
    res_molnet[2:3] <- lapply(res_molnet[2:3], dplyr::filter, class.name %in% keep)
    summary_molnet <- lapply(res_molnet, summarise)
  })
)

s2.583 <- new_section2(
  c("Visualize the results."),
  rblock({
    vis_lst <- lapply(
      res_molnet, visualize_stat, mcn = mcnHier.,
      weight = c(pl = 1, pm = .8, pr = .7)
    )
    pdfs2 <- list()
    for (i in names(vis_lst)) {
      pdfs2[[ i ]] <- paste0(tmp, "/", i, "_classified_accuracy_Molnet.pdf")
      w <- if (i == "origin") 15 else 11
      pdf(pdfs2[[ i ]], w, 11)
      draw(vis_lst[[ i ]])
      dev.off()
    }
    vis <- visualize_statComplex(
      res_molnet, mcnHier.,
      y_cut_left = c(50, 700),
      y_cut_right = c(800, 2000),
      y_cut_left_breaks = c(50, seq(100, 700, by = 200)),
      y_cut_right_breaks = c(1200, 1600),
    )
    pdf(f2.583 <- paste0(tmp, "/gather_classified_accuracy_Molnet.pdf"), 18, 13)
    draw(vis)
    dev.off()
  })
)

s2.584 <- include_figure(pdfs2[[ 1 ]], "originMolnet",
  "Classified accuracy (MolnetEnhancer) of origin dataset")
s2.585 <- include_figure(pdfs2[[ 2 ]], "mediumMolnet",
  "Classified accuracy (MolNetEnhancer) of medium noise dataset")
s2.586 <- include_figure(pdfs2[[ 3 ]], "highMolnet",
  "Classified accuracy (MolNetEnhancer) of high noise dataset")

s2.587 <- c(
  "For the origin dataset, the average false rate of MolNetEnhancer classifying is",
  paste0(summary_molnet$origin$false, "% ", ref(s2.584), ";"),
  "the average classified number of 'features' is",
  paste0(summary_molnet$origin$sum, "."),
  "For the medium noise dataset, the average false rate of MolNetEnhancer classifying is",
  paste0(summary_molnet$medium_noise$false, "% ", ref(s2.585), ";"),
  "the average classified number of 'features' is",
  paste0(summary_molnet$medium_noise$sum, "."),
  "For the high noise dataset, the average false rate of MolNetEnhancer classifying is",
  paste0(summary_molnet$high_noise$false, "% ", ref(s2.586), ";"),
  "the average classified number of 'features' is",
  paste0(summary_molnet$high_noise$sum, ".")
)

s2.588 <- new_section2(
  c("See following:"),
  rblock({
    summary_molnet
  }, args = list(eval = T, echo = T))
)

s2.589 <- include_figure(f2.583, "gatherMol",
  "Classified accuracy (MolNetEnhancer) of three levels dataset")

s2.5891 <- c("See results", paste0(ref(s2.589), "."))

s2.59 <- new_section2(
  c("Compare classified number of features for two methods in three levels:"),
  rblock({
    vis <- visualize_comparison(dominant_res, res_molnet)
    pdf(f2.8 <- paste0(tmp, "/classified_number_comparison.pdf"))
    draw(vis)
    dev.off()
  })
)

s2.591 <- include_figure(f2.8, "numberComp",
  "Comparison of classified number for MCnebula2 and MolNetEnhancer"
)

s2.592 <- c("As shown", paste0(ref(s2.591), ","),
  "MCnebula2 has a higher noise tolerance than MolNetEnhancer.",
)

s2.59 <- new_heading("Summary", 3)

s2.591 <- new_section2(
  c("The following data is available."),
  rblock({
    res_parallel <- attr(vis, "data")
    comman_class <- unique(res_parallel[[1]]$class.name)
    res_mcnebula <- dominant_res
    res_molnet
    summary_mcnebula <- summary
    summary_molnet
  })
)

s2.592 <- new_section2(
  c("Also required:"),
  rblock({
    res_mcnebula_common <- lapply(
      res_mcnebula, dplyr::filter, class.name %in% comman_class
    )
    res_molnet_common <- lapply(
      res_molnet, dplyr::filter, class.name %in% comman_class
    )
    summary_mcnebula_common <- lapply(res_mcnebula_common, summarise)
    summary_molnet_common <- lapply(res_molnet_common, summarise)
  })
)

s2.593 <- new_section2(
  c("Draw the figure:"),
  rblock({
    lst_molnet <- visualize_summary(summary_molnet_common)
    lst_mcnebula <- visualize_summary(summary_mcnebula_common)
    vis <- lapply(namel(lst_mcnebula, lst_molnet),
      function(lst) {
        bar <- as_grob(lst$p.num)
        ring <- list(
          p.ratioRelFal = as_grob(lst$p.ratioRelFal$p.m),
          p.ratioSt = as_grob(lst$p.ratioSt$p.m)
        )
        frame <- frame_col(c(p.ratioSt = 2, p.ratioRelFal = 3), ring)
        frame_col(c(bar = 2, frame = 5), c(namel(bar), namel(frame)))
      })
    vis <- frame_row(c(lst_mcnebula = 1, lst_molnet = 1), vis)
    label_mcnebula <- gtext90("MCnebula", "#4DBBD5FF")
    label_molnet <- gtext90("GNPS", "#E64B35FF")
    group_labels <- frame_row(c(label_mcnebula = 1, null = .1, label_molnet = 1),
      namel(label_mcnebula, label_molnet, null = nullGrob()))
    vis2 <- frame_col(c(group_labels = .1, vis = 5), namel(group_labels, vis))
    legend <- frame_col(c(null = 2, p.ratioSt = 2, p.ratioRelFal = 3),
      list(null = nullGrob(), p.ratioRelFal = lst_molnet$p.ratioRelFal$p.l,
        p.ratioSt = lst_molnet$p.ratioSt$p.l))
    vis3 <- frame_row(c(vis2 = 5, legend = .5), namel(vis2, legend))
    label_sum <- gtext90("Sum number", "#709AE1FF", 0)
    label_false <- gtext90("Relative false rate", "#FED439FF", 0)
    label_stab <- gtext90("Stability", "#91D1C2", 0)
    title_labels <- frame_col(
      c(null = .2, label_sum = 2, null = .1, label_stab = 2, null = .1, label_false = 3), 
      namel(label_sum, label_false, label_stab, null = nullGrob())
    )
    vis4 <- frame_row(
      c(title_labels = 1, null = .5, vis3 = 15), 
      namel(title_labels, vis3, null = nullGrob())
    )
    vis4 <- ggather(vis4, vp = viewport(, , .95, .95))
    pdf(f2.59 <- paste0(tmp, "/evaluation_summary.pdf"), 15, 5)
    draw(vis4)
    dev.off()
  })
)

s2.594 <- include_figure(f2.59, "summary", "Evaluation sumary for MCnebula and benchmark method")

s2.595 <- c("Under the same conditions (common chemical classes), MCnebula outperforms the",
  "benchmark method (Fig. \\@ref(fig:summary)).",
  "The formula for Relative false rate:",
  "- RelativeFalseRate = 1 - (1 - FalseRate) * (1 - AverageLostRate)"
)

notShow1 <- new_section2(
  c("Combine figure..."),
  rblock({
    vis_summary <- into(grecta("a"), vis4)
    vis_classify <- visualize_comparison(dominant_res, res_molnet,
      from = c("MCnebula", "GNPS"))
    vis_classify <- into(grecta("b"), vis_classify)
    vis_id <- visualize_idRes(list(`No cut-off` = idRes, `0.5 cut-off` = idRes.5))
    vis_id <- into(grecta("c"), vis_id)
    frame1 <- frame_col(
      c(vis_classify = 1.5, vis_id = 1),
      namel(vis_classify, vis_id)
    )
    frame2 <- frame_row(
      c(vis_summary = 1, frame1 = 2),
      namel(vis_summary, frame1)
    )
    frame2 <- ggather(frame2, vp = viewport(, , .95, .95))
    pdf(paste0(tmp, "/compare_accuracy.pdf"), 14, 14)
    draw(frame2)
    dev.off()
  })
)

s2.7 <- new_heading("Evaluate accuracy of identification", 2)

s2.8 <- new_section2(
  c("Evaluate the accuracy of identification with origin dataset:"),
  rblock({
    identified <- dplyr::select(
      lst[[1]]$features_annotation, .features_id, inchikey2d, tani.score
    )
    identified <- dplyr::filter(identified, !is.na(tani.score))
    # nrow(identified) # [1] 6610
    index <- dplyr::filter(lst[[1]]$nebula_index, .features_id %in% identified$.features_id)
    index <- split(index, ~ class.name)
    index <- index[names(index) %in% dominant_res[[1]]$class.name]
    ref_inchikey2d <- dplyr::mutate(
      mgf_metadata, inchikey2d = stringr::str_extract(INCHIKEY, "^[A-Z]*")
    )
    ref_inchikey2d <- dplyr::select(ref_inchikey2d, .features_id = .id, inchikey2d)
    idRes <- stat_identification(index, identified, ref_inchikey2d)
    idRes.summary <- summarise(idRes)
  })
)

s2.81 <- new_section2(
  c("Set a cut-off for 'tani.score' (Tanimoto similarity)."),
  rblock({
    identified.5 <- dplyr::filter(identified, tani.score >= .5)
    index.5 <- lapply(index, dplyr::filter, .features_id %in% identified.5$.features_id)
    idRes.5 <- stat_identification(index.5, identified.5, ref_inchikey2d)
    idRes.5.summary <- summarise(idRes.5)
  })
)

s2.82 <- c(
  "For all identified compounds, the average false rate was",
  paste0(idRes.summary$false, "%."), "Set 0.5 for 'tani.score' as threshold,
  the average false rate was", paste0(idRes.5.summary$false, "%.")
)

s2.83 <- new_section2(
  c("Visualize the results:"),
  rblock({
    vis <- visualize_idRes(list(`No cut-off` = idRes, `0.5 cut-off` = idRes.5))
    pdf(f2.83 <- paste0(tmp, "/identified_accuracy.pdf"), 6, 9)
    draw(vis)
    dev.off()
  })
)

s2.84 <- include_figure(f2.83, "idres",
  "Identified accuracy of compounds in each classified chemical class")

s2.85 <- c("See results", paste0(ref(s2.84), "."))

s100 <- new_heading("Session infomation", 1)

s100.1 <- rblock({
  sessionInfo()
}, args = list(eval = T))

# ==========================================================================
# gather
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sections <- gather_sections()
report <- do.call(new_report, sections)
yaml(report)[1] <- c("title: Evaluation of MCnebula2")
render_report(report, file <- paste0(tmp, "/report.Rmd"))
rmarkdown::render(file)
