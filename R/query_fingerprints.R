# ==========================================================================
# Convert PubChem CID to tanimoto based chemical similarity
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases query_fingerprints
#' 
#' @title Convert PubChem CID to tanimoto based chemical similarity
#' 
#' @description bulk search for tanimoto via pubchem API
#' (http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/.../SDF)
#' @family queries
#' 
#' @name query_fingerprints
NULL
#> NULL

#' @export query_fingerprints
#' @aliases query_fingerprints
#' @description
#' \code{query_fingerprints}: ...
#' @rdname query_fingerprints
query_fingerprints <- function(
    cid,
    dir, 
    rdata.name = "fingerprints.rdata", 
    curl_cl = NULL,
    extra = NULL,
    gather_as_rdata = T, 
    group_number = 50,
    ...
) {
  rdata <- file.path(dir, rdata.name)
  cid_set <- extract_rdata_list(rdata)
  if (!is.null(cid_set)) {
    cid_set <- unname(cid_set)
    cid_set <- lapply(rapply(cid_set, enquote, how = "unlist"), eval)
    # cid_set <- data.table::rbindlist(cid_set)
    if (length(cid_set) > 0)
      extra <- cid_set
    else
      extra <- NULL
    cid <- cid[!cid %in% names(cid_set)]
    if(length(cid) == 0)
      return(paste0(dir, "/", rdata.name))
  } else {
    extra <- NULL
  }
  group <- grouping_vec2list(cid, group_number = group_number)
  pbapply::pblapply(
    group, pubchem_get_fingerprints,
    dir = dir, ..., cl = curl_cl
  )
  if (gather_as_rdata) {
    cat("## gather data\n")
    packing_as_rdata_list(
      dir, 
      pattern = "^G[0-9]+\\.qs$",
      dedup = F,
      rdata = rdata.name,
      extra = extra
    )
  }
  return(paste0(dir, "/", rdata.name))
}

#' @export pubchem_get_fingerprints
#' @aliases pubchem_get_fingerprints
#' @description \code{pubchem_get_fingerprints}: ...
#' @rdname query_fingerprints
pubchem_get_fingerprints <- function(
  cid, 
  dir,
  ...
) {
  name_cid <- cid
  savename <- attr(cid, "name")
  file <- paste0(dir, "/", savename)
  cid <- paste0(cid, collapse = ",")
  url_start <- "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  url_end <- "/SDF"
  url <- paste0(url_start, cid, url_end)
  
  text <- httr::GET(url)
  if (text$status != 404) {
    text <- rawToChar(httr::content(text))
    text <- unlist(strsplit(text, "\n"))
  } else return(list(NA))
  
  y <- regexpr("^\\${4,4}", text, perl = TRUE)
  index <- which(y != -1)
  indexDF <- data.frame(
    start = c(1, index[-length(index)] + 1), 
    end = index
  )
  sdf_list <- lapply(
    seq(along = indexDF[, 1]), 
    function(x) {
      text[seq(indexDF[x, 1], indexDF[x, 2])]
    }
  )
  names(sdf_list) <- name_cid
  
  qs::qsave(sdf_list, file = paste0(file, ".qs"))
  return(sdf_list)
}

#' @export packing_as_rdata_list
#' @aliases packing_as_rdata_list
#' @description \code{packing_as_rdata_list}: gather table as .rdata
#' @rdname query_synonyms
packing_as_rdata_list <- function(
    path,
    pattern,
    rdata,
    extra = NULL,
    rm_files = T,
    dedup = T
) {
    file_set <- list.files(path, pattern = pattern)
    if (length(file_set) == 0)
      return()
    read_fun <- if (all(grepl("\\.qs$", file_set))) qs::qread else read_tsv
    list <- pbapply::pblapply(file.path(path, file_set), read_fun)
    names(list) <- file_set
    list <- c(extra, list)
    if(dedup){
      df <- data.table::data.table(name = names(list), n = 1:length(list))
      df <- dplyr::distinct(df, name, .keep_all = T)
      list <- list[df$n]
    }
    if(rm_files){
      lapply(file.path(path, file_set), file.remove)
    }
    save(list, file = file.path(path, rdata))
}

#' @export sdf_convert_tanimoto
#' @import chemmineR
#' @aliases sdf_convert_tanimoto
#' @description \code{sdf_convert_tanimoto}: ...
#' @rdname query_fingerprints
sdf_convert_tanimoto <- function(
    sdf_list = NULL,
    dir, 
    rdata.name = "fingerprints.rdata"
) {
  if (is.null(sdf_list)) {
    rdata <- file.path(dir, rdata.name)
    if (file.exists(rdata)) {
      sdf_list <- exMCnebula2::extract_rdata_list(rdata)
      sdf_list <- unname(sdf_list)
      sdf_list <- lapply(rapply(sdf_list, enquote, how = "unlist"), eval)
    } else stop("Error:Missing sdf_list or rdata file")
  } 
  cmpd.sdf.list <- new("SDFstr", a = sdf_list)
  sd.list <- as(cmpd.sdf.list, "SDFset")
  cid(sd.list) <- sdfid(sd.list)
  fpset <- fp2bit(sd.list, type = 2)
  out <- sapply(rownames(fpset), function(x) {
    ChemmineR::fpSim(x = fpset[x,], fpset,sorted = FALSE) }) 
  mat <- as.matrix(out)
  id <- is.na(mat) # used to allow missing
  mat[id] <- "nna"
  mat[lower.tri(mat)] <- "na" # use to allow missing values
  diag(mat) <- "na"
  obj <- reshape2::melt(mat)
  colnames(obj) <- c("source", "target", "value")
  obj <- obj[!obj$value == "na", ]
  obj$value[obj$value == "nna"] <- NA
  obj$value <- as.numeric(as.character(obj$value)) 
  return(obj)
}
