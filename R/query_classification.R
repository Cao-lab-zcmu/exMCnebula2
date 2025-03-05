# ==========================================================================
# query classification for compounds using classyfire API
# run after query_inchikey
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases query_classification
#'
#' @title Query classification of compounds via packgage of 'classifyerR'
#'
#' @description The used function is:
#' \code{classyfireR::get_classification(inchikey)}
#' @family queries
#'
#' @name query_classification
NULL
#> NULL

#' @export query_classification
#' @aliases query_classification
#' @description \code{query_classification}: ...
#' @rdname query_classification
query_classification <- 
  function(
    inchikey2d,
    dir,
    inchikey.rdata = paste0(dir, "/inchikey.rdata"),
    rdata.name = "classification.rdata",
    classyfire_cl = NULL,
    planB = FALSE, 
    gather_as_rdata = T,
    ...
    ){
    rdata <- paste0(dir, "/", rdata.name)
    classes <- extract_rdata_list(rdata)
    if (!is.null(classes))
      inchikey2d <- inchikey2d[!inchikey2d %in% names(classes)]
    if(length(inchikey2d) == 0)
      return(paste0(dir, "/", rdata.name))
    inchikey_set <- extract_rdata_list(inchikey.rdata, inchikey2d)
    if (is.null(inchikey_set))
      stop("is.null(inchikey_set) == T. File `inchikey.rdata` may not exists.")
    sets <- lapply(inchikey_set, function(df){
      if("InChIKey" %in% colnames(df))
        return(df)
    })
    sets <- data.table::rbindlist(sets)
    sets <- dplyr::mutate(sets, inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]{1,}"))
    l <- classyfire_get_classification(sets, dir, classyfire_cl = classyfire_cl, ...)
    if (is.logical(l))
      return(paste0(dir, "/", rdata.name))
    if (gather_as_rdata) {
      cat("## gather data\n")
      packing_as_rdata_list(dir, pattern = "^[A-Z]{14}$",
        rdata = rdata.name, extra = classes)
    }
    return(paste0(dir, "/", rdata.name))
  }

#' @export classyfire_get_classification
#' @aliases classyfire_get_classification
#' @description \code{classyfire_get_classification}: ...
#' @rdname query_classification
classyfire_get_classification <-
  function(
    sets,
    dir,
    classyfire_cl = NULL,
    log_file = paste0(dir, "/classyfire.log"),
    planB = FALSE,
    ...
    ){
    if (file.exists(log_file)){
      log_df <- data.table::fread(log_file)
      sets <- dplyr::filter(sets, !InChIKey %in% log_df$log)
      if(nrow(sets) == 0)
        return(F)
    }
    sets <- split(data.frame(sets), ~ inchikey2d)
    log <- pbapply::pblapply(names(sets), cl = classyfire_cl,
      function(inchikey2d) {
        set <- sets[[ inchikey2d ]]
        unlist(lapply(set[["InChIKey"]], .get_classification,
            file = paste0(dir, "/", inchikey2d)),
          use.names = F)
      })
    log <- unlist(log, use.names = F)
    log <- data.frame(log = log)
    if (exists("log_df"))
      log <- dplyr::bind_rows(log_df, log)
    write_tsv(log, file = log_file)
  }

.get_classification <-
  function(inchikey, file, planB = FALSE) {
    if (!file.exists(file)) {
      if (planB == FALSE){ch <- classyfireR::get_classification(inchikey)} 
      else {ch <- get_classification_planB(inchikey)}
    } else{
      return()
    }
    if (is.null(ch)) {
      return(inchikey)
    } else{
      ch <- classyfireR::classification(ch)
      if (nrow(ch) != 0) write_tsv(ch, file) else return(inchikey)
    }
  }

#' @export
#' @import RSQLite

get_classification_planB <- function(inchi_key, conn=NULL)
{
  cache_hits <- 0
  if (! is.null(conn))  {
    qry <-
      RSQLite::dbSendQuery(conn,
                           "SELECT InChiKey,Classification FROM classyfire WHERE InChiKey=?")
    RSQLite::dbBind(qry, inchi_key)
    key <- RSQLite::dbFetch(qry)
    cache_hits <-RSQLite::dbGetRowCount(qry)
    RSQLite::dbClearResult(qry)
  }
  
  if (cache_hits==1) {
    object <- unserialize(charToRaw(key$Classification))
    message(crayon::green(clisymbols::symbol$tick, 'cached: ', inchi_key))
    return(object)
  } else {
    entity_url <- 'http://cfb.fiehnlab.ucdavis.edu/entities/'
    
    entity_query <- paste0(entity_url, inchi_key, '.json')
    
    response <- httr::RETRY(
      verb = "GET",
      url = entity_query,
      times = 10,
      terminate_on = c(404),
      quiet = T
    )
    
    if (response$status_code == 429) {
      stop('Request rate limit exceeded!')
    }
    
    if (response$status_code == 404) {
      message(crayon::red(clisymbols::symbol$cross, inchi_key))
    }
    
    if (response$status_code == 200) {
      text_content <- httr::content(response, 'text')
      
      if (text_content == '{}') {
        message(crayon::red(clisymbols::symbol$cross, inchi_key))
        return(invisible(NULL))
      } else{
        message(crayon::green(clisymbols::symbol$tick, inchi_key))
      }
      
      json_res <- jsonlite::fromJSON(text_content)
      
      classification <- classyfireR:::parse_json_output(json_res)
      
      
      object <- methods::new('ClassyFire')
      
      
      object@meta <-
        list(
          inchikey = json_res$inchikey,
          smiles = json_res$smiles,
          version = json_res$classification_version
        )
      
      object@classification <- classification
      
      if (length(json_res$direct_parent) > 0) {
        object@direct_parent <- json_res$direct_parent
      }
      
      if (length(json_res$alternative_parents) > 0) {
        object@alternative_parents <-
          tibble::tibble(
            name = json_res$alternative_parents$name,
            description = json_res$alternative_parents$description,
            chemont_id = json_res$alternative_parents$chemont_id,
            url = json_res$alternative_parents$url
          )
      } else{
        object@alternative_parents <- tibble::tibble()
      }
      
      if (length(json_res$predicted_chebi_terms) > 0) {
        object@predicted_chebi <- json_res$predicted_chebi_terms
      } else{
        object@predicted_chebi <- vector(mode = 'character')
      }
      
      
      if (length(json_res$external_descriptors) > 0) {
        object@external_descriptors <-
          parse_external_desc(json_res)
      } else{
        object@external_descriptors <- tibble::tibble()
      }
      
      if (length(json_res$description) > 0) {
        object@description <- json_res$description
      }
      
      
      if (! is.null(conn))  {
        qry2 <-
          RSQLite::dbSendQuery(
            conn,
            "INSERT INTO classyfire (InChiKey,InChi,Classification) VALUES(?,?,?)"
          )
        RSQLite::dbBind(qry2, list(inchi_key, "NULL", rawToChar(
          serialize(object, connection = NULL, ascii = TRUE)
        )))
        RSQLite::dbClearResult(qry2)
      }
      return(object)
    }
  }
}
