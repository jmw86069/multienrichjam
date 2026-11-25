
#' Define cache key for enrichment results
#' 
#' Define cache key for enrichment results
#' 
#' ## Todo
#' 
#' * Consider sorting `genedata`, `pathwaydata`, `universe`
#' to make each robust to different sort order at input.
#' 
#' @param genedata atomic vector of `character` genes or
#'    `numeric` gene scores.
#' @param pathwaydata `list` named by pathway, containing
#'    a `character` vector of genes.
#' @param enrich_method `character` string
#' @param universe `character` vector of genes
#' @param minGSSize,maxGSSize `numeric` minimum and maximum gene set
#'    size.
#' @param params `list` with any additional parameters as needed.
#' 
#' @returns `character` string output from `digest::digest()`
#' 
#' @keywords internal
#' @noRd
get_enrich_cache_key <- function(
   genedata=NULL,
   pathwaydata=NULL,
   enrich_method=NULL,
   universe=NULL,
   minGSSize=NULL,
   maxGSSize=NULL,
   params=NULL
) {
   #
   if (length(genedata) > 1 && is.atomic(genedata)) {
      genedata <- sort(genedata)
   }
   if (length(pathwaydata) >= 1) {
      if (inherits(pathwaydata, "list")) {
         pathwaydata <- jamba::uniques(jamba::mixedSorts(pathwaydata))
      } else if (is.atomic(pathwaydata)) {
         pathwaydata <- sort(pathwaydata)
      }
   }
   if (length(universe) > 1) {
      universe <- sort(universe)
   }
   cache_digest <- digest::digest(
      list(
         genedata=genedata,
         pathwaydata=pathwaydata,
         enrich_method=enrich_method,
         universe=universe,
         minGSSize=minGSSize,
         maxGSSize=maxGSSize,
         params=params
      )
   )
   return(cache_digest)
}

#' Create an enrichment cache file
#' 
#' Create an enrichment cache file, given a cache path and key.
#' 
#' @param cache_path `character` file path
#' @param cache_key `character` cache key from `get_enrich_cache_key()`
#' @param ... additional arguments are ignored.
#' 
#' @returns `character` full path to cache file.
#' 
#' @keywords internal
#' @noRd
get_enrich_cache_file <- function(
   cache_path,
   cache_key,
   ...
) {
   #
   cache_file <- file.path(
      cache_path,
      paste0(cache_key, ".rds")
   )
   return(cache_file)
}
