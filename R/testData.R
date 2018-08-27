#' Test counts data.
#'
#' @title Counts data used to demonstrate sp.scRNAseq package.
#' @docType data
#' @name testCounts
#' @format matrix
#' \describe{
#'     \item{rownames}{Gene names}
#'     \item{colnames}{Samples. Singlets prefixed with "s" and multiplets "m".}
#' }
#' @usage testCounts
#' @return Matrix of counts.
#' @examples
#' data(testCounts)
#'
NULL

#' Test counts meta data.
#'
#' @title Meta data for the testCounts dataset.
#' @docType data
#' @name testMeta
#' @format tibble
#' \describe{
#'     \item{sample}{Sample ID}
#'     \item{cellNumber}{Indicates if samples are singlets or multiplets.}
#'     \item{cellTypes}{Indicates the cell types included in the sample.}
#' }
#' @usage testMeta
#' @return Tibble with meta data.
#' @examples
#' data(testMeta)
#'
NULL
