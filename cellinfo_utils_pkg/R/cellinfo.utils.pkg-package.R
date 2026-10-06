#' Utilities for cellinfo objects
#'
#' @description
#' **cellinfo.utils.pkg** provides helpers to create, refresh, annotate, subset,
#' seriate, and visualise `cellinfo` list objects used with Seurat single-cell
#' workflows. The functions were extracted from project `mlutils` setup scripts
#' and packaged with standard roxygen2 documentation.
#'
#' A `cellinfo` object is typically a named list containing at least:
#' \itemize{
#'   \item `cell.list`: named list of cell barcodes by group
#'   \item `markerlist`: named list of marker genes by group
#'   \item `cell.annotation` / `marker.annotation`: data frames of annotations
#'   \item derived fields such as `cells`, `markers`, `gaps.cells`, `gaps.markers`
#' }
#'
#' Some workflows still call companion helpers from the parent mlutils codebase
#' (for example `join_meta_exp2()`, `AddModuleScore3()`, `seriatecells()`,
#' `tpdf()`). Those symbols are not re-exported here.
#'
#' @keywords internal
"_PACKAGE"
