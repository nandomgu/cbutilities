#' @keywords internal
#' @noRd
fcat <- function(...) {
  cat(paste(..., sep = " "), "\n")
}

#' Set names on an object
#'
#' @param x Object to name.
#' @param nam Character vector of names.
#'
#' @return `x` with `names(x)` set to `nam`.
#' @keywords internal
#' @noRd
givename <- function(x, nam) {
  names(x) <- nam
  x
}

#' Set names on an object (alias)
#'
#' @inheritParams givename
#' @keywords internal
#' @noRd
givenames <- function(x, nam) {
  givename(x, nam)
}

#' Set column names (optionally at selected indices)
#'
#' @param x A matrix-like object.
#' @param ind Optional integer indices of columns to rename. If `NULL`, all
#'   columns are renamed.
#' @param nms Character vector of new names.
#'
#' @return `x` with updated column names.
#' @keywords internal
#' @noRd
givecolnames <- function(x, ind = NULL, nms) {
  if (is.null(ind)) {
    colnames(x) <- nms
  } else {
    colnames(x)[ind] <- nms
  }
  x
}

#' Set row names
#'
#' @param x A matrix-like object.
#' @param nms Character vector of row names.
#'
#' @return `x` with updated row names.
#' @keywords internal
#' @noRd
giverownames <- function(x, nms) {
  rownames(x) <- nms
  x
}

#' Remove `NA` values from a vector
#'
#' @param x A vector.
#'
#' @return `x` without `NA` entries.
#' @keywords internal
#' @noRd
removenas <- function(x) {
  x[!is.na(x)]
}

#' Move rownames into a column
#'
#' @param df A data frame.
#' @param colname Name of the new column that will store rownames.
#'
#' @return A data frame with rownames stored in `colname` and rownames cleared.
#' @keywords internal
#' @noRd
names2col <- function(df, colname) {
  df[[colname]] <- rownames(df)
  rownames(df) <- NULL
  df
}

#' Move a column into rownames
#'
#' @param df A data frame.
#' @param colname Column to use as rownames (removed from columns).
#'
#' @return A data frame with rownames taken from `colname`.
#' @keywords internal
#' @noRd
col2names <- function(df, colname) {
  rn <- df[[colname]]
  df[[colname]] <- NULL
  rownames(df) <- rn
  df
}

#' Subset rows by name, preserving missing rows as `NA`
#'
#' @param df A data frame.
#' @param rows Character vector of row names to select.
#'
#' @return A data frame with rows ordered as in `rows`.
#' @keywords internal
#' @noRd
getrows <- function(df, rows) {
  missing <- setdiff(rows, rownames(df))
  if (length(missing) > 0) {
    filler <- df[rep(NA_integer_, length(missing)), , drop = FALSE]
    rownames(filler) <- missing
    df <- rbind(df, filler)
  }
  df[rows, , drop = FALSE]
}

#' Jaccard / set overlap of two character vectors
#'
#' @param vector1 First vector.
#' @param vector2 Second vector.
#'
#' @return Numeric overlap in \\eqn{[0, 1]}.
#' @keywords internal
#' @noRd
zsoverlap <- function(vector1, vector2) {
  vector1 <- unique(as.character(vector1))
  vector2 <- unique(as.character(vector2))
  if (length(vector1) == 0 && length(vector2) == 0) {
    return(1)
  }
  length(intersect(vector1, vector2)) / length(union(vector1, vector2))
}

#' Identify rows with non-zero variance
#'
#' @param mat A numeric matrix.
#'
#' @return Logical vector, `TRUE` for variant rows.
#' @keywords internal
#' @noRd
get.variant.rows <- function(mat) {
  apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
}

#' Random distinct-ish colors
#'
#' @param n Number of colors.
#'
#' @return Character vector of hex colors.
#' @keywords internal
#' @noRd
randomcolors <- function(n) {
  grDevices::rainbow(n)
}

#' Trim each character vector in a named marker list to length `n`
#'
#' @param markerlist Named list of marker gene vectors.
#' @param n Maximum markers to keep per group.
#'
#' @return Trimmed marker list.
#' @keywords internal
#' @noRd
adjust.markerlist <- function(markerlist, n) {
  lapply(markerlist, function(x) {
    if (length(x) > n) {
      x[seq_len(n)]
    } else {
      x
    }
  })
}
