#' Get all column/row attributes
#'
#' The method get.attributes.df only returns attributes with one dimension. This function
#' will return a list with attributes with both 1 and 2 dimensions.
#'
#' @param obj The loom object whose attributes will be retrieved. Only 1 or 2 dimensions are
#' supported.
#' @param MARGIN Margin; 1 means row (genes), and 2 means column (cells)
#'
#' @return A list with all attributes of the margin

get_all_attributes <- function(obj, MARGIN) {
  if (!(MARGIN %in% c(1,2))) {stop("MARGIN must be either 1 or 2")}
  if (MARGIN == 1) {
    attr_names <- names(obj$row.attrs)
    attr_paths <- paste0("row_attrs/", attr_names)
  } else {
    attr_names <- names(obj$col.attrs)
    attr_paths <- paste0("col_attrs/", attr_names)
  }
  ndims <- sapply(attr_paths,
                  function(x) {
                    length(obj[[x]]$dims)
                  })
  res <- lapply(seq_along(attr_paths),
                function(i) {
                  if (ndims[i] == 2) {
                    obj[[attr_paths[i]]][,]
                  } else {
                    obj[[attr_paths[i]]][]
                  }
                })
  names(res) <- attr_names
  res
}

#' Extend column attributes when adding cells to loom object
#'
#' The purpose of this function is to make sure that the doublet loom object
#' can be properly combined with the copy of the real cell loom object.
#' In oroder to combine properly, the two loom objects must have the same column and row
#' attributes and the same layers. This function takes care of the column attributes by taking
#' \code{n_cells} rows from column attributes of the loom object with real cells, so the data
#' type will match. I can simply copy the row attribute (for genes) from the object with real cells
#' to that with doublets.
#'
#' @param n_cells Number of cells to add.
#' @param obj The loom object to which cells are to be added. Only column attributes with 1
#' or 2 dimensions are supported.
#' @param cell.names Name of dataset in loom object that stores cell names
#'
#' @return A data frame with dummy attributes of the right column names, right
#' data types, and the right dimensions.

extend_col_attrs <- function(n_cells, obj, cell.names = "cell_names") {
  # Get all the column attributes from obj, including those with 2 dimensions
  col_attrs_all <- get_all_attributes(obj, 2)
  # Randomly get n_cells
  n_total <- obj$matrix$dims[1]
  inds <- sample(1:n_total, n_cells)
  res <- lapply(col_attrs_all,
         function(x) {
           if (length(dim(x)) == 2) {
             x[,inds]
           } else {
             x[inds]
           }
         })
  res[[cell.names]] <- paste0("X", 1:n_cells)
  res
}

#' Create loom object for the simulated doublets
#'
#' This function creates a new loom object for the simulated doublets. Here I assume that you have
#' already run \code{NormalizeData}, \code{ScaleData}, and \code{FindVariableGenes} for the
#' loom object with real cells.
#'
#' @param n_cells Number of doublets.
#' @param mat The matrix with artificial doublet gene expression values
#' @param obj The loom object with real cells, used to get dummy column attributes to make
#' sure that the dummy has the right data type and dimensions.
#' @param overwrite Whether to overwrite if there's an existing file \code{doublets.loom}; this
#' is useful in case a previous run of \code{doubletFinder.loom} threw an error and failed to
#' remove the temporary loom files before exiting.
#'
#' @return A connection to the doublet loom file
#' @seealso \code{\link{extend_col_attrs}}
#'
create_doublet_loom <- function(n_cells, mat, obj, overwrite = FALSE) {
  n_genes <- obj[["matrix"]]$dims[2]
  row_attrs <- get_all_attributes(obj, 1)
  col_attrs <- extend_col_attrs(n_cells, obj)
  create("doublets.loom", do.transpose = FALSE,
         data = mat, gene.attrs = row_attrs, cell.attrs = col_attrs,
         layers = list(norm_data = mat, scale_data = mat),
         overwrite = overwrite)
}
