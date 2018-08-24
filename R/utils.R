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

#' Add cells while bypassing checks
#'
#' The \code{add.cells} method of \code{loom} class does a lot of checks on the dimensions of
#' data for cells to be added to a \code{loom} object. However, here, while I do need to add cells,
#' I already know the dimensions and that they work, so I can bypass the checks and save some time,
#' which can be a lot if there are many iterations.
#' Here I assume that the matrices will not be transposed. This function is basically a lightweight
#' version of the method \code{loom}.
#'
#' @param obj The \code{loom} object to which cells are to be added.
#' @param n_cells Number of cells to be added; the dimensions of the data are assumed to be
#' compatible with this number.
#' @param matrix.data A Matrix with \code{n_cells} rows and as many columns as there're genes;
#' data for new cells to be added to the \code{matrix} dataset of the \code{loom} object.
#' @param layers.data A named list of matrices with the same dimension as \code{matrix.data}
#' for layer data for the new cells.
#'
add_cells_lt <- function(obj, n_cells, matrix.data, layers.data) {
  dims.fill <- obj[['matrix']]$dims[1]
  dims.fill <- (dims.fill + 1L):(dims.fill + n_cells)
  # Add matrix data
  obj[['matrix']][dims.fill, ] <- matrix.data
  # Add layers
  layers.names <- names(x = obj[['layers']])
  for (i in layers.names) {
    obj[['layers']][[i]][dims.fill, ] <- layers.data[[i]]
  }
}

#' @importFrom progress progress_bar
NULL

#' Add column attributes for all the added cells
#'
#' In case there are many column attributes, and many chunks of new cells to be added, iterating
#' over all those attributes for each chunk would be very inefficient. This function adds the
#' column attributes to the new cells after they have all been added to the \code{loom} object.
#'
#' @param obj The \code{loom} object to which cells are to be added.
#' @param n_cells Number of cells to be added; the dimensions of the data are assumed to be
#' compatible with this number.
#' @param attributes.data A named list of matrices (with \code{n_cells} columns) or vectors
#' (with length \code{n_cells}) for column attributes for the new cells.
#'
add_col_attrs_lt <- function(obj, n_cells, attributes.data) {
  n_tot <- obj[['matrix']]$dims[1]
  dims.fill <- (n_tot - n_cells + 1):n_tot
  # Add column attributes
  attrs.names <- names(x = obj[['col_attrs']])
  matrices <- vapply(attributes.data, is.matrix, FUN.VALUE = logical(1L))
  pb <- progress_bar$new(format = "|:bar| :percent :eta :elapsed")
  counter <- 0; tot <- length(attributes.data)
  for (i in attrs.names) {
    if (matrices[i]) {
      if (is.infinite(obj[["col_attrs"]][[i]]$maxdims[2])) {
        obj[["col_attrs"]][[i]][,dims.fill] <- attributes.data[[i]]
      } else {
        this_attr <- obj[["col_attrs"]][[i]][,]
        this_attr <- list(cbind(this_attr, attributes.data[[i]]))
        names(this_attr) <- i
        obj$add.col.attribute(this_attr, overwrite = TRUE)
      }
    } else {
      if (is.infinite(obj[["col_attrs"]][[i]]$maxdims)) {
        obj[['col_attrs']][[i]][dims.fill] <- attributes.data[[i]]
      } else {
        this_attr <- obj[["col_attrs"]][[i]][]
        this_attr <- list(c(this_attr, attributes.data[[i]]))
        names(this_attr) <- i
        obj$add.col.attribute(this_attr, overwrite = TRUE)
      }
    }
    counter <- counter + 1
    pb$update(counter / tot)
  }
  pb$terminate()
  # Update shape
  obj$update.shape()
}
