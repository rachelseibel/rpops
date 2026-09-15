#' Delineate management units
#'
#' Groups raster cells into farm-level management units. A single grower
#' decision applies per unit. Three methods are supported: contiguous host
#' patches, polygon field boundaries (e.g. CROME), and a regular grid.
#'
#' @param host_raster SpatRaster of host density. Cells with value 0 or NA
#'   are treated as non-host and assigned NA in the output.
#' @param method One of \code{"patch"} (default), \code{"polygon"}, or
#'   \code{"grid"}.
#' @param polygon_layer SpatVector of field boundary polygons. Required when
#'   \code{method = "polygon"}.
#' @param grid_size Numeric. Side length (metres) of each grid cell. Required
#'   when \code{method = "grid"}.
#'
#' @return Integer SpatRaster where each unique positive value is a management
#'   unit ID. Non-host cells are NA.
#'
#' @importFrom terra patches rasterize is.na values<- crs res
#' @export
delineate_management_units <- function(host_raster,
                                       method        = c("patch", "polygon", "grid"),
                                       polygon_layer = NULL,
                                       grid_size     = NULL) {
  method <- match.arg(method)

  # Mask: only host cells (value > 0 and not NA) are delineated
  host_mask <- host_raster
  terra::values(host_mask)[terra::values(host_mask) == 0] <- NA

  unit_rast <- switch(
    method,
    patch   = delineate_patches(host_mask),
    polygon = delineate_from_polygons(host_mask, polygon_layer),
    grid    = delineate_grid(host_mask, grid_size)
  )

  unit_rast
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @importFrom terra patches
delineate_patches <- function(host_mask) {
  terra::patches(host_mask, directions = 8, zeroAsNA = TRUE)
}

#' @importFrom terra rasterize patches
delineate_from_polygons <- function(host_mask, polygon_layer) {
  if (is.null(polygon_layer))
    stop("polygon_layer is required when method = 'polygon'")

  # Rasterise polygon IDs onto the host raster extent/resolution
  poly_id_rast <- terra::rasterize(polygon_layer, host_mask,
                                   field = 1:nrow(polygon_layer))

  # Mask to host cells only
  unit_rast <- terra::mask(poly_id_rast, host_mask)

  # Any host cells not covered by a polygon fall back to patch delineation
  uncovered <- is.na(terra::values(unit_rast)) & !is.na(terra::values(host_mask))
  if (any(uncovered)) {
    fallback       <- delineate_patches(host_mask)
    max_poly_id    <- max(terra::values(unit_rast), na.rm = TRUE)
    fb_vals        <- terra::values(fallback)
    fb_vals[!is.na(fb_vals)] <- fb_vals[!is.na(fb_vals)] + max_poly_id
    unit_rast_vals <- terra::values(unit_rast)
    unit_rast_vals[uncovered] <- fb_vals[uncovered]
    terra::values(unit_rast) <- unit_rast_vals
  }

  unit_rast
}

#' @importFrom terra rast ext crs res values<-
delineate_grid <- function(host_mask, grid_size) {
  if (is.null(grid_size) || grid_size <= 0)
    stop("grid_size must be a positive number when method = 'grid'")

  cell_res <- terra::res(host_mask)[1]
  cells_per_side <- max(1L, round(grid_size / cell_res))

  nr <- terra::nrow(host_mask)
  nc <- terra::ncol(host_mask)

  grid_row <- ceiling(seq_len(nr) / cells_per_side)
  grid_col <- ceiling(seq_len(nc) / cells_per_side)
  n_grid_cols <- max(grid_col)

  # Build cell-level grid ID matrix (row-major)
  id_mat <- outer(grid_row, grid_col, function(r, c) (r - 1L) * n_grid_cols + c)

  unit_rast <- host_mask
  # id_mat is [row, col]; terra values are row-major, so flatten transposed
  # (as.integer(id_mat) would be column-major and transpose the grid -> striping).
  terra::values(unit_rast) <- as.integer(t(id_mat))
  terra::values(unit_rast)[is.na(terra::values(host_mask))] <- NA

  unit_rast
}
