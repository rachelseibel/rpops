#' Assign grower types to management units
#'
#' Assigns each management unit to one of three grower types — early adopter
#' (1), responsive (2), or non-adopter (3) — using one of three spatial
#' structures: random, spatially correlated (clustered via a Gaussian random
#' field), or empirically constrained from a survey-derived adoption
#' probability raster.
#'
#' @param unit_raster Integer SpatRaster of management unit IDs (output of
#'   \code{\link{delineate_management_units}}).
#' @param type_probs Named numeric vector summing to 1 with elements
#'   \code{early_adopter}, \code{responsive}, and \code{non_adopter}.
#' @param spatial_structure One of \code{"random"} (default),
#'   \code{"clustered"}, or \code{"empirical"}.
#' @param cluster_range Numeric. Spatial range parameter (metres) for the
#'   Gaussian random field. Required when \code{spatial_structure = "clustered"}.
#' @param empirical_raster SpatRaster of adoption probability values [0, 1],
#'   co-registered with \code{unit_raster}. Required when
#'   \code{spatial_structure = "empirical"}.
#' @param random_seed Integer seed for reproducibility.
#'
#' @return Integer SpatRaster with the same extent/resolution as
#'   \code{unit_raster}: 1 = early adopter, 2 = responsive, 3 = non-adopter,
#'   NA = non-host buffer.
#'
#' @importFrom terra values rast ext res crs nrow ncol
#' @export
assign_grower_types <- function(unit_raster,
                                type_probs,
                                spatial_structure = c("random", "clustered", "empirical", "graded"),
                                cluster_range     = NULL,
                                cluster_strength  = 1,
                                empirical_raster  = NULL,
                                random_seed       = NULL,
                                unit_aggregation  = c("field_rank", "modal")) {
  spatial_structure <- match.arg(spatial_structure)
  unit_aggregation  <- match.arg(unit_aggregation)

  if (!is.null(random_seed)) set.seed(random_seed)

  # Validate type_probs
  required_names <- c("early_adopter", "responsive", "non_adopter")
  if (!all(required_names %in% names(type_probs)))
    stop("type_probs must have names: early_adopter, responsive, non_adopter")
  if (abs(sum(type_probs) - 1) > 1e-6)
    stop("type_probs must sum to 1")

  unit_vals <- terra::values(unit_raster, mat = FALSE)
  unit_ids  <- sort(unique(unit_vals[!is.na(unit_vals) & unit_vals > 0]))

  type_assignment <- switch(
    spatial_structure,
    random   = assign_types_random(unit_ids, type_probs),
    clustered = assign_types_clustered(unit_raster, unit_ids, type_probs,
                                       cluster_range, unit_aggregation),
    graded    = assign_types_graded(unit_raster, unit_ids, type_probs,
                                    cluster_range, cluster_strength),
    empirical = assign_types_empirical(unit_raster, unit_ids, type_probs,
                                       empirical_raster)
  )

  # Build output raster: same template, fill from assignment vector
  type_rast <- unit_raster
  out_vals  <- rep(NA_integer_, terra::ncell(unit_raster))
  cell_ids  <- which(!is.na(unit_vals) & unit_vals > 0)

  for (idx in cell_ids) {
    uid <- unit_vals[idx]
    out_vals[idx] <- type_assignment[[as.character(uid)]]
  }

  terra::values(type_rast) <- out_vals
  type_rast
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Random (spatially uncorrelated) type assignment
#'
#' Each unit draws its type independently from a multinomial.
assign_types_random <- function(unit_ids, type_probs) {
  n      <- length(unit_ids)
  types  <- sample(c(1L, 2L, 3L), size = n, replace = TRUE,
                   prob = type_probs[c("early_adopter", "responsive", "non_adopter")])
  stats::setNames(as.list(types), as.character(unit_ids))
}

#' Clustered type assignment via a Gaussian random field
#'
#' Simulates a zero-mean Gaussian random field with exponential covariance and
#' range parameter \code{cluster_range} (metres), thresholds at the quantiles
#' implied by \code{type_probs}, then assigns each unit the majority type among
#' its constituent cells.
#'
#' Requires the \code{fields} package.
#'
#' @importFrom terra values res nrow ncol ext
assign_types_clustered <- function(unit_raster, unit_ids, type_probs, cluster_range,
                                   unit_aggregation = c("field_rank", "modal")) {
  unit_aggregation <- match.arg(unit_aggregation)
  if (is.null(cluster_range) || cluster_range <= 0)
    stop("cluster_range must be a positive number when spatial_structure = 'clustered'")
  if (!requireNamespace("fields", quietly = TRUE))
    stop("Package 'fields' is required for clustered type assignment. ",
         "Install with: install.packages('fields')")

  nr  <- terra::nrow(unit_raster)
  nc  <- terra::ncol(unit_raster)

  # Build grid in CELL-INDEX coordinates (spacing = 1 cell) so that
  # cluster_range is expressed in CELLS: cluster_range = 8 gives an exponential
  # correlation range of 8 cells (= 8 km at 1 km resolution). NOTE: an earlier
  # version laid the grid out in metres (spacing = res_m) while passing
  # cluster_range straight to aRange, so a range of "8" meant 8 metres against
  # 1000 m cells -- i.e. white noise and no real clustering. Working in cell
  # units removes that resolution dependence.
  row_centres <- seq_len(nr) - 0.5
  col_centres <- seq_len(nc) - 0.5

  # Simulate one realisation of a Gaussian random field
  obj   <- fields::circulantEmbeddingSetup(
    grid = list(x = col_centres, y = row_centres),
    cov.args = list(Covariance = "Exponential", aRange = cluster_range)
  )
  field_vals <- as.vector(fields::circulantEmbedding(obj))

  # ---------------------------------------------------------------------------
  # Reduce the smooth cell-level field to a UNIT-level assignment.
  #
  # Both spatial structures must realise the SAME marginal composition
  # (type_probs); only the spatial ARRANGEMENT differs (random = independent
  # per unit; clustered = contiguous patches). The default "field_rank" method
  # takes one field value per unit (its mean), then cuts the ranked units at
  # the exact counts implied by type_probs. This preserves the unit-level
  # marginal by construction while keeping spatial autocorrelation, because
  # neighbouring units have similar mean field values and so share a type.
  #
  # unit_aggregation = "modal" reproduces the earlier behaviour (each unit
  # takes the majority cell type). That does NOT preserve the marginal -- the
  # minority class is hollowed out and the plurality class inflated -- so it is
  # retained for backward compatibility only, not recommended for the
  # random-vs-clustered arrangement experiment.
  unit_vals <- terra::values(unit_raster, mat = FALSE)

  if (identical(unit_aggregation, "modal")) {
    q1 <- stats::quantile(field_vals, probs = type_probs["non_adopter"],
                          na.rm = TRUE)
    q2 <- stats::quantile(field_vals, probs = type_probs["non_adopter"] +
                            type_probs["responsive"], na.rm = TRUE)
    cell_types <- ifelse(field_vals <= q1, 3L,
                  ifelse(field_vals <= q2, 2L, 1L))
    assignment <- list()
    for (uid in unit_ids) {
      idx   <- which(unit_vals == uid)
      modal <- as.integer(names(sort(table(cell_types[idx]),
                                     decreasing = TRUE)[1]))
      assignment[[as.character(uid)]] <- modal
    }
    return(assignment)
  }

  # field_rank (default): one field value per unit, cut ranked units at the
  # exact counts implied by type_probs.
  unit_field <- vapply(unit_ids, function(uid)
                       mean(field_vals[which(unit_vals == uid)], na.rm = TRUE),
                       numeric(1))
  n      <- length(unit_ids)
  n_non  <- min(round(type_probs["non_adopter"]  * n), n)
  n_resp <- min(round(type_probs["responsive"]   * n), n - n_non)
  ord    <- order(unit_field)               # ascending: low field -> non-adopter
  utypes <- integer(n)
  if (n_non  > 0)          utypes[ord[seq_len(n_non)]]          <- 3L
  if (n_resp > 0)          utypes[ord[n_non + seq_len(n_resp)]] <- 2L
  if (n_non + n_resp < n)  utypes[ord[(n_non + n_resp + 1L):n]] <- 1L
  stats::setNames(as.list(utypes), as.character(unit_ids))
}

#' Cut a per-unit covariate field into types at the exact composition counts
#'
#' Shared field_rank logic: rank units by \code{unit_field} and assign the
#' lowest-ranked to non-adopter, next to responsive, highest to early-adopter,
#' cutting at the counts implied by \code{type_probs}. Preserves the marginal
#' composition exactly, whatever the spatial pattern of \code{unit_field}.
.cut_field_rank <- function(unit_field, unit_ids, type_probs) {
  n      <- length(unit_ids)
  n_non  <- min(round(type_probs["non_adopter"] * n), n)
  n_resp <- min(round(type_probs["responsive"]  * n), n - n_non)
  ord    <- order(unit_field)               # ascending: low field -> non-adopter
  utypes <- integer(n)
  if (n_non  > 0)          utypes[ord[seq_len(n_non)]]          <- 3L
  if (n_resp > 0)          utypes[ord[n_non + seq_len(n_resp)]] <- 2L
  if (n_non + n_resp < n)  utypes[ord[(n_non + n_resp + 1L):n]] <- 1L
  stats::setNames(as.list(utypes), as.character(unit_ids))
}

#' Graded type assignment: continuous control of clustering intensity
#'
#' Unifies the random and clustered arrangements on a single continuous axis so
#' that clustering can be swept as a dose-response rather than a two-level
#' contrast. Builds two per-unit covariate fields -- a spatially correlated one
#' (Gaussian random field, exponential range \code{cluster_range} in CELLS,
#' identical construction to \code{assign_types_clustered}) and an independent
#' random one -- then blends them by the strength parameter
#' \code{cluster_strength} = rho in [0, 1]:
#'
#'   z(rho) = rho * z_clust + sqrt(1 - rho^2) * z_rand
#'
#' The sqrt(1 - rho^2) weight holds the blended field's variance constant, so
#' rho is exactly the correlation between the blend and the pure clustered
#' field. Limiting cases:
#'   rho = 0 -> pure random arrangement (equivalent to spatial_structure="random")
#'   rho = 1 -> fully clustered at cluster_range (equivalent to "clustered")
#' Intermediate rho gives graded patch INTENSITY at fixed patch SIZE. Varying
#' cluster_range at rho = 1 instead gives graded patch SIZE (correlation-length
#' sweep). The marginal composition (type_probs) is preserved for every rho by
#' the shared field_rank cut, so only the ARRANGEMENT changes across the sweep.
#'
#' @importFrom terra values nrow ncol
assign_types_graded <- function(unit_raster, unit_ids, type_probs,
                                cluster_range, cluster_strength = 1) {
  rho <- cluster_strength
  if (is.null(rho) || rho < 0 || rho > 1)
    stop("cluster_strength (rho) must be in [0, 1] for spatial_structure = 'graded'")
  if (is.null(cluster_range) || cluster_range <= 0)
    stop("cluster_range must be a positive number when spatial_structure = 'graded'")
  if (!requireNamespace("fields", quietly = TRUE))
    stop("Package 'fields' is required for graded type assignment.")

  unit_vals <- terra::values(unit_raster, mat = FALSE)

  # Independent random covariate, one value per unit
  z_rand <- stats::rnorm(length(unit_ids))

  # Spatially correlated covariate: GRF over cell-index grid (as clustered),
  # then averaged to one value per unit.
  nr <- terra::nrow(unit_raster); nc <- terra::ncol(unit_raster)
  # circulantEmbedding can fail ("Weight function has negative values") when the
  # correlation range approaches the grid extent; enlarge the circulant grid M
  # and retry, escalating the padding factor until it embeds. M must be >= the
  # data grid; use the next power of two above (factor * dim).
  nextpow2 <- function(x) 2^ceiling(log2(x))
  field_vals <- NULL
  for (fac in c(2, 3, 4, 6, 8)) {
    Mx <- nextpow2(fac * nc); My <- nextpow2(fac * nr)
    obj <- tryCatch(
      fields::circulantEmbeddingSetup(
        grid = list(x = seq_len(nc) - 0.5, y = seq_len(nr) - 0.5),
        M = c(Mx, My),
        cov.args = list(Covariance = "Exponential", aRange = cluster_range)),
      error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(obj)) {
      fv <- tryCatch(as.vector(fields::circulantEmbedding(obj)),
                     error = function(e) NULL, warning = function(w) NULL)
      if (!is.null(fv) && all(is.finite(fv))) { field_vals <- fv; break }
    }
  }
  if (is.null(field_vals))
    stop("graded GRF embedding failed even with padding; ",
         "reduce cluster_range relative to the grid extent")
  z_clust <- vapply(unit_ids, function(uid)
                    mean(field_vals[which(unit_vals == uid)], na.rm = TRUE),
                    numeric(1))

  # Standardise both to unit variance, then blend at constant variance.
  zsc <- function(z) { s <- stats::sd(z); if (is.na(s) || s == 0) z - mean(z) else (z - mean(z)) / s }
  z_clust <- zsc(z_clust); z_rand <- zsc(z_rand)
  unit_field <- rho * z_clust + sqrt(1 - rho^2) * z_rand

  .cut_field_rank(unit_field, unit_ids, type_probs)
}

#' Empirical type assignment from a survey-derived adoption probability raster
#'
#' For each unit, extracts the mean adoption probability from
#' \code{empirical_raster}, draws an early-adopter Bernoulli trial, then
#' splits non-early-adopter units between responsive and non-adopter in
#' proportion to \code{type_probs}.
#'
#' @importFrom terra extract values
assign_types_empirical <- function(unit_raster, unit_ids, type_probs,
                                   empirical_raster) {
  if (is.null(empirical_raster))
    stop("empirical_raster is required when spatial_structure = 'empirical'")

  unit_vals   <- terra::values(unit_raster, mat = FALSE)
  emp_vals    <- terra::values(empirical_raster, mat = FALSE)

  # Proportion of responsive among non-early-adopters
  p_resp_given_not_ea <- type_probs["responsive"] /
                         (type_probs["responsive"] + type_probs["non_adopter"])

  assignment <- list()
  for (uid in unit_ids) {
    idx      <- which(unit_vals == uid)
    mean_p   <- mean(emp_vals[idx], na.rm = TRUE)
    mean_p   <- max(0, min(1, mean_p))  # clamp

    if (stats::runif(1) < mean_p) {
      assignment[[as.character(uid)]] <- 1L  # early adopter
    } else if (stats::runif(1) < p_resp_given_not_ea) {
      assignment[[as.character(uid)]] <- 2L  # responsive
    } else {
      assignment[[as.character(uid)]] <- 3L  # non-adopter
    }
  }
  assignment
}
