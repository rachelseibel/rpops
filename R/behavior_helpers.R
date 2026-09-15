#' Parse and validate a behavior module configuration
#'
#' Reads the \code{behavior:} block from a YAML config file (or list), injects
#' it into the base PoPS config, and calls \code{validate_behavior_config()}.
#'
#' @param behavior_config_file Path to a YAML file that contains a top-level
#'   \code{behavior:} key, or a list already parsed from such a file.
#' @param pops_config Optional base PoPS config list (output of
#'   \code{configuration()}). If supplied, the behavior parameters are merged
#'   into it and the merged list is returned. If NULL, only the behavior
#'   sub-list is returned.
#'
#' @return A list with class \code{"behavior_config"} (or merged with pops
#'   config). On validation failure, \code{$failure} is set to an error string.
#'
#' @importFrom yaml read_yaml
#' @export
behavior_configuration <- function(behavior_config_file, pops_config = NULL) {
  if (is.list(behavior_config_file)) {
    raw <- behavior_config_file
  } else {
    raw <- read_yaml(behavior_config_file)
  }

  if (!is.null(raw[["behavior"]])) {
    bcfg <- raw[["behavior"]]
  } else {
    bcfg <- raw
  }

  check <- validate_behavior_config(bcfg)
  if (!check$checks_passed) {
    if (!is.null(pops_config)) {
      pops_config$behavior_failure <- check$failed_check
      return(pops_config)
    }
    return(list(failure = check$failed_check))
  }

  if (!is.null(pops_config)) {
    pops_config$behavior <- bcfg
    return(pops_config)
  }
  structure(bcfg, class = c("behavior_config", "list"))
}

# ---------------------------------------------------------------------------

#' Validate a behavior configuration list
#'
#' Checks all required fields are present and within valid ranges.
#'
#' @param bcfg Named list of behavior parameters (the \code{behavior:} block).
#' @param start_date,end_date Optional ISO-8601 schedule bounds (the
#'   simulation's \code{start_date} / \code{end_date}). When supplied, each
#'   \code{decision_date} is checked to fall within
#'   \code{[start_date, end_date]}, and the last decision date plus
#'   \code{pesticide_duration} is checked to stay on or before
#'   \code{end_date}. This prevents the cryptic C++ "Date is outside of
#'   schedule" error.
#' @param time_step Optional time-step unit (\code{"day"}, \code{"week"}, or
#'   \code{"month"}); used to convert \code{pesticide_duration} into a date
#'   offset for the buffer check. Defaults to \code{"day"}.
#' @return List with \code{$checks_passed} (logical) and
#'   \code{$failed_check} (character or NULL).
#' @export
validate_behavior_config <- function(bcfg, start_date = NULL, end_date = NULL,
                                     time_step = "day") {
  fail <- function(msg) list(checks_passed = FALSE, failed_check = msg)

  # Required top-level fields
  required <- c("type_probs", "decision_dates", "grower_params",
                 "treatment_efficacy", "pesticide_duration")
  missing_fields <- required[!required %in% names(bcfg)]
  if (length(missing_fields))
    return(fail(paste("Missing behavior config fields:", paste(missing_fields, collapse = ", "))))

  # type_probs
  tp <- bcfg$type_probs
  if (!all(c("early_adopter", "responsive", "non_adopter") %in% names(tp)))
    return(fail("type_probs must contain: early_adopter, responsive, non_adopter"))
  if (abs(sum(unlist(tp)) - 1) > 1e-6)
    return(fail("type_probs must sum to 1"))

  # Probability fields
  prob_fields <- c("treatment_efficacy")
  for (f in prob_fields) {
    v <- bcfg[[f]]
    if (!is.numeric(v) || v < 0 || v > 1)
      return(fail(paste(f, "must be numeric in [0, 1]")))
  }

  # pesticide_duration
  if (!is.numeric(bcfg$pesticide_duration) || bcfg$pesticide_duration <= 0)
    return(fail("pesticide_duration must be a positive number"))

  # grower_params: check per-type probability fields
  for (tp_name in c("early_adopter", "responsive", "non_adopter")) {
    if (!tp_name %in% names(bcfg$grower_params))
      return(fail(paste("grower_params missing type:", tp_name)))
    gp <- bcfg$grower_params[[tp_name]]
    for (pf in c("willingness_to_treat", "detection_probability")) {
      v <- gp[[pf]]
      if (is.null(v) || !is.numeric(v) || v < 0 || v > 1)
        return(fail(paste("grower_params$", tp_name, "$", pf,
                          "must be numeric in [0, 1]", sep = "")))
    }
    if (is.null(gp$perception_radius) || gp$perception_radius <= 0)
      return(fail(paste("grower_params$", tp_name,
                        "$perception_radius must be a positive integer", sep = "")))
  }

  # Spatial-structure-specific checks
  ss <- bcfg$spatial_structure
  if (!is.null(ss)) {
    if (ss == "clustered" && (is.null(bcfg$cluster_range) || bcfg$cluster_range <= 0))
      return(fail("cluster_range must be a positive number when spatial_structure = 'clustered'"))
    if (ss == "empirical" && is.null(bcfg$empirical_raster_file))
      return(fail("empirical_raster_file must be specified when spatial_structure = 'empirical'"))
  }

  # decision_dates must be parseable, non-empty, and (when bounds are known)
  # fall within the simulation schedule with room for the treatment duration.
  dd <- tryCatch(as.Date(unlist(bcfg$decision_dates)),
                 error = function(e) NA)
  if (length(dd) == 0 || any(is.na(dd)))
    return(fail("decision_dates must be a non-empty vector of ISO-8601 dates"))

  if (!is.null(start_date) && !is.null(end_date)) {
    sd <- as.Date(start_date)
    ed <- as.Date(end_date)
    if (any(dd < sd) || any(dd > ed))
      return(fail(sprintf(
        "decision_dates must fall within the simulation schedule [%s, %s]; out-of-range: %s",
        sd, ed, paste(as.character(dd[dd < sd | dd > ed]), collapse = ", "))))

    # Buffer: last decision + pesticide_duration must stay within the schedule.
    # Use the longest duration across types (per-type overrides may exceed the
    # global value) so the check is conservative.
    step_days <- switch(time_step %||% "day",
                        day = 1, week = 7, month = 30, 1)
    per_type_dur <- vapply(bcfg$grower_params,
                           function(g) g$pesticide_duration %||% bcfg$pesticide_duration,
                           numeric(1))
    max_dur <- max(bcfg$pesticide_duration, per_type_dur)
    buffer_end <- max(dd) + max_dur * step_days
    if (buffer_end > ed)
      return(fail(sprintf(
        paste("last decision_date (%s) + max pesticide_duration (%d %s steps) extends to %s,",
              "past end_date (%s); shorten the duration, move the last decision",
              "earlier, or extend end_date"),
        max(dd), max_dur, time_step %||% "day", buffer_end, ed)))
  }

  list(checks_passed = TRUE, failed_check = NULL)
}

# ---------------------------------------------------------------------------

#' Compute total epidemic size
#'
#' @param infected_raster SpatRaster of infected host counts.
#' @param host_raster SpatRaster of total host counts.
#' @param area Logical. If TRUE returns area in hectares; if FALSE (default)
#'   returns number of infected cells.
#'
#' @return Numeric scalar.
#' @importFrom terra values res
#' @export
compute_epidemic_size <- function(infected_raster, host_raster, area = FALSE) {
  inf_vals  <- terra::values(infected_raster, mat = FALSE)
  host_vals <- terra::values(host_raster, mat = FALSE)
  n_cells   <- sum(inf_vals > 0 & !is.na(host_vals) & host_vals > 0, na.rm = TRUE)
  if (!area) return(n_cells)
  cell_ha <- prod(terra::res(infected_raster)) / 1e4
  n_cells * cell_ha
}

# ---------------------------------------------------------------------------

#' Compute outbreak probability across stochastic replicates
#'
#' @param infected_list List of infected-host matrices (one per replicate),
#'   as returned in \code{$host_pools[[1]]$infected[[last_step]]} by
#'   \code{pops_behavior_simulate()}.
#' @param host_raster SpatRaster of total host counts.
#' @param threshold_fraction Numeric in (0, 1]. Minimum fraction of host area
#'   that must be infected to count as an "outbreak".
#'
#' @return Numeric in [0, 1]: proportion of replicates classed as outbreaks.
#' @importFrom terra values
#' @export
compute_outbreak_probability <- function(infected_list, host_raster,
                                         threshold_fraction = 0.1) {
  host_vals   <- terra::values(host_raster, mat = FALSE)
  total_host  <- sum(host_vals > 0, na.rm = TRUE)
  threshold_n <- threshold_fraction * total_host

  is_outbreak <- sapply(infected_list, function(m) {
    sum(as.vector(m) > 0, na.rm = TRUE) >= threshold_n
  })
  mean(is_outbreak)
}

# ---------------------------------------------------------------------------

#' Build a grower_params list suitable for pops_model_cpp()
#'
#' Converts the behavior config's grower_params block into the ordered list
#' of lists expected by the C++ interface (early_adopter first, then
#' responsive, then non_adopter).
#'
#' @param behavior_config Behavior config list (output of
#'   \code{behavior_configuration()}).
#'
#' @return List of three named lists, each with fields
#'   \code{willingness_to_treat}, \code{decision_threshold},
#'   \code{perception_radius}, \code{detection_probability},
#'   \code{treatment_efficacy}, and \code{pesticide_duration}. The last two
#'   may be set per type in \code{grower_params}; when omitted they fall back
#'   to the global \code{treatment_efficacy} / \code{pesticide_duration}.
#' @export
build_grower_params_list <- function(behavior_config) {
  bcfg <- if (!is.null(behavior_config$behavior)) behavior_config$behavior
          else behavior_config

  efficacy <- bcfg$treatment_efficacy
  duration <- bcfg$pesticide_duration

  lapply(c("early_adopter", "responsive", "non_adopter"), function(tp) {
    gp <- bcfg$grower_params[[tp]]
    list(
      willingness_to_treat  = gp$willingness_to_treat,
      decision_threshold    = gp$decision_threshold %||% 0.05,
      perception_radius     = as.integer(gp$perception_radius),
      detection_probability = gp$detection_probability,
      # per-type overrides, falling back to the global behavior-block values
      treatment_efficacy    = gp$treatment_efficacy %||% efficacy,
      pesticide_duration    = as.integer(gp$pesticide_duration %||% duration)
    )
  })
}

# NULL-coalescing helper (not exported)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------

#' Rasterize county-level adoption rates to host raster resolution
#'
#' Converts a SpatVector of county polygons with adoption rate attributes
#' to a SpatRaster at the same extent and resolution as the host raster.
#' Each cell receives the adoption rate from the county it falls in.
#'
#' @param county_polygons SpatVector (vect()) with polygon geometries.
#' @param adoption_field Character. Name of the column in county_polygons
#'   containing adoption rates (numeric, 0-1).
#' @param host_raster SpatRaster. Target for extent and resolution.
#' @param smooth Logical. If TRUE, apply spatial smoothing to the rasterized
#'   adoption rates (optional, default FALSE).
#' @param smooth_range Numeric. Range (in map units) for spatial smoothing.
#'   Only used if smooth = TRUE.
#'
#' @return SpatRaster with adoption rates, same extent and resolution as
#'   host_raster. NA cells are those outside all polygons.
#'
#' @importFrom terra rasterize res ext
#' @export
rasterize_county_adoption <- function(county_polygons, adoption_field,
                                      host_raster,
                                      smooth = FALSE, smooth_range = NULL) {
  # Rasterize: each cell gets the value from the polygon it intersects
  adoption_raster <- terra::rasterize(
    county_polygons, host_raster,
    field = adoption_field,
    fun = "mean"  # in case a cell straddles boundaries, take mean
  )

  if (smooth && !is.null(smooth_range)) {
    # Optional spatial smoothing (Gaussian focal mean)
    # For now, a stub; more sophisticated smoothing can be added later
    # w <- terra::focalMat(adoption_raster, smooth_range, type = "Gauss")
    # adoption_raster <- terra::focal(adoption_raster, w = w, na.rm = TRUE)
  }

  adoption_raster
}

# ---------------------------------------------------------------------------

#' Disaggregate county adoption rates by management unit
#'
#' Given a rasterized county adoption rate and a management unit raster,
#' computes the mean adoption rate per management unit. This is used as the
#' probability for type assignment when spatial_structure = "empirical".
#'
#' @param unit_raster SpatRaster of management unit IDs.
#' @param adoption_raster SpatRaster of adoption rates (0-1 per county).
#'
#' @return Named numeric vector: names are unit IDs, values are mean
#'   adoption rates for that unit (weighted by host cell count).
#'
#' @importFrom terra values mask
#' @export
unit_adoption_from_raster <- function(unit_raster, adoption_raster) {
  # Mask adoption rates to host cells only (where unit_raster is not NA)
  adoption_masked <- terra::mask(adoption_raster, unit_raster,
                                  maskvalues = NA)

  # Get cell-level data
  unit_vals <- as.vector(terra::values(unit_raster, mat = FALSE))
  adopt_vals <- as.vector(terra::values(adoption_masked, mat = FALSE))

  # Compute mean adoption per unit
  unit_ids <- unique(unit_vals[!is.na(unit_vals)])
  unit_adoption <- setNames(
    sapply(unit_ids, function(uid) {
      mean(adopt_vals[unit_vals == uid], na.rm = TRUE)
    }),
    as.character(unit_ids)
  )

  unit_adoption
}
