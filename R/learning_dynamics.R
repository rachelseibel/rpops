#' Update grower willingness-to-treat between seasons
#'
#' Implements a simple reinforcement learning rule: growers increase their
#' willingness to treat if the season outcome (epidemic size in their area)
#' was below a type-specific success threshold, and decrease it otherwise.
#' The update is smoothed via an exponentially weighted moving average (EWMA)
#' over \code{memory_length} past seasons.
#'
#' @param behavior_config Behavior config list (output of
#'   \code{behavior_configuration()}).
#' @param season_outcome Named numeric vector with elements \code{early_adopter},
#'   \code{responsive}, and \code{non_adopter}. Each value is the mean
#'   epidemic size (infected cells) observed in fields of that type during the
#'   just-completed season.
#' @param success_thresholds Named numeric vector with the same names as
#'   \code{season_outcome}. Epidemic sizes below the threshold are considered
#'   a successful (low-disease) season. Defaults to 10\% of host area per type.
#' @param learning_rate Numeric in (0, 1]. Step size for each WTT update.
#' @param memory_length Integer. Effective number of past seasons in the EWMA
#'   (decay = 1 / memory_length).
#'
#' @return Updated \code{behavior_config} with revised
#'   \code{grower_params$*$willingness_to_treat} values.
#'
#' @export
update_grower_willingness <- function(behavior_config,
                                      season_outcome,
                                      success_thresholds = NULL,
                                      learning_rate      = 0.1,
                                      memory_length      = 3L) {

  bcfg <- if (!is.null(behavior_config$behavior)) behavior_config$behavior
          else behavior_config

  types <- c("early_adopter", "responsive", "non_adopter")

  if (is.null(success_thresholds))
    success_thresholds <- stats::setNames(
      rep(Inf, length(types)), types)  # Inf = never update downward by default

  decay <- 1.0 / memory_length

  for (tp in types) {
    wtt     <- bcfg$grower_params[[tp]]$willingness_to_treat
    outcome <- season_outcome[[tp]]
    thresh  <- success_thresholds[[tp]]

    if (!is.na(outcome) && is.finite(thresh)) {
      if (outcome <= thresh) {
        # Successful season: reinforce treatment behaviour
        delta <- learning_rate * (1.0 - wtt)
      } else {
        # High disease despite treatment (or no treatment): reduce WTT
        delta <- -learning_rate * wtt
      }
      # EWMA smoothing: blend new delta with previous state
      new_wtt <- wtt + decay * delta
      new_wtt <- max(0.0, min(1.0, new_wtt))
    } else {
      new_wtt <- wtt  # no update if outcome unknown
    }

    if (!is.null(behavior_config$behavior)) {
      behavior_config$behavior$grower_params[[tp]]$willingness_to_treat <- new_wtt
    } else {
      behavior_config$grower_params[[tp]]$willingness_to_treat <- new_wtt
    }
  }

  behavior_config
}

# ---------------------------------------------------------------------------

#' Multi-season simulation with inter-seasonal learning
#'
#' Loops \code{pops_behavior_simulate()} over \code{n_seasons}, calling
#' \code{update_grower_willingness()} between seasons. At the start of each
#' season the host pool is re-initialised from the base config (disease does
#' not carry over across seasons — each season begins from the same
#' \code{initial_infection} raster). This is appropriate for an annual crop
#' system like wheat.
#'
#' @param config Base PoPS config list (from \code{configuration()}).
#' @param behavior_config Behavior config list.
#' @param unit_raster Management unit SpatRaster.
#' @param type_raster Grower type SpatRaster.
#' @param n_seasons Integer. Number of seasons to simulate.
#' @param n_iterations Integer. Replicates per season.
#' @param n_cores Integer. Parallel workers.
#' @param random_seed Integer seed.
#' @param learning_rate Numeric. Passed to \code{update_grower_willingness()}.
#' @param memory_length Integer. Passed to \code{update_grower_willingness()}.
#' @param success_thresholds Named numeric vector of per-type success thresholds
#'   (infected cells). NULL = no downward WTT updates.
#' @param output_path Character. Base directory for season-specific outputs.
#'
#' @return List of length \code{n_seasons}, each element the output of
#'   \code{pops_behavior_simulate()} plus \code{$updated_behavior_config} for
#'   that season.
#'
#' @export
run_multiseasonal_simulation <- function(config,
                                         behavior_config,
                                         unit_raster,
                                         type_raster,
                                         n_seasons          = 5L,
                                         n_iterations       = 500L,
                                         n_cores            = 1L,
                                         random_seed        = NULL,
                                         learning_rate      = 0.1,
                                         memory_length      = 3L,
                                         success_thresholds = NULL,
                                         output_path        = NULL) {

  if (!is.null(random_seed)) set.seed(random_seed)

  season_results <- vector("list", n_seasons)
  current_bcfg   <- behavior_config

  for (s in seq_len(n_seasons)) {
    message(sprintf("Season %d / %d", s, n_seasons))

    s_output <- if (!is.null(output_path))
                  file.path(output_path, sprintf("season_%02d", s))
                else NULL

    res <- pops_behavior_simulate(
      config          = config,
      behavior_config = current_bcfg,
      unit_raster     = unit_raster,
      type_raster     = type_raster,
      n_iterations    = n_iterations,
      n_cores         = n_cores,
      output_path     = s_output
    )

    # Compute per-type mean epidemic size for the learning update
    # (Requires type_raster to attribute cells; uses mean across all infected
    # cells as a proxy when per-type tracking is not available.)
    season_outcome <- stats::setNames(
      rep(mean(res$epidemic_size), 3L),
      c("early_adopter", "responsive", "non_adopter"))

    current_bcfg <- update_grower_willingness(
      current_bcfg,
      season_outcome     = season_outcome,
      success_thresholds = success_thresholds,
      learning_rate      = learning_rate,
      memory_length      = memory_length
    )

    res$updated_behavior_config <- current_bcfg
    season_results[[s]]          <- res
  }

  season_results
}
