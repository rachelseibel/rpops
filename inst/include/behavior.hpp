/*
 * PoPS behavior module — grower fungicide decision action
 *
 * Implements a spatially explicit grower decision model that fires at
 * scheduled spray windows (T1, T2, T3). At each decision step each
 * management unit observes local disease prevalence, draws a stochastic
 * treatment decision conditioned on its grower type, and injects the
 * resulting treatment map directly into the Treatments object so it is
 * applied within the same simulation step.
 *
 * Three grower types are supported:
 *   type_code 1 — early adopter  (prophylactic + reactive)
 *   type_code 2 — responsive     (reactive only)
 *   type_code 3 — non-adopter    (near-zero baseline)
 *
 * Authors: Rachel Seibel
 * License: GPL-3
 */

#ifndef POPS_BEHAVIOR_HPP
#define POPS_BEHAVIOR_HPP

#include "scheduling.hpp"
#include "date.hpp"

#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <utility>
#include <stdexcept>

namespace pops {

// ---------------------------------------------------------------------------
// Per-type parameter struct
// ---------------------------------------------------------------------------

/**
 * @brief Parameters for one grower type.
 *
 * All probability fields must be in [0, 1].
 * perception_radius is in raster cells (integer).
 */
struct GrowerTypeParams {
    double willingness_to_treat;    ///< Base treatment probability
    double decision_threshold;      ///< Local prevalence above which treatment is triggered
    int    perception_radius;       ///< Focal window radius (cells)
    double detection_probability;   ///< Probability of detecting one infected plant
    double treatment_efficacy;      ///< Fraction of infection removed on treatment
    int    pesticide_duration;      ///< Steps this type's treatment remains effective
};

// ---------------------------------------------------------------------------
// Prevalence perception mode
// ---------------------------------------------------------------------------

/**
 * @brief How a grower perceives local disease prevalence.
 *
 * Window: detected infected / total host within a circular window of
 *   perception_radius around the unit centroid. Sweeps in neighbouring-unit
 *   host and non-host buffer, so a larger radius DILUTES perceived prevalence.
 * Unit:   detected infected / total host over the grower's OWN management-unit
 *   cells only (host-only denominator). Prevalence reflects the grower's own
 *   field and is independent of perception_radius and unit spacing.
 */
enum class PrevalenceMode { Window = 0, Unit = 1 };

// ---------------------------------------------------------------------------
// GrowerDecisionAction
// ---------------------------------------------------------------------------

/**
 * @brief Grower decision action fired at scheduled spray windows.
 *
 * @tparam IntegerRaster  Raster type for integer counts (e.g. Rcpp::IntegerMatrix)
 * @tparam FloatRaster    Raster type for floating-point maps  (e.g. Rcpp::NumericMatrix)
 *
 * Construct once before the simulation loop, then call should_act(step) each
 * iteration and, when true, call fill_treatment_map() + retrieve date_for_step()
 * to pass the treatment to Treatments::add_treatment().
 */
template<typename IntegerRaster, typename FloatRaster>
class GrowerDecisionAction
{
public:
    /**
     * @brief Construct and precompute unit metadata.
     *
     * @param unit_map     Integer raster of management unit IDs (positive = host unit,
     *                     0/negative = non-host buffer).
     * @param type_map     Integer raster of type codes (1 = early adopter,
     *                     2 = responsive, 3 = non-adopter). Must be co-registered
     *                     with unit_map.
     * @param type_params  Parameters for each type, ordered early_adopter (index 0),
     *                     responsive (index 1), non_adopter (index 2).
     * @param decision_date_strings  ISO-8601 date strings (YYYY-MM-DD) for each
     *                     decision window (T1, T2, T3, …).
     * @param scheduler    Scheduler used to convert dates to step indices.
     * @param rows         Raster row count.
     * @param cols         Raster column count.
     * @param pesticide_duration  Number of simulation steps the treatment lasts.
     * @param behavior_seed  Seed for the behavior-module RNG (kept separate from
     *                     the main dispersal/establishment generators).
     */
    GrowerDecisionAction(
        const IntegerRaster& unit_map,
        const IntegerRaster& type_map,
        std::vector<GrowerTypeParams> type_params,
        const std::vector<std::string>& decision_date_strings,
        const Scheduler& scheduler,
        int rows,
        int cols,
        int pesticide_duration,
        unsigned behavior_seed,
        PrevalenceMode prevalence_mode = PrevalenceMode::Window)
        : type_params_(std::move(type_params)),
          rows_(rows),
          cols_(cols),
          pesticide_duration_(pesticide_duration),
          prevalence_mode_(prevalence_mode),
          generator_(behavior_seed)
    {
        for (const auto& ds : decision_date_strings) {
            pops::Date d(ds);
            unsigned step = scheduler.schedule_action_date(d);
            decision_steps_.push_back(step);
            step_to_date_.emplace(step, d);
        }
        std::sort(decision_steps_.begin(), decision_steps_.end());
        precompute_units(unit_map, type_map);
    }

    /** Returns true if a grower decision should fire at this simulation step. */
    bool should_act(unsigned step) const
    {
        return std::binary_search(
            decision_steps_.begin(), decision_steps_.end(), step);
    }

    /** Returns the Date associated with a decision step (for treatments.add_treatment). */
    pops::Date date_for_step(unsigned step) const
    {
        auto it = step_to_date_.find(step);
        if (it == step_to_date_.end())
            throw std::out_of_range(
                "GrowerDecisionAction::date_for_step: step not a decision step");
        return it->second;
    }

    /** Fallback treatment duration (used when a type has none set). */
    int pesticide_duration() const { return pesticide_duration_; }

    /**
     * @brief Compute treatment decisions, bucketed by pesticide duration.
     *
     * For each management unit: observe local prevalence and draw a stochastic
     * treatment decision. Treating units contribute their cells (set to the
     * unit type's treatment_efficacy) to the treatment map for that type's
     * pesticide_duration. Because PoPS schedules a treatment's decay by a
     * single duration, units with different durations must be applied as
     * separate treatments -- hence one map per distinct duration.
     *
     * The decision for each unit is drawn exactly once here, so RNG
     * consumption is independent of how many distinct durations exist.
     *
     * @param infected      Current infected-host raster (first host species).
     * @param total_hosts   Total host raster.
     * @return Map from pesticide_duration (steps) to a FloatRaster whose
     *         treated cells hold the relevant unit's treatment_efficacy. Empty
     *         if no unit decided to treat.
     */
    std::map<int, FloatRaster> fill_treatment_maps(
        const IntegerRaster& infected,
        const IntegerRaster& total_hosts)
    {
        std::map<int, FloatRaster> maps_by_duration;

        for (auto& kv : units_) {
            UnitInfo& ui = kv.second;
            const GrowerTypeParams& params = type_params_[ui.type_idx];

            double prevalence =
                prevalence_mode_ == PrevalenceMode::Unit
                    ? compute_prevalence_unit(
                          infected, total_hosts, ui,
                          params.detection_probability)
                    : compute_prevalence(
                          infected, total_hosts,
                          ui.centroid_row, ui.centroid_col,
                          params.perception_radius,
                          params.detection_probability);

            if (draw_decision(prevalence, ui.type_idx, params)) {
                int dur = params.pesticide_duration > 0
                              ? params.pesticide_duration
                              : pesticide_duration_;
                auto it = maps_by_duration.find(dur);
                if (it == maps_by_duration.end())
                    it = maps_by_duration.emplace(
                             dur, FloatRaster(rows_, cols_)).first;
                for (const auto& cell : ui.cells) {
                    it->second(cell.first, cell.second) =
                        params.treatment_efficacy;
                }
            }
        }
        return maps_by_duration;
    }

private:
    // -----------------------------------------------------------------------
    // Internal unit metadata
    // -----------------------------------------------------------------------

    struct UnitInfo {
        int centroid_row;
        int centroid_col;
        int type_idx;  ///< 0-based index into type_params_
        std::vector<std::pair<int, int>> cells;
    };

    std::unordered_map<int, UnitInfo>       units_;
    std::unordered_map<unsigned, pops::Date> step_to_date_;
    std::vector<unsigned>                   decision_steps_;
    std::vector<GrowerTypeParams>           type_params_;
    int                                     rows_;
    int                                     cols_;
    int                                     pesticide_duration_;
    PrevalenceMode                          prevalence_mode_;
    std::mt19937                            generator_;

    // -----------------------------------------------------------------------
    // Precomputation
    // -----------------------------------------------------------------------

    void precompute_units(const IntegerRaster& unit_map,
                          const IntegerRaster& type_map)
    {
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                int uid       = unit_map(i, j);
                int type_code = type_map(i, j);
                if (uid <= 0 || type_code <= 0) continue;

                auto& ui      = units_[uid];
                ui.type_idx   = type_code - 1;  // 1-based → 0-based
                ui.cells.emplace_back(i, j);
            }
        }

        for (auto& kv : units_) {
            UnitInfo& ui = kv.second;
            if (ui.cells.empty()) continue;
            double sum_r = 0, sum_c = 0;
            for (const auto& cell : ui.cells) {
                sum_r += cell.first;
                sum_c += cell.second;
            }
            double n = static_cast<double>(ui.cells.size());
            ui.centroid_row = static_cast<int>(std::round(sum_r / n));
            ui.centroid_col = static_cast<int>(std::round(sum_c / n));
        }
    }

    // -----------------------------------------------------------------------
    // Prevalence estimation
    // -----------------------------------------------------------------------

    /**
     * @brief Compute perceived local disease prevalence for one unit.
     *
     * Iterates over a circular focal window of radius `perception_radius`
     * centred on the unit centroid. Each infected plant is independently
     * detected with probability `detection_probability` (Binomial draw).
     * Prevalence = sum(detected infected plants) / sum(total host plants).
     */
    double compute_prevalence(
        const IntegerRaster& infected,
        const IntegerRaster& total_hosts,
        int centroid_row,
        int centroid_col,
        int perception_radius,
        double detection_probability)
    {
        int total_infected_detected = 0;
        int total_host_plants       = 0;

        int r_min = std::max(0,         centroid_row - perception_radius);
        int r_max = std::min(rows_ - 1, centroid_row + perception_radius);
        int c_min = std::max(0,         centroid_col - perception_radius);
        int c_max = std::min(cols_ - 1, centroid_col + perception_radius);

        for (int i = r_min; i <= r_max; ++i) {
            for (int j = c_min; j <= c_max; ++j) {
                int dr = i - centroid_row;
                int dc = j - centroid_col;
                if (dr * dr + dc * dc > perception_radius * perception_radius)
                    continue;

                int hosts = total_hosts(i, j);
                if (hosts <= 0) continue;
                total_host_plants += hosts;

                int inf = infected(i, j);
                if (inf > 0 && detection_probability > 0.0) {
                    int n_detect = std::min(inf, hosts);
                    std::binomial_distribution<int> bd(n_detect, detection_probability);
                    total_infected_detected += bd(generator_);
                }
            }
        }

        if (total_host_plants == 0) return 0.0;
        return static_cast<double>(total_infected_detected)
               / static_cast<double>(total_host_plants);
    }

    /**
     * @brief Compute perceived prevalence over a unit's OWN cells only.
     *
     * Host-only denominator: detected infected / total host across the cells
     * belonging to this management unit. Independent of perception_radius and
     * of how units are spaced, so it reflects the grower's own field rather
     * than a geometric neighbourhood that mixes units and buffer.
     */
    double compute_prevalence_unit(
        const IntegerRaster& infected,
        const IntegerRaster& total_hosts,
        const UnitInfo& ui,
        double detection_probability)
    {
        int total_infected_detected = 0;
        int total_host_plants       = 0;

        for (const auto& cell : ui.cells) {
            int hosts = total_hosts(cell.first, cell.second);
            if (hosts <= 0) continue;
            total_host_plants += hosts;

            int inf = infected(cell.first, cell.second);
            if (inf > 0 && detection_probability > 0.0) {
                int n_detect = std::min(inf, hosts);
                std::binomial_distribution<int> bd(n_detect, detection_probability);
                total_infected_detected += bd(generator_);
            }
        }

        if (total_host_plants == 0) return 0.0;
        return static_cast<double>(total_infected_detected)
               / static_cast<double>(total_host_plants);
    }

    // -----------------------------------------------------------------------
    // Treatment decision
    // -----------------------------------------------------------------------

    /**
     * @brief Stochastic treatment decision for one management unit.
     *
     * - Early adopter  (type_idx 0): treat if prevalence > threshold OR with
     *   base probability willingness_to_treat (prophylactic component).
     * - Responsive     (type_idx 1): treat if prevalence > threshold AND
     *   Bernoulli(willingness_to_treat) is true.
     * - Non-adopter    (type_idx 2): treat with near-zero Bernoulli baseline.
     */
    bool draw_decision(double prevalence,
                       int type_idx,
                       const GrowerTypeParams& params)
    {
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        switch (type_idx) {
            case 0:  // early adopter
                return (prevalence > params.decision_threshold)
                    || (unif(generator_) < params.willingness_to_treat);
            case 1:  // responsive
                return (prevalence > params.decision_threshold)
                    && (unif(generator_) < params.willingness_to_treat);
            default: // non-adopter
                return unif(generator_) < params.willingness_to_treat;
        }
    }
};

}  // namespace pops

#endif  // POPS_BEHAVIOR_HPP
