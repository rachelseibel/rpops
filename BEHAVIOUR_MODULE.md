# Grower behaviour module

This branch adds an endogenous grower decision-making layer to PoPS. It was
developed by Rachel Seibel for Chapter 5 of her PhD thesis (University of
Warwick, MathSys CDT), which applies it to Septoria tritici blotch in English
winter wheat.

It is kept on a branch of a fork so that the additions can be reviewed as a
diff against `ncsu-landscape-dynamics/rpops` without anything landing in the
upstream repository.

## What it does

Upstream PoPS applies treatment where and when a treatment raster says to.
This module generates that treatment map instead of requiring it as input.
The landscape is divided into management units, each unit is assigned a grower
type, and at each decision date a unit forms an imperfect perception of local
disease and draws a stochastic decision on whether to spray. The units that
decide to spray are collected into a treatment map, which is then handed to
the PoPS treatment module rather than removing infection directly. The
existing spread loop is unchanged when grower behaviour is switched off.

## Files added

| File | Purpose |
|---|---|
| `R/management_units.R` | Delineate management units from the host raster: contiguous patches, polygon field boundaries (e.g. CROME), or a regular grid. |
| `R/behavior_type_assignment.R` | Assign grower types to units under a chosen spatial structure: random, clustered via a Gaussian random field, graded between the two, or empirically constrained from a survey-derived adoption surface. |
| `R/behavior_helpers.R` | Parse and validate the `behavior:` block of a YAML config, build the per-type parameter list, and compute epidemic-size and outbreak-probability summaries. |
| `R/learning_dynamics.R` | Between-season updating of willingness-to-treat, and the multi-season driver. |

## Files modified

| File | Change |
|---|---|
| `R/configuration.R` | Accept and validate the `behavior:` configuration block. |
| `R/pops_model.R`, `R/pops_simulate.R`, `R/pops.r` | Route the generated treatment map into the existing treatment process. |
| `src/pops.cpp`, `src/RcppExports.cpp`, `R/RcppExports.R` | Expose the per-decision-date hooks the R layer needs. |
| `DESCRIPTION`, `NAMESPACE` | Package metadata and exports (see the note below). |

Nothing in the upstream spread, dispersal, or treatment machinery is replaced.

## Note on the package name

`DESCRIPTION` renames the package from `PoPS` to `PoPSbehaviour`, so that this
build can be installed alongside upstream PoPS without a collision. This is
deliberate but it does mean the fork is not a drop-in replacement: code that
does `library(PoPS)` will need `library(PoPSbehaviour)`. If the module is ever
proposed upstream, that rename should be dropped first, since it accounts for a
large share of the metadata diff and none of the behavioural logic.

Authorship in `DESCRIPTION` has been adjusted so that the maintainer field
points at the fork author rather than at the upstream maintainer. Upstream
authors are retained as authors. This is a convention for the fork only.

## Reviewing the changes

```bash
git remote add upstream https://github.com/ncsu-landscape-dynamics/rpops.git
git fetch upstream
git diff upstream/main...grower-behaviour
```

Or on GitHub, use "Compare across forks" with `ncsu-landscape-dynamics/rpops`
as the base and this branch as the head.

## Citing

If you use this module, please cite the PoPS model itself alongside the thesis
chapter:

- Jones, C. et al. (2021). *Iteratively forecasting biological invasions with
  PoPS and a little help from our friends.* Frontiers in Ecology and the
  Environment 19(8): 411–418.
- Petrasova, A. et al. (2020). *Fast simulation of a pest or pathogen spread
  with PoPS.*
