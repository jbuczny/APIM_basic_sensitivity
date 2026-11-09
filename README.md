# APIM_basic_sensitivity

Simulation-based power and sensitivity analysis for the Actor–Partner
Interdependence Mediation Model (APIMeM) with **two parallel mediators**
and distinguishable dyads, at fixed sample sizes.

For each sample size in `datasets` (here 130 and 174 dyads), the script
answers: *with this many dyads, how much power is there to detect the
actor, partner and indirect effects of a given size?* Each data set is
analysed on its own, as if it were the only sample.

## Model

Variables, for dyad members 1 and 2 (e.g. men and women):

| Role        | Member 1 | Member 2 |
|-------------|----------|----------|
| Predictor   | X1       | X2       |
| Mediator M  | M1       | M2       |
| Mediator W  | W1       | W2       |
| Outcome     | Y1       | Y2       |

The population and analysis model is the saturated distinguishable-dyad
APIMeM (Ledermann, Macho & Kenny, 2011): all actor and partner a paths
(X → M, X → W), b paths (M → Y, W → Y) and direct c paths (X → Y), with
correlated residuals within and between partners. The lavaan syntax is in
section 3 of the script. All variables have unit variance, so the
coefficients are standardized betas.

Simple indirect effects (actor–actor, partner–partner, actor–partner,
partner–actor, per mediator and outcome member) are tested with a Monte
Carlo confidence interval (Preacher & Selig, 2012); joint significance is
available as an option. Total indirect effects and M-versus-W contrasts
always use the Monte Carlo interval.

## Outputs

For each data set:

- **(A)** power at the hypothesized betas for all 20 paths, the 16 simple
  indirect effects, the total indirect effects and the M − W contrasts;
- **(B)** power curves: each focal simple indirect effect swept across a
  grid of betas, other paths held at their hypothesized values;
- the smallest beta reaching the target power (80%), by linear
  interpolation;
- a plot of the power curves, exported to PNG.

## Files

| File | Contents |
|------|----------|
| `apimem_power_2med_fixedN_parallel_multiN.R` | Analysis script |
| `apimem2_power_A.csv` | (A) Power and Monte Carlo SE at the hypothesized betas, all effects, both data sets |
| `apimem2_power_B.csv` | (B) Power by beta for each focal simple indirect effect, both data sets |
| `apimem2_power_by_dataset.csv` | (A) with the data sets side by side |
| `power_plots/apimem2_power_D1_N130.png` | Power curves, N = 130 dyads |
| `power_plots/apimem2_power_D2_N174.png` | Power curves, N = 174 dyads |

The results in this repository were produced with the settings in the
script as committed: 1,000 replications per scenario, 2,000 Monte Carlo
draws per confidence interval, α = .05, seed 2024.

## Requirements

R (≥ 4.0) with the packages `lavaan` and `future.apply`:

```r
install.packages(c("lavaan", "future.apply"))
```

## Running the analysis

1. Open `apimem_power_2med_fixedN_parallel_multiN.R`.
2. Set the working directory to the repository folder (in RStudio:
   *Session → Set Working Directory → To Source File Location*). The CSV
   files and the `power_plots/` folder are written there.
3. Adjust the settings in section 1 if needed (see below).
4. Source the script. With the default settings a run takes several
   minutes per data set.

### Main settings (section 1)

| Setting | Meaning |
|---------|---------|
| `datasets` | Named vector of sample sizes (number of dyads), one analysis each |
| `paths_hyp` | Hypothesized standardized betas for all 20 paths |
| `rxx`, `rmm`, `rww`, `rmw_in`, `rmw_x`, `ryy` | Predictor correlation and residual correlations |
| `beta_grid`, `focal_ind`, `sweep_component` | Grid and effects for the power curves (B) |
| `ind_test`, `mc_draws` | Test for simple indirect effects; Monte Carlo draws per interval |
| `n_reps`, `alpha`, `target`, `seed` | Replications per scenario, significance level, target power, seed |
| `parallel`, `n_workers`, `parallel_over` | Parallel processing |
| `export_png`, `out_dir`, `png_width`, `png_height`, `png_res` | Plot export |

## Computation and reproducibility

- **Estimation.** With `engine = "fast"` (default) each replication uses
  the closed-form maximum-likelihood solution, which for this saturated
  model is identical to fitting the lavaan model; this is checked against
  lavaan before the run starts. `engine = "lavaan"` fits the SEM in every
  replication.
- **Parallel processing** uses `future.apply`, with forked workers where
  available and separate R sessions otherwise (e.g. on Windows).
- **Reproducibility.** Replication *i* always starts from
  `set.seed(seed + i)`. Results are therefore identical across runs, for
  any number of workers, and with `parallel <- FALSE`. The same
  replications are used in every scenario (common random numbers), so
  differences between grid points reflect effect size rather than
  simulation noise.

## References

Ledermann, T., Macho, S., & Kenny, D. A. (2011). Assessing mediation in
dyadic data using the actor-partner interdependence model. *Structural
Equation Modeling, 18*(4), 595–612.
https://doi.org/10.1080/10705511.2011.607099

Preacher, K. J., & Selig, J. P. (2012). Advantages of Monte Carlo
confidence intervals for indirect effects. *Communication Methods and
Measures, 6*(2), 77–98. https://doi.org/10.1080/19312458.2012.679848
