# Synthetic vertical-profile experiments

Synthetic (OSSE) retrievals of H2O, CH4 and CO2 vertical profiles from a DCS spectrum. Each
experiment is a configuration of the SpectralFits **ProfileOSSE** app (in `configs/`). The app
builds a layered atmosphere and simulates a noisy spectrum on the wavenumbers of a real DCS
measurement (15-minute average 10 of `../data/20160921.h5`). It then retrieves the layer columns
with an optimal-estimation solver and writes a run folder. The scripts in this folder only
analyse those run folders or build cross-section tables.

## How to run

Run everything from this directory, with the sibling SpectralFit.jl checkout as the project.

To run one configuration (all of its variants):

```bash
julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/us_standard_co2.toml
```

To run only some variants:

```bash
julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/us_standard_co2.toml profile_tables isobar_tables
```

To check a configuration without running it (this loads the cross sections):

```bash
julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/crosslist_co2.toml --dry-run
```

To run every configuration and then the analysis (several hours, most of it line-by-line):

```bash
./run_profiles.sh
```

From Julia: `SpectralFits.Apps.run_config("configs/legacy_co2.toml"; variants = ["isobar"])`
returns one `ProfileOSSESummary` per variant.

## Configurations

`base.toml` holds the settings every experiment shares: the template record, template noise,
seed 1, a 0.002 cm⁻¹ fine grid, 5 Legendre baseline coefficients with σ 1, prior nugget 1e-6,
Levenberg-Marquardt, and `output_root = "../runs"`. `us_standard.toml` holds the US Standard
atmosphere. Neither file is run on its own. Experiment names are `<experiment>_<variant>`.

| Config | Variants | What it is |
|---|---|---|
| `us_standard_ch4.toml` | `profile`, `isotherm`, `isobar`, and the same three with `_tables` | Former `run_profiles.jl` CH4 runs. US Standard 1976 atmosphere, 35 layers 0–10 km. H2O is 0.006 exp(−z/2 km), CH4 2000→1900 ppb, CO2 410→400 ppm, ideal-gas layer columns. Window 6050–6120 cm⁻¹ in log space, γ 40, template noise. Prior correlated in height (H2O 0.01 km, CH4 and CO2 0.5 km). HITRAN 2008 CH4, HITRAN 2020 H2O and CO2. `_tables` variants use `../spectra/profile_tables/CH4_<gas>.jld2`. |
| `us_standard_co2.toml` | same six | Former `run_profiles.jl` CO2 runs. Same atmosphere. Window 6180–6260 cm⁻¹ in linear space, γ 50. HITRAN 2020 H2O and CO2; `_tables` variants use `CO2_<gas>.jld2`. |
| `legacy_ch4.toml` | `profile`, `isotherm`, `isobar` | Former `ch4_profile.jl`, `CH4_isotherm_test.jl` and `CH4_isobar_test.jl`, with the hand-port fixes. p = 1000 hPa exp(−z/8.5 km), T = 270 K − 6.5 K/km z, H2O 0.005→0.001, SNR 2000, log space. Prior correlated in pressure, 50 hPa (isotherm: diagonal). γ is 40 / 80 / 75. |
| `legacy_co2.toml` | `profile`, `isotherm`, `isobar` | Former `CO2_profile.jl`, `CO2_isotherm_test.jl` and `CO2_isobar_test.jl`, linear space. `profile`: T₀ 290 K, SNR 2000, prior correlated 50 hPa. `isotherm`: layers at 300 K, columns from T₀ 280 K, SNR 1500, diagonal prior. `isobar`: T₀ 290 K, columns from T₀ 270 K (`column_temperature`), CO2 σ 80 ppm. |
| `ch4_profile2.toml` | none | Former `CO2_profile2.jl`, a CH4-window run despite its name. 44 layers, hydrostatic dry-air columns, T₀ 300 K, H2O 0.001→0.01, CH4 2000→1800 ppb. HITRAN 2020 for all gases, 100 baseline coefficients, diagonal prior, Gauss-Newton, linear space, SNR 2000. |
| `crosslist_ch4.toml` | `truth-HIT08_retr-{HIT08,HIT16,HIT20,TCCON}_{profile,isotherm,isobar}` (12) | Cross-line-list runs on `us_standard_ch4`: the truth spectrum uses HITRAN 2008 CH4 (`truth_cross_section`), and the retrieval uses another CH4 list. |
| `crosslist_co2.toml` | `truth-ABSCO51_retr-{ABSCO51,HIT16,HIT20,TCCON}_*` and `truth-HIT20_retr-{HIT20,HIT16,TCCON}_*` (21) | The same on `us_standard_co2`. The ABSCO 5.1 truth is the dry table from `make_absco_table.jl`. |

The modes collapse the atmosphere as follows:

- `profile` keeps T(z) and p(z).
- `isotherm` sets every layer to the column-weighted mean T, or to `isotherm_temperature_K`.
- `isobar` sets every layer to the column-weighted mean p.

Layer columns always come from the uncollapsed profile. A correlated prior in `pressure_hPa`
uses the profile pressures. The full key reference is the docstring of
`profile_osse_config` (`SpectralFit.jl/src/apps/profile_osse/profile_config.jl`).

## Outputs

Each run writes `runs/<yyyy-mm-dd>_<experiment>_<SpectralFits git hash>[_k]/`, containing:

- `configuration.toml`: the fully merged configuration, with absolute paths.
- `summary.txt`: status, χ²/dof, noise σ and SNR, plus one CSV row per gas (`gas,dfs,column_error_percent,column_sigma_percent`) and the solver report.
- `results/<experiment>_results.jld2`: `load_results` format. `results[1]` is the retrieval (`x`, `K`, `averaging_kernel`, `dfs_per_element`, `posterior_covariance`, ...). The attributes hold `x_true` and `prior_sigma` (with `_labels`), `height_m`, layer `pressure_hPa` and `temperature_K`, `dry_air_column_cm2`, `profile_pressure_hPa`, `mode`, `noise_sigma`, `snr`, `seed` and `target`.

`find_run("runs", experiment)` returns the latest folder of an experiment, and `load_run` in
`util.jl` loads it. The `.jld2` files in this folder (`CH4_profile.jld2`, ...,
`*_tables.jld2`) and `linelist_results/` are outputs of the v0.2 hand-port scripts. They are
kept for comparison; nothing writes to them any more.

## Scripts

| Script | Purpose |
|---|---|
| `util.jl` | Shared helpers: `profile_config(name, variant)`, `load_run(experiment)`, `retrieved_vmr`, `true_vmr`, `column_ak` (v0.1 convention (H'A)./H'), `layer_dfs`, `rms`. |
| `load_profile.jl` | Example of loading one run and printing its target profile. Usage: `julia --project=../../SpectralFit.jl -i load_profile.jl us_standard_co2_isobar_tables`. |
| `error_analysis.jl` | Column averaging kernels and retrieved profiles of the six runs of a set. Usage: `error_analysis.jl [legacy\|us_standard] [suffix]`, for example `us_standard _tables`. Writes `runs/error_analysis_<set><suffix>.csv`, and the two PDF figures when Plots is installed. |
| `linelist_errors.jl` | Analysis of the crosslist runs. Writes `runs/linelist_analysis/linelist_summary_<stem>.csv` and `linelist_profiles_<stem>.csv` (same columns as the hand-port's). Δ columns are relative to the control, the first retrieval of each truth. Missing runs are skipped, and the stem gets the modes found. |
| `linelist_tables_md.jl` | Merges the summary CSVs and prints markdown tables (`julia linelist_tables_md.jl [dir]`). |
| `make_profile_tables.jl` | Rebuilds `../spectra/profile_tables/*.jld2` (the smallest p/T tables within 0.1σ), reading atmosphere, windows and line lists from the line-by-line variants of `us_standard_*.toml`. |
| `make_absco_table.jl` | Builds `../spectra/linelist_tables/ABSCO51_CO2_dry.jld2` from the OCO ABSCO 5.1 file. |
| `run_profiles.sh` | Runs the configs (all, or the names given) and then the analysis scripts. |
| `inversion.jl` | v0.1 code, not ported. |

## Changes in results, v0.1 → v0.2

- **DFS and posterior.** SpectralFits now computes `posterior_covariance`, `averaging_kernel` and `dfs_per_element` with a Jacobi-scaled inverse of the normal matrix. These states mix layer columns near 1e22 molecules/cm² with baseline coefficients near 1, so the unscaled normal matrix has a condition number near 1e54. The DFS of v0.1, and of early v0.2 before this fix, could be wrong by up to about 1.5. The DOF values below match the σ-scaled diagnostics of the hand-port. Old v0.1 DOF numbers should not be compared with them.
- **Prior nugget.** Gaussian-kernel prior correlations on these layer spacings are numerically singular: for 50 hPa on layers 10–30 hPa apart the condition number is near 1e16. v0.1 inverted the matrix anyway, so its Sₐ⁻¹ was mostly round-off. v0.2 adds `nugget × σᵢ²` to the diagonal (default 1e-6, set in `base.toml`).
- **US Standard constants.** The app's US Standard pressure uses the package constants (dry-air molar mass 28.9647 g/mol and R = N_A k_B). The hand-port used 28.9644 g/mol and R = 8.31446. Layer pressures differ in about the fifth digit.
- **Legacy fixes** (kept from the hand-port, listed in the headers of `legacy_*.toml` and `ch4_profile2.toml`): the H2O profile direction, the CH4 isobar H2O prior, noise σ = mean intensity / SNR, line-by-line everywhere, and dry-air hydrostatic columns.

The app runs of 2026-10-20 (SpectralFits a93f406 plus uncommitted changes) compared with the
hand-port's `run_profiles.jl tables` runs, which used the hand-written σ-scaled diagnostics:

```csv
window,mode,column_error_pct_app,column_error_pct_handport,dof_app,dof_handport,chi2_obs_per_dof_app
CO2,profile,-0.0277,-0.028,3.979,3.98,0.991
CO2,isotherm,-0.0230,-0.023,4.155,4.15,0.992
CO2,isobar,-0.0135,-0.013,2.798,2.80,0.992
CH4,profile,-0.0018,-0.002,3.268,3.27,0.996
CH4,isotherm,+0.0005,+0.001,3.447,3.45,0.996
CH4,isobar,-0.0082,-0.008,2.274,2.27,0.995
```

`legacy_co2_isobar` converged in 3 iterations, with χ²/dof 0.991, CO2 DOF 2.93 and column
error −0.006 %.

For the crosslist runs with the HIT20 CO2 truth in the isobar mode, the column errors match the
hand-port's `linelist_results/` (HIT20 −0.014 vs −0.010 %, HIT16 +0.749 vs +0.740 %, TCCON
+0.865 vs +0.865 %). The DOF match as well (2.79). The shape RMS against truth differs more
(control 1.32 vs 0.78 %, HIT16 6.0 vs 11.8 %, TCCON 15.4 vs 17.7 %). With about 2.8 DOF over
35 layers the profile shape is weakly constrained, so it is sensitive to the small changes above
and to where the solver stops (Rodgers tolerance 0.01).
