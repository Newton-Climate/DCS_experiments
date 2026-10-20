# DCS retrievals, 2016 DCS A campaign

Retrievals of the 2016 DCS A open-path data (`../data/*.h5`) with the SpectralFits v0.2 batch
app. Each experiment is a TOML config in [configs/](configs/); the scripts here only analyse
the results.

## Running a config

From this directory, with the SpectralFit.jl checkout next to `DCS_experiments`:

```bash
julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/hit24.toml
```

Name variants after the config to run only those (`... configs/hit24.toml hit24_lm hit24_sdvlm`),
and add `--dry-run` to build and check every configuration without retrieving. From Julia:
`SpectralFits.Apps.run_config("configs/hit24.toml"; variants = ["hit24_lm"], dry_run = true)`.

[run_retrievals.sh](run_retrievals.sh) starts configs in the background (default
`hit08 hit16 hit20 tccon`; `./run_retrievals.sh hit24 merlin`), logging to `<config>.out`.

Every config extends [configs/campaign_2016.toml](configs/campaign_2016.toml): every `.h5` file
in `../data`, 15-minute averages, 4 worker processes, DCS A path 195017 cm, pressure scale
0.9938, shot noise from 6255.1–6255.4 cm⁻¹, 100 Legendre baseline coefficients, p and T
retrieved, the CH4 line masks (de-weighted by 1e-5), Gauss-Newton, and wide priors (v0.1 fitted
by unconstrained least squares). A variant's experiment is `<experiment>_<variant name>`.

## Results

`../data/runs/<yyyy-mm-dd>_<experiment>_<SpectralFits git hash>[_k]/`, with the configuration
(`configuration.txt`, `manifest.toml`) and `results/<input>_results.jld2` per input file. Find
the latest run with `find_run("../data/runs", "ch4_hit24")`, list experiments with
`run_experiments("../data/runs")`, read with `load_results(path).results`. Each result has
`x` (e.g. `r.x.CH4`, columns in molecules/cm²), `status`, `y`, `y_model`, `K`,
`posterior_covariance`, `chi2_obs_per_dof`, `rmse` (full-weight channels), `coordinate`
(wavenumbers, cm⁻¹), `time` (start of the average) and `metadata` (`machine_time`,
`first_sample_time`, `conditions` = the p, T and path the retrieval used).

## Configs and experiments

Windows: CH4 6050–6115 cm⁻¹, CO2 6180–6260 cm⁻¹ (v0.1 retrieved them separately, so a
two-window v0.1 script is one config with a variant per window). Tables are in `../spectra`.

| config | experiments | was | lookup tables needed | built by (`../spectra/tables/`) |
|---|---|---|---|---|
| [hit08.toml](configs/hit08.toml) | `hit08_ch4`, `hit08_co2` (masks off) | run_hit08.jl | tccon_H2O, tccon_HDO, hit20_12CO2, hit08_12CH4 | full_range.toml ¹ |
| [hit16.toml](configs/hit16.toml) | `hit16_ch4`, `hit16_co2` | run_hit16.jl ² | tccon_H2O, tccon_HDO, hit16_12CO2, hit16_12CH4 | full_range.toml ¹ |
| [hit20.toml](configs/hit20.toml) | `hit20_ch4`, `hit20_co2` | run_hit20.jl ² | tccon_H2O, tccon_HDO, hit20_12CO2, hit20_12CH4 | full_range.toml ¹ |
| [tccon.toml](configs/tccon.toml) | `tccon_ch4`, `tccon_co2` | run_tccon.jl ² | tccon_H2O, tccon_HDO, tccon_12CO2, tccon_12CH4 | full_range.toml ¹ |
| [ch4_broadening.toml](configs/ch4_broadening.toml) | `ch4_broadening_ch4`, `_co2` (H2O broadening of CH4, 0.34); `ch4_broadening_test_off`, `_test_on` (H2O + CH4 only, 6050–6120 cm⁻¹) | run_CH4broadening.jl, test_ch4_broadening.jl | as hit08; tests: tccon_H2O, hit20_12CH4 | full_range.toml ¹ |
| [oco.toml](configs/oco.toml) | `oco_co2` (ABSCO dry-air slice) | run_OCO.jl | tccon_H2O, tccon_HDO, hit08_12CH4; OCO_spectra.hdf (present) | full_range.toml ¹ |
| [oco_h2o.toml](configs/oco_h2o.toml) | `oco_h2o_co2` (ABSCO with H2O axis); `oco_h2o_test_absco`, `_test_hit20` | run_OCO_H2O.jl, test_OCO_tables.jl | tccon_H2O, tccon_HDO, hit08_12CH4 (tests: hit20_12CH4, hit20_12CO2); OCO_spectra.hdf | full_range.toml ¹ |
| [h2o_test.toml](configs/h2o_test.toml) | `h2o_oco_co2`, `h2o_tccon_co2`, `h2o_hit16_co2`, `h2o_hit20_co2` (CO2 window) | h2o_test/run_*.jl | OCO_H2O ³, tccon_H2O, hit16_H2O, hit20_H2O, tccon_HDO, hit08_12CH4; OCO_spectra.hdf | full_range.toml ¹ |
| [hit24.toml](configs/hit24.toml) | `ch4_hit24`, `ch4_hit24_lm`, `ch4_hit24_sdvlm`, `ch4_hit08`, and each `..._fixpT` | run_hit24.jl (CH4) | tccon_H2O_ch4win, tccon_HDO_ch4win, hit20_12CO2_ch4win, hit24_12CH4{,_lm,_sdvlm}_ch4win (present), hit08_12CH4_ch4win | hit24_ch4win.toml; line mixing: ../make_hit24_lm_lookup_tables.jl (branch) then hit24_lm_ch4win.toml |
| [hit24_co2.toml](configs/hit24_co2.toml) | `co2_oco`, `co2_hit24`, `co2_oco_fixpT`, `co2_hit24_fixpT` | run_hit24.jl (CO2) | tccon_H2O_co2win, tccon_HDO_co2win, hit08_12CH4_co2win, OCO_12CO2_co2win, hit24_12CO2_co2win | hit24_co2win.toml |
| [merlin.toml](configs/merlin.toml) | `r6_merlin`, `r6_hit08_narrow`, `r6_merlin_HT` (6076.05–6077.95 cm⁻¹, 5 coefficients, p/T not retrieved, residual noise, `fine_grid_margin_cm = 0.04`) | run_merlin.jl, run_hit08_narrow.jl, run_merlin_HT.jl, data_driven_noise.jl | hit20_H2O_narrow, merlin_12CH4, merlin_12CH4_HT, hit08_12CH4_narrow | merlin_window.toml; merlin_12CH4, merlin_12CH4_HT: ../make_merlin_lookup_tables.jl (vSmartMOM line-mixing branch) |

All tables are `.jld2`. Only the `hit24_12CH4*_ch4win`, `tccon_H2O_ch4win`, `tccon_HDO_ch4win`
and `hit20_12CO2_ch4win` tables exist now, so only `ch4_hit24`, `ch4_hit24_lm`,
`ch4_hit24_sdvlm` and their `_fixpT` variants pass `--dry-run` today; every other config fails
only on the missing table files listed above (checked by building each with an existing table
in their place).

1. Tables are built by the SpectralFits `lookup_table` app from the configs in
   `../spectra/tables/` (one variant per table; `julia -t auto --project=../../SpectralFit.jl
   ../../SpectralFit.jl/bin/run_config.jl ../spectra/tables/full_range.toml [table ...]`). Each
   build writes `<table>.toml` next to the table with its configuration, git hash and the
   interpolation error between nodes. The full-range tables are about 340 MB each. Tables
   written by older Julia/JLD2 versions do not load under Julia 1.12; rebuild them.
2. Commit 11ece7c copied run_hit08.jl over run_hit16.jl, run_hit20.jl and run_tccon.jl, so
   those experiments silently used the HITRAN 2008 CH4 and HITRAN 2020 CO2 tables; the configs
   use the intended line lists.
3. Nothing builds OCO_H2O.jld2 (OCO_spectra.hdf holds only CO2). v0.1 also took the h2o_test CO2 from OCO_12CO2_weakband.jld2, which
   nothing builds now; h2o_test.toml reads the ABSCO dry-air slice instead, as oco.toml does.

The `_fixpT` variants fix p and T to the tower record: `[conditions] file =
"../../data/tower_pT_20160921.csv"`, written by [export_tower_pT.jl](export_tower_pT.jl) from
`../fignewtons/timeseries/data/obs.jld2` (the v0.1-era times there are one day late; the CSV has
true UTC times). Each 15-minute average takes the row nearest its first sample time, within 60 s.

The test variants (`ch4_broadening_test_*`, `oco_h2o_test_*`) retrieve the whole of 20160921,
where the v0.1 test scripts fitted only its first 15-minute average, and start p and T from the
campaign prior (850 hPa, 285 K) rather than from the mean measured p and T.

## Analysis scripts

Run from this directory with `julia --project=../../SpectralFit.jl <script>`. They share
[common.jl](common.jl) (paths, `OLD_T_OFFSET_S`, `load_run`, `pT`, `dry_column`, tower and
v0.1-record loaders), read the latest run of each experiment from `../data/runs`, and skip
experiments that have not been run with a printed note.

| script | what |
|---|---|
| [compare_hit24.jl](compare_hit24.jl) `[ch4] [co2] [variants]` | HITRAN 2024 vs the stored campaign results (HITRAN 2008 CH4, OCO ABSCO CO2) and Picarro; `variants`: one row per hit24 / hit24_co2 experiment. Per-record CSVs in `../data/runs/hit24_comparison/` |
| [rmse_hit24.jl](rmse_hit24.jl) | spectral RMSE per hit24 / hit24_co2 experiment, OCO ABSCO RMSE from its stored χ² |
| [rmse_lines_hit24.jl](rmse_lines_hit24.jl) | RMSE on absorption-line channels only, for experiment pairs (HITRAN 2008 / 2024, 2024 / 2024 LM, 2024 / 2024 SDVLM, OCO / HITRAN 2024 CO2) |
| [compare_merlin_vs_hit08.jl](compare_merlin_vs_hit08.jl) | `r6_merlin` vs `r6_hit08_narrow` CH4 columns and χ² |
| [compare_all.jl](compare_all.jl) | `r6_hit08_narrow`, `r6_merlin`, `r6_merlin_HT` three-way comparison |
| [calibrate_noise.jl](calibrate_noise.jl) | assumed vs empirical (residual RMS) noise for the R(6) and hit24 runs |
| [export_tower_pT.jl](export_tower_pT.jl) | writes `../data/tower_pT_20160921.csv` for the `_fixpT` variants |

Times: `old_t(r) = machine_time − OLD_T_OFFSET_S` (2082758400 s) matches the `t` of the
v0.1-era results in `../fignewtons/timeseries/data`; the true LabVIEW-to-Unix offset is
2082844800 s.
