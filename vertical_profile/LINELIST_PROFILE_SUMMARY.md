# Line-list errors in the synthetic vertical-profile retrievals

Run on 2026-10-19 (Julia 1.12.7). These are synthetic vertical-column DCS spectra on the US Standard 0–10 km atmosphere with 35 layers. The empirical SNR is 709 for CH4 and 935 for CO2. The prior σ is 25% with a 0.5 km correlation length. Every run uses the same noise realization (seed 1).

For the target gas, the truth spectrum and the retrieval use different line lists. H2O, and CO2 in the CH4 window, use HITRAN 2020 in both. The **control** row retrieves with the truth's own line list. "Δ vs control" is the cross-list retrieval minus the control. It isolates the line-list error from noise and prior smoothing.

Diagnostics are computed in prior-σ-scaled coordinates: K̃ = K D, S̃ₐ⁻¹ = D Sₐ⁻¹ D, Ã = (K̃ᵀSₑ⁻¹K̃ + S̃ₐ⁻¹)⁻¹K̃ᵀSₑ⁻¹K̃, with D = diag(σ_prior). The condition number of the normal matrix is about 1e53 in raw units and 1.9e10 in scaled units. The remaining 1e10 comes from the Gaussian prior correlation matrix.

| Species | Truth | Retrieved with |
|---|---|---|
| CH4 | HITRAN 2008 | HITRAN 2008 (control), HITRAN 2016, HITRAN 2020, TCCON 2020 |
| CO2 | OCO ABSCO 5.1, dry | ABSCO (control), HITRAN 2016, HITRAN 2020, TCCON 2020 |
| CO2 | HITRAN 2020 | HITRAN 2020 (control), HITRAN 2016, TCCON 2020 |

ABSCO caveat: ABSCO contains all CO2 isotopologues. The HITRAN and TCCON models here are ¹²CO₂ (iso 1) only. In this window ABSCO's band strength is 1.8% above HITRAN 2020 iso 1, and ¹³CO₂ alone accounts for 1.1%. So the ABSCO-truth runs mix the isotopologue mismatch with the broadening differences. That is why the HITRAN 2020-truth CO2 set was added as a clean broadening-only case.

## Answer

**1. In profile mode, yes: line-list errors mostly distort the profile shape and leave the column nearly correct.**
- CH4, HITRAN 2008 truth retrieved with HITRAN 2016:
  - The column moves by only +0.48%.
  - The profile is distorted by 58% RMS, which is 3× the ~20% posterior layer σ.
  - Over 0–2 km the layers average −81%, going negative at the surface.
  - Over 2–5 km they average +57%.
- HITRAN 2020 behaves similarly: the column moves +0.61% and the shape by 95% RMS.
- TCCON is closest to HITRAN 2008: the column moves −0.31% and the shape by 15% RMS, with χ²/dof 1.26.
- CO2, HITRAN 2020 truth retrieved with HITRAN 2016: the column moves only +0.007%, and the shape moves 13% RMS (−9% at 0–2 km, +7% at 2–5 km).
- The column averaging kernel stays at 0.999–1.002 in every band and every run. So the column is protected by the retrieval's structure, not by luck.

**2. No: isobar mode does not reduce the distortion. It moves the line-list error into the column.**
Without a vertical pressure gradient, the retrieval can no longer absorb a broadening mismatch by moving gas between high- and low-pressure layers.

| Case | Profile mode: Δcol / shape RMS | Isobar mode: Δcol / shape RMS |
|---|---|---|
| CH4, HITRAN 2016 | +0.48% / 58% | **+3.98% / 230%** |
| CH4, HITRAN 2020 | +0.61% / 95% | +2.06% / 190% |
| CH4, TCCON | −0.31% / 15% | −0.54% / 34% |
| CO2 (HITRAN 2020 truth), HITRAN 2016 | +0.007% / 13% | +0.75% / 12% |
| CO2 (HITRAN 2020 truth), TCCON | +0.43% / 3.9% | +0.88% / 18% |

Isotherm mode looks like profile mode: the column changes by the same amount and the shape is somewhat larger. This is consistent with pressure broadening carrying the vertical information and temperature being redundant with it.

**3. Every cross-list run shows a spectral misfit.** χ²/dof is 4.5–8.6 for CH4 with HITRAN 2016 or 2020, and 3.1–3.9 for CO2 with ABSCO truth. It is only 1.06–1.9 for TCCON or HITRAN 2016 against HITRAN 2020. So the large shape errors appear together with residuals the retrieval cannot fit.

For the paper: in a profile retrieval, a pressure-broadening error between line lists mainly shows up as a spurious vertical gradient, typically low near the surface and high at 2–5 km. The column is biased only at the ~0.5% level. Forcing a single effective pressure (isobar) removes that outlet, and the same error then lands in the column at 1–4%.

## Tables

The per-run numbers below come from `linelist_results/linelist_summary_all.csv`, with column definitions under the table.

| Species | Truth | Retrieval | Mode | Column err % | Δcol vs control % | Shape RMS vs control % | Max shape vs control % | Col σ % | DOF total | DOF 0–2 | DOF 2–5 | DOF 5–10 | χ²/dof |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CH4 | HIT08 | HIT08 | profile | +0.027 | +0.000 | 0.000 | 0.000 | 0.041 | 3.32 | 0.99 | 0.73 | 1.60 | 1.029 |
| CH4 | HIT08 | HIT16 | profile | +0.506 | +0.479 | 57.793 | 135.542 | 0.043 | 3.33 | 0.99 | 0.73 | 1.61 | 6.584 |
| CH4 | HIT08 | HIT20 | profile | +0.634 | +0.606 | 94.656 | 191.769 | 0.042 | 3.34 | 0.99 | 0.73 | 1.61 | 4.621 |
| CH4 | HIT08 | TCCON | profile | -0.286 | -0.313 | 15.475 | 27.304 | 0.041 | 3.33 | 0.99 | 0.73 | 1.61 | 1.261 |
| CH4 | HIT08 | HIT08 | isotherm | +0.036 | +0.000 | 0.000 | 0.000 | 0.039 | 3.48 | 1.04 | 0.77 | 1.67 | 1.029 |
| CH4 | HIT08 | HIT16 | isotherm | +0.504 | +0.467 | 82.628 | 193.914 | 0.041 | 3.46 | 1.03 | 0.77 | 1.66 | 6.222 |
| CH4 | HIT08 | HIT20 | isotherm | +0.535 | +0.499 | 119.474 | 229.502 | 0.040 | 3.46 | 1.03 | 0.77 | 1.66 | 4.474 |
| CH4 | HIT08 | TCCON | isotherm | -0.140 | -0.177 | 15.687 | 31.086 | 0.039 | 3.48 | 1.04 | 0.77 | 1.67 | 1.254 |
| CH4 | HIT08 | HIT08 | isobar | +0.011 | +0.000 | 0.000 | 0.000 | 0.030 | 2.24 | 0.76 | 0.47 | 1.02 | 1.029 |
| CH4 | HIT08 | HIT16 | isobar | +3.990 | +3.979 | 230.477 | 428.541 | 0.031 | 2.19 | 0.74 | 0.45 | 0.99 | 8.565 |
| CH4 | HIT08 | HIT20 | isobar | +2.073 | +2.062 | 190.002 | 331.731 | 0.031 | 2.19 | 0.74 | 0.45 | 0.99 | 5.124 |
| CH4 | HIT08 | TCCON | isobar | -0.529 | -0.540 | 34.336 | 54.203 | 0.030 | 2.22 | 0.75 | 0.46 | 1.01 | 1.242 |
| CO2 | ABSCO51 | ABSCO51 | profile | +0.029 | +0.000 | 0.000 | 0.000 | 0.025 | 3.94 | 1.11 | 0.86 | 1.98 | 1.017 |
| CO2 | ABSCO51 | HIT16 | profile | +0.517 | +0.488 | 17.733 | 47.972 | 0.025 | 3.96 | 1.11 | 0.86 | 1.99 | 3.352 |
| CO2 | ABSCO51 | HIT20 | profile | +0.504 | +0.475 | 14.327 | 33.733 | 0.025 | 3.97 | 1.11 | 0.87 | 1.99 | 3.142 |
| CO2 | ABSCO51 | TCCON | profile | +0.941 | +0.912 | 14.261 | 27.409 | 0.025 | 3.95 | 1.11 | 0.86 | 1.98 | 3.268 |
| CO2 | ABSCO51 | ABSCO51 | isotherm | +0.022 | +0.000 | 0.000 | 0.000 | 0.023 | 4.12 | 1.16 | 0.91 | 2.06 | 1.017 |
| CO2 | ABSCO51 | HIT16 | isotherm | +0.486 | +0.464 | 17.971 | 48.961 | 0.023 | 4.14 | 1.16 | 0.91 | 2.07 | 3.309 |
| CO2 | ABSCO51 | HIT20 | isotherm | +0.474 | +0.452 | 17.486 | 35.961 | 0.023 | 4.14 | 1.16 | 0.91 | 2.07 | 3.096 |
| CO2 | ABSCO51 | TCCON | isotherm | +0.927 | +0.905 | 16.151 | 41.809 | 0.023 | 4.14 | 1.16 | 0.91 | 2.07 | 3.164 |
| CO2 | ABSCO51 | ABSCO51 | isobar | -0.008 | +0.000 | 0.000 | 0.000 | 0.015 | 2.80 | 0.91 | 0.61 | 1.28 | 1.017 |
| CO2 | ABSCO51 | HIT16 | isobar | +1.008 | +1.016 | 12.958 | 21.030 | 0.015 | 2.79 | 0.90 | 0.61 | 1.28 | 3.850 |
| CO2 | ABSCO51 | HIT20 | isobar | +0.251 | +0.260 | 19.338 | 36.745 | 0.015 | 2.79 | 0.90 | 0.61 | 1.28 | 3.070 |
| CO2 | ABSCO51 | TCCON | isobar | +1.130 | +1.138 | 12.197 | 20.006 | 0.015 | 2.79 | 0.90 | 0.61 | 1.28 | 3.315 |
| CO2 | HIT20 | HIT20 | profile | +0.029 | +0.000 | 0.000 | 0.000 | 0.025 | 3.97 | 1.11 | 0.87 | 1.99 | 1.017 |
| CO2 | HIT20 | HIT16 | profile | +0.035 | +0.007 | 13.058 | 35.423 | 0.025 | 3.96 | 1.11 | 0.86 | 1.99 | 1.422 |
| CO2 | HIT20 | TCCON | profile | +0.460 | +0.431 | 3.755 | 7.430 | 0.025 | 3.95 | 1.11 | 0.86 | 1.98 | 1.087 |
| CO2 | HIT20 | HIT20 | isotherm | +0.021 | +0.000 | 0.000 | 0.000 | 0.023 | 4.14 | 1.16 | 0.91 | 2.07 | 1.017 |
| CO2 | HIT20 | HIT16 | isotherm | +0.028 | +0.008 | 8.263 | 18.344 | 0.023 | 4.14 | 1.16 | 0.91 | 2.07 | 1.410 |
| CO2 | HIT20 | TCCON | isotherm | +0.469 | +0.448 | 5.683 | 13.842 | 0.023 | 4.14 | 1.16 | 0.91 | 2.07 | 1.062 |
| CO2 | HIT20 | HIT20 | isobar | -0.010 | +0.000 | 0.000 | 0.000 | 0.015 | 2.79 | 0.90 | 0.61 | 1.28 | 1.017 |
| CO2 | HIT20 | HIT16 | isobar | +0.740 | +0.750 | 11.969 | 28.636 | 0.015 | 2.79 | 0.90 | 0.61 | 1.28 | 1.873 |
| CO2 | HIT20 | TCCON | isobar | +0.865 | +0.875 | 18.193 | 32.791 | 0.015 | 2.79 | 0.90 | 0.61 | 1.28 | 1.173 |

Column definitions:
- **Shape RMS vs control**: the RMS over layers of (Δvmr_j/vmr_true_j − Δcol), in percent. It measures the profile change with the column change removed.
- **Col σ**: the posterior column 1σ.
- **DOF bands**: the sums of diag(Ã) over 0–2, 2–5 and 5–10 km.

Files in `linelist_results/`:
- `linelist_summary_{CH4,CO2,all}.csv`: one row per run. Besides the columns above, it also has the shape RMS against truth, the band-mean column AK, the condition numbers and the run time.
- `linelist_profiles_{CH4,CO2}.csv`: one row per run and layer. Columns: z, p, T, true and retrieved VMR, error, Δ vs control, posterior σ, layer DOF and column AK.
- `<species>_<mode>_truth-<A>_retr-<B>.jld2`: the full InversionResults for each run.

## How to run

The numbers above come from the v0.2 hand-port (`run_profiles.jl crosslist`, since removed;
outputs kept in `linelist_results/`). The experiments are now ProfileOSSE configs; see
README.md. From `vertical_profile/`:

```bash
julia --project=../../SpectralFit.jl make_absco_table.jl
```

```bash
julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/crosslist_ch4.toml
```

```bash
julia --project=../../SpectralFit.jl ../../SpectralFit.jl/bin/run_config.jl configs/crosslist_co2.toml
```

```bash
julia --project=../../SpectralFit.jl linelist_errors.jl
```

```bash
julia linelist_tables_md.jl > runs/linelist_analysis/summary_tables.md
```

`make_absco_table.jl` writes `../spectra/linelist_tables/ABSCO51_CO2_dry.jld2`. Single
variants (`truth-<A>_retr-<B>_<mode>`) can be named after the config file. Each retrieval
takes about 0.5–2 min. `linelist_errors.jl` writes the same CSV columns to
`runs/linelist_analysis/`, using the package's stored averaging kernel and posterior
covariance instead of the scaled diagnostics described above (same quantities).
