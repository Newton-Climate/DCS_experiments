#=

Empirical noise calibration from post-fit residual.

The noise model used in the retrieval (ShotNoise, or the residual noise of configs/merlin.toml
for the MERLIN-window runs) may not describe the residual: on this dataset the h5 `snr` and
SNRavg fields are all NaN, and the flat-region estimate at 6255.1–6255.4 cm⁻¹ over-predicts noise (χ² < 1). This
script measures the *actual* noise from the retrieval residual and rescales χ² to be
physically meaningful.

For each record and each retrieval:
    σ_empirical = r.rmse, RMS of measurement − model over the full-weight (unmasked) channels
    σ_SF        = 1 / sqrt(median diag(Sₑ⁻¹))       (SpectralFits' assumption)
    ratio       = σ_empirical / σ_SF

With σ_empirical taken over the same channels as χ², rescaling χ² by (σ_SF / σ_empirical)²
gives n / (n − m) for every record (n channels, m state elements), so the ratio is the
calibration: ratio ≈ 1 means the residual is explained by the assumed noise, ratio < 1 that
the noise model over-estimates σ (χ² < 1), ratio > 1 that the model leaves structure.
Rescaling the noise does not change the retrieved x, only the reported χ².

Experiments that have not been run are skipped with a note.

Usage: julia --project=../../SpectralFit.jl calibrate_noise.jl
=#

include(joinpath(@__DIR__, "common.jl"))

experiments = [("Voigt HIT08", "r6_hit08_narrow"), ("Voigt MERLIN", "r6_merlin"), ("HT    MERLIN", "r6_merlin_HT"),
               ("Voigt HIT24", "ch4_hit24"), ("HIT24 LM", "ch4_hit24_lm"), ("HIT24 SDVLM", "ch4_hit24_sdvlm")]

function calibrate(results)
    σ_emp = [r.rmse for r in results]                                   # post-fit residual RMS
    σ_SF  = [1 / sqrt(median(diag(r.obs_precision))) for r in results]  # Sₑ⁻¹ is Diagonal
    return (; σ_SF, σ_emp, ratio = σ_emp ./ σ_SF, χ_orig = [r.chi2_obs_per_dof for r in results])
end

println()
@printf "%-14s | %-13s | %-13s | %-11s | %-13s | %-13s\n" "" "σ_SF med" "σ_emp med" "SNR_emp" "χ² orig med" "σ_emp/σ_SF med"
println(repeat("-", 90))
for (name, experiment) in experiments
    rs = load_run(experiment)
    rs === nothing && continue
    c = calibrate(rs)
    # SNR uses continuum intensity ≈ 1 for DCS transmission fits, so SNR ≈ 1/σ
    snr_emp = 1 ./ c.σ_emp
    @printf "%-14s | %-13.5f | %-13.5f | %-11.1f | %-13.4f | %-13.4f\n" name  median(c.σ_SF) median(c.σ_emp) median(snr_emp) median(c.χ_orig) median(c.ratio)
end
