#=

CO2 cross-section table from the OCO ABSCO v5.1 file (../spectra/OCO_spectra.hdf),
used as the "truth" line list for the CO2 cross-line-list profile experiments
(configs/crosslist_co2.toml). Standalone: needs nothing from the other scripts here.

ABSCO is kept on its native nodes (0.01 cm⁻¹, 25 hPa, 10 K), cropped to the
6175–6265 cm⁻¹ window and to the p/T box of the US-Standard 0–10 km profile, with
the H2O broadener VMR set to 0 (dry), so the comparison with the HITRAN/TCCON
retrievals isolates the CO2 line list. A cubic B-spline over (ν, p, T) is stored as a
vSmartMOM InterpolationModel, loadable with load_cross_section_model(LookupTable, path).

Output: ../spectra/linelist_tables/ABSCO51_CO2_dry.jld2

Run:  julia --project=../../SpectralFit.jl make_absco_table.jl

=#

using SpectralFits, Interpolations, Printf
using vSmartMOM.Absorption: InterpolationModel, save_interpolation_model
import HDF5

const SPECTRA_DIR = joinpath(@__DIR__, "../spectra/")
const OUT_DIR     = mkpath(joinpath(SPECTRA_DIR, "linelist_tables"))
const OUT_FILE    = joinpath(OUT_DIR, "ABSCO51_CO2_dry.jld2")

const ν_RANGE = (6175.0, 6265.0)     # cm⁻¹, covers the 6180–6260 retrieval window
const P_RANGE = (250.0, 1026.0)      # hPa, covers 264–1013 hPa of the profile
const T_RANGE = (220.0, 290.0)       # K, covers 223–288 K of the profile
const BROADENER_INDEX = 1            # H2O VMR 0.00 (dry); 2 → 0.03, 3 → 0.06

"""Indices of `v` inside [lo, hi] as a unit range."""
inside(v, lo, hi) = findfirst(>=(lo - 1e-6), v):findlast(<=(hi + 1e-6), v)

function make_absco_table()
    HDF5.h5open(joinpath(SPECTRA_DIR, "OCO_spectra.hdf"), "r") do fid
        ν = HDF5.read(fid, "Wavenumber")
        p = HDF5.read(fid, "Pressure") .* 0.01          # Pa → hPa
        T = HDF5.read(fid, "Temperature")               # Julia order (T, p)
        broad = HDF5.read(fid, "Broadener_01_VMR")

        # the T grid is the same at every pressure in this file; check rather than assume
        Tgrid = size(T, 1) == length(p) ? T[1, :] : T[:, 1]
        @assert all(t -> t ≈ Tgrid, eachrow(size(T, 1) == length(p) ? T : permutedims(T))) "ABSCO T grid varies with p"

        iν = inside(ν, ν_RANGE...)
        ip = inside(p, P_RANGE...)
        iT = inside(Tgrid, T_RANGE...)
        @printf "ABSCO %s: ν %.2f–%.2f (%d), p %.2f–%.2f hPa (%d), T %.0f–%.0f K (%d), broadener H2O VMR %.2f\n" HDF5.read_attribute(fid, "version") ν[iν[1]] ν[iν[end]] length(iν) p[ip[1]] p[ip[end]] length(ip) Tgrid[iT[1]] Tgrid[iT[end]] length(iT) broad[BROADENER_INDEX]

        # Julia reads the C-ordered {p, T, broadener, ν} dataset as (ν, broadener, T, p)
        σ = fid["Gas_02_Absorption"][iν, BROADENER_INDEX, iT, ip]   # (ν, T, p)
        σ = permutedims(σ, (1, 3, 2))                                 # (ν, p, T)
        @assert all(isfinite, σ) && all(>=(0), σ)

        # uniform grids as ranges (InterpolationModel needs AbstractRange)
        νr = range(ν[iν[1]], ν[iν[end]], length=length(iν))
        pr = range(p[ip[1]], p[ip[end]], length=length(ip))
        Tr = range(Tgrid[iT[1]], Tgrid[iT[end]], length=length(iT))
        @assert maximum(abs, νr .- ν[iν]) < 1e-6 && maximum(abs, pr .- p[ip]) < 1e-6 && maximum(abs, Tr .- Tgrid[iT]) < 1e-6

        itp = interpolate(Float64.(σ), BSpline(Cubic(Line(OnGrid()))))
        model = InterpolationModel(itp, 2, 1, νr, pr, Tr)
        save_interpolation_model(model, OUT_FILE)
        @printf "wrote %s (%.1f MB)\n" OUT_FILE filesize(OUT_FILE) / 1e6
        return model
    end
end

"""Integrated band strength of ABSCO vs HITRAN 2020 12CO2 (iso 1) and HITRAN 2016 13CO2."""
function compare_band_strength(model)
    ν = 6180.0:0.002:6260.0
    absco = LookupTable(model)
    hitran(file, iso) = load_cross_section_model(HitranVoigt, joinpath(SPECTRA_DIR, file); molecule = 2,
                                                 isotope = iso, range_cm = (first(ν), last(ν)))
    σ(m, p, T) = cross_section(m, ν, GasConditions(p, T, 0.0))
    hit20 = hitran("hit20_12CO2.par", 1)
    for (pp, TT) in ((1000.0, 288.0), (500.0, 250.0))
        a = σ(absco, pp, TT)
        h = σ(hit20, pp, TT)
        @printf "p=%4.0f hPa T=%3.0f K: ∫σ ABSCO / HITRAN20 12CO2 = %.4f, max|Δσ|/max σ = %.3g\n" pp TT sum(a) / sum(h) maximum(abs, a .- h) / maximum(h)
    end
    for iso in 2:3
        try
            h13 = hitran("hit16_13CO2.par", iso)
            @printf "HITRAN16 CO2 iso %d / HITRAN20 iso 1 band strength in window: %.4f\n" iso sum(σ(h13, 1000.0, 288.0)) / sum(σ(hit20, 1000.0, 288.0))
        catch e
            println("HITRAN16 CO2 iso $iso: none in window ($(typeof(e)))")
        end
    end
end

model = make_absco_table()
compare_band_strength(model)
