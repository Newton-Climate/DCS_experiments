#=

Build the smallest cross-section interpolation tables that cover the US Standard synthetic
vertical-profile experiments (configs/us_standard_ch4.toml, configs/us_standard_co2.toml).

The tables only span the p/T box of the US-Standard 0–10 km profile used there
(profile, isotherm and isobar modes all fall inside it) and only the two
retrieval windows. The p and T resolution is chosen by brute force: line-by-line
cross-sections are computed on a fine master grid, every nested coarser grid is
tried, and the one with the fewest nodes is kept whose interpolation error in
the forward model stays below `TOL` × the data-driven 1σ noise, for every layer
of every mode. Coefficients are stored as Float32 to halve the file size.

Atmosphere, windows, measurement spaces, line lists and the noise template are read from the
line-by-line variants (profile, isotherm, isobar) of those configs, so the tables follow
the configs. (The tables in ../spectra/profile_tables were built on 2026-10-19 from the
run_profiles.jl atmosphere, whose US Standard constants differ from the package's in the
fourth digit; the p/T box is rounded outward and still covers it.)

Output: ../spectra/profile_tables/<window>_<molecule>.jld2 and a summary in
../spectra/profile_tables/README.txt (overwritten)

Run:  julia --project=../../SpectralFit.jl make_profile_tables.jl

=#

include(joinpath(@__DIR__, "util.jl"))
using Dates
using vSmartMOM.Absorption: InterpolationModel, save_interpolation_model
using Interpolations: interpolate, BSpline, Cubic, Line, OnGrid

const OUT_DIR = joinpath(SPECTRA_DIR, "profile_tables")
const TOL = 0.1   # allowed max table error per gas, in units of the 1σ measurement noise
const MODES = (:profile, :isotherm, :isobar)

# line-by-line configs per window and mode
const CONFIGS = Dict(w => Dict(m => profile_config("us_standard_$(lowercase(String(w)))", m) for m in MODES)
                     for w in (:CH4, :CO2))

"""(p, T) of every layer in every mode of window `w`."""
mode_layers(w) = Dict(m => (a = ProfileOSSE.profile_atmosphere(c); (a.pressure_hPa, a.temperature_K))
                      for (m, c) in CONFIGS[w])

# table box: the profile's p/T extremes, rounded outward so the end nodes are not
# hit exactly by floating-point round-off
const ATM = ProfileOSSE.profile_atmosphere(CONFIGS[:CH4][:profile])
const P_RANGE = (floor(minimum(ATM.pressure_hPa)), ceil(maximum(ATM.pressure_hPa)))
const T_RANGE = (floor(minimum(ATM.temperature_K)), ceil(maximum(ATM.temperature_K)))

# nested candidate sizes: every coarser grid's nodes are a subset of the master grid
const NP_MASTER, NT_MASTER = 33, 17
const NP_CANDIDATES = (3, 5, 9, 17, 33)
const NT_CANDIDATES = (3, 5, 9, 17)

"""Cubic B-spline table on the given p/T nodes, built like vSmartMOM.make_interpolation_model."""
function table_model(cs, mol, iso, ν, p_grid, T_grid; FT=Float64)
    itp = interpolate(FT.(cs), BSpline(Cubic(Line(OnGrid()))))
    return InterpolationModel(itp, mol, iso, ν, p_grid, T_grid)
end

table_bytes(nν, np, nT; bytes=4) = nν * (np + 2) * (nT + 2) * bytes   # prefiltered coefficients are padded by 1 on each side

"""Data-driven noise and instrument for a window (the configs' template record)."""
function window_noise(config, ν, space)
    meas = ProfileOSSE.template_measurement(config)
    σ_I  = sqrt(meas.noise_variance[1])
    # noise in the observable: log-intensity for LogSpace, intensity otherwise
    σ_obs = space isa LogSpace ? σ_I / mean(meas.y) : σ_I
    instrument = bind_fine_grid(prepare(DCS(), meas), ν)
    return (; instrument, σ_obs, baseline = maximum(meas.y), snr = mean(meas.y) / σ_I)
end

"""Forward-model observable for optical depth τ alone: log T (LogSpace) or baseline × T."""
function observable(τ, ν, noise, space)
    t = apply_instrument(noise.instrument, ν, exp.(-τ), nothing)
    return space isa LogSpace ? log.(t) : noise.baseline .* t
end

σ_at(model, ν, p, T) = cross_section(model, ν, GasConditions(p, T, 0.0))

function build_window(species::Symbol; io=stdout)
    wname = String(species)
    config = CONFIGS[species][:profile]
    lo, hi = config.window_cm
    ν, space = lo:config.fine_grid_spacing_cm:hi, config.measurement_space
    noise = window_noise(config, ν, space)
    layers = mode_layers(species)
    vcd = ProfileOSSE.profile_atmosphere(config).dry_air_column_cm2
    p_master = range(P_RANGE..., length=NP_MASTER)
    T_master = range(T_RANGE..., length=NT_MASTER)
    @printf io "\n== %s window %.0f–%.0f cm⁻¹: %d points, σ_obs = %.3g (SNR %.0f)\n" wname ν[1] ν[end] length(ν) noise.σ_obs noise.snr

    chosen = Dict{Symbol,Any}()
    for name in config.gases
        molnum = config.table["gases"][String(name)]["cross_section"]["molecule"]
        lbl_model = config.absorbers[name].model
        lbl(p, T) = σ_at(lbl_model, ν, p, T)

        # line-by-line master grid and reference spectra at every mode layer
        t = @elapsed begin
            master = zeros(length(ν), NP_MASTER, NT_MASTER)
            for (i, p) in enumerate(p_master), (j, T) in enumerate(T_master)
                master[:, i, j] .= lbl(p, T)
            end
        end
        col = config.truth_vmr[name] .* vcd
        refs = Dict(m => [lbl(p[k], T[k]) for k in eachindex(p)] for (m, (p, T)) in layers)
        obs_ref = Dict(m => observable(sum(col .* refs[m]), ν, noise, space) for m in keys(layers))
        @printf io "%s: master %d×%d line-by-line grid in %.0f s\n" name NP_MASTER NT_MASTER t

        # max table error over all modes, in units of σ_obs
        function table_error(model)
            table = LookupTable(model)
            err = 0.0
            for (m, (p, T)) in layers
                τ = sum(col[k] .* σ_at(table, ν, p[k], T[k]) for k in eachindex(p))
                err = max(err, maximum(abs, observable(τ, ν, noise, space) .- obs_ref[m]) / noise.σ_obs)
            end
            return err
        end

        # try every nested grid, smallest first
        candidates = sort([(np, nT) for np in NP_CANDIDATES, nT in NT_CANDIDATES][:], by=c -> (c[1] + 2) * (c[2] + 2))
        best = nothing
        for (np, nT) in candidates
            ip = 1:(NP_MASTER - 1) ÷ (np - 1):NP_MASTER
            iT = 1:(NT_MASTER - 1) ÷ (nT - 1):NT_MASTER
            model = try
                table_model(master[:, ip, iT], molnum, 1, ν, p_master[ip], T_master[iT]; FT=Float32)
            catch e
                @printf io "   %2d p × %2d T nodes: cannot build a cubic table (%s)\n" np nT typeof(e)
                continue
            end
            err = table_error(model)
            @printf io "   %2d p × %2d T nodes (Δp %5.1f hPa, ΔT %4.1f K): max error %.3g σ, %.1f MB\n" np nT step(p_master[ip]) step(T_master[iT]) err table_bytes(length(ν), np, nT) / 1e6
            if err <= TOL
                best = (; np, nT, err, model)
                break
            end
        end
        best === nothing && error("$name: no candidate grid reaches $(TOL)σ; extend NP/NT_MASTER")

        path = joinpath(OUT_DIR, "$(wname)_$(name).jld2")
        save_interpolation_model(best.model, path)
        @printf io "-> %s: %d p × %d T, max error %.3g σ, file %.1f MB\n" basename(path) best.np best.nT best.err filesize(path) / 1e6
        chosen[name] = (; best.np, best.nT, best.err, path)
    end

    # combined check: all of this window's tables together, loaded from disk as the *_tables variants will
    names = config.gases
    tables = Dict(name => load_cross_section_model(LookupTable, chosen[name].path) for name in names)
    lbls = Dict(name => config.absorbers[name].model for name in names)
    worst = 0.0
    for (m, (p, T)) in layers
        τ_tab = sum(config.truth_vmr[name][k] * vcd[k] .* σ_at(tables[name], ν, p[k], T[k]) for name in names for k in eachindex(p))
        τ_lbl = sum(config.truth_vmr[name][k] * vcd[k] .* σ_at(lbls[name], ν, p[k], T[k]) for name in names for k in eachindex(p))
        e = maximum(abs, observable(τ_tab, ν, noise, space) .- observable(τ_lbl, ν, noise, space)) / noise.σ_obs
        @printf io "   all gases, %-8s: max error %.3g σ\n" m e
        worst = max(worst, e)
    end
    return (; chosen, worst, noise)
end

"""Write to the terminal and to the README log at once."""
struct Tee <: IO; a::IO; b::IO; end
Base.unsafe_write(t::Tee, p::Ptr{UInt8}, n::UInt) = (unsafe_write(t.a, p, n); unsafe_write(t.b, p, n); flush(t.a); flush(t.b); n)

function main()
    open(joinpath(OUT_DIR, "README.txt"), "w") do file
        io = Tee(stdout, file)
        println(io, "Cross-section tables for the *_tables variants of vertical_profile/configs/us_standard_*.toml, made by vertical_profile/make_profile_tables.jl on $(Dates.today())")
        @printf io "p box %.0f–%.0f hPa, T box %.0f–%.0f K; tolerance %.2g σ per gas\n" P_RANGE... T_RANGE... TOL
        for species in (:CH4, :CO2)
            build_window(species; io=io)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    mkpath(OUT_DIR)
    main()
end
