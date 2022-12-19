using SpectralFits, LinearAlgebra, ForwardDiff, Plots, LaTeXStrings


height_coords(i) = 1e3*i*(0.4+0.02*i)

function temperature_profile(T₀, z; tropopause_height=11000)
    # make the standard atmospheric profile
    # info found here:

    
    trop_ind = findlast(z .< tropopause_height)
    #isotherm_ind = findlast(z .< 20000)
    T = zeros(length(z))
    
    T[1:trop_ind] =   T₀ .- 6.5e-3 .* z[1:trop_ind]
    T[trop_ind+1:end] .= T[trop_ind]
#    T[isotherm_ind+1:end] =   T[isotherm_ind] .+ 1.0e-3 .* (z[isotherm_ind+1:end] .- z[isotherm_ind])
    return T
end

    

function DOF(in::AbstractResults, ind)

    # calculate gain matrix
    G = calc_gain_matrix(in)

    # Averaging Kernal
    A = G * in.K

    degrees = tr(A[ind, ind])
    return degrees
end


function averaging_kernel(in::AbstractResults, vcd, ind)

    # calculate gain matrix
    G = calc_gain_matrix(in)

    H = vcd ./ sum(vcd)
    
    A = G * in.K
    

    cK = (H'*A[ind, ind]) ./ H'
    return cK'
end


function DOF_at_prior(f, xₐ, σ², inversion_results, ind)
        y = f(xₐ)
    if typeof(xₐ) <: AbstractDict
        xₐ = assemble_state_vector!(xₐ)
    end
    
    K = ForwardDiff.jacobian(f, xₐ)
    Sₑ⁻¹ = SpectralFits.make_obs_error(y, σ²=σ², linear=true)
                          

    G = inv(K'*inversion_results.Sₑ⁻¹*K + inversion_results.Sₐ⁻¹)*K'*inversion_results.Sₑ⁻¹
    A = G * K
    return tr(A[ind, ind])
end

function error_stats(results, vcd, ind, ϵ, xₐ, x_true; conv=1e6, outfile="error_stats.pdf")

    xₐ = assemble_state_vector!(xₐ)
    x_retrieved = assemble_state_vector!(result.x)
    x_true = assemble_state_vector!(x_true)

    
    H = vcd ./ sum(vcd)
    G = calc_gain_matrix(result)
    retrieval_noise = conv*G*ϵ
    retrieval_noise = retrieval_noise[ind] ./ vcd 
    A = G*result.K
    smoothing_error = conv*(A-I) * (x_retrieved-xₐ)
    smoothing_error = smoothing_error[ind] ./ vcd 
    retrieval_error = conv*(x_true - x_retrieved)[ind] ./ vcd

    posterior_uncertainty = sqrt.([results.S[i,i] for i in ind])
    posterior_uncertainty = conv*posterior_uncertainty ./ vcd
    

    out = Dict("smoothing_error" => smoothing_error,
               "retrieval_noise" => retrieval_noise,
               "A" => A,
               "retrieval_error" => retrieval_error,
               "posterior_uncertainty" => posterior_uncertainty)

    # plot(A[ind,:])
    return out
end


function plot_lhr_stats(stats_in, z)
    width = 2
    
    p = plot(stats_in["retrieval_error"], z, label="retrieval_error", lw=width, color=:red)
    plot!(stats_in["smoothing_error"], z, label="smoothing_error", lw=width, color=:blue)
    plot!(stats_in["retrieval_noise"], z, label="error from instrument noise", lw=width, color=:green)
    plot!(ylabel="height (m)", xlabel="ppm")

    return p
end

function plot_profile(result, z, vcd, xₐ, x_true, species, conv, xlabel)
    width = 2
    p = plot(conv*result.x[species] ./ vcd, z, label="retrieved", lw=width, color=:red)
    plot!(conv*x_true[species] ./ vcd, z, label="truth", lw=width, color=:black)
    plot!(conv*xₐ[species] ./ vcd, z, label="a priori", lw=width, color=:yellow)
    plot!(xlabel=xlabel, ylabel="height (m)")
    return p
end


function plot_fit(_result; outfile="fit.pdf")
    p1 = plot(result.grid, result.measurement, label="Observed", lw=2, color=:black)
    plot!(result.grid, result.model, label="modelled", lw=2, color=:orange)
    plot!(ylabel="intensity", title="fitted spectrum")

    p2 = plot(result.grid, result.model - result.measurement, label=:false, lw=2, color=:red)
    plot!(xlabel=L"wave-numbers (cm^{-1}", ylabel="model - measurement")

    p3 = plot(p1, p2, layout=(2,1))
    plot!(fontfamily="serif-roman", legendfont=font("Computer Modern", 7))
    savefig(outfile)
    return p3
end


"""Make custom prior error covarience matrix"""
function make_prior_error(σ::Union{Array{<:Real,1}, OrderedDict}, nGas, pr, pcorr)
    FT = Float64 # Should be dynamic in the future!
    # First get Diagonal elements (also to get the total size)
    out = FT[]
    nKey = zeros(Int,length(σ))
    i = 1
    for key in keys(σ)
        nKey[i] = length(σ[key]);
        i += 1
        out = append!(out, σ[key])
    end
    n = length(out)

    # Create prior covariance matrix:
    Sₐ = zeros(FT,n,n)
    # First fill the diagonal (will be done later as well but not for nonGases)
    for i=1:n
        Sₐ[i,i] = out[i]^2
    end

    # now start populating the Sa matrix for Gases (with off-diagonal elements)
    for i = 1:nGas
        offset = sum(nKey[1:i]) - nKey[i] 
        for i_1 = 1:nKey[i], i_2 = 1:nKey[i] 
            i_11 = i_1 + offset
            i_22 = i_2 + offset
            
            # Add correlation length scale here:
            Sₐ[i_11,i_22] = out[i_11] * out[i_22] *
                        exp(-(pr[i_1]-pr[i_2])^2/(2pcorr[i]^2))
        end
    end
    return inv(Sₐ)
end
