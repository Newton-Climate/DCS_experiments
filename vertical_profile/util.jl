using SpectralFits, LinearAlgebra, ForwardDiff
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

    
