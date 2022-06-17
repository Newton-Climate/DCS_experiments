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
