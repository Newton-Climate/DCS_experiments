using SpectralFits, LinearAlgebra
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
