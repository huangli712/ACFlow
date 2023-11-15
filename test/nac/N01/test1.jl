using SparseIR

T = Float64 #BigFloat
#setprecision(128)

#function generate_input_data(rho::Function, beta::Float64)
function generate_input_data(beta::Float64)
    lambda = 1e+4
    wmax = lambda / beta
    basis = FiniteTempBasisSet(beta, wmax, 1e-15)

    hnw = length(basis.smpl_wn_f.sampling_points)÷2
    
    input_smpl = Array{Complex{T}}(undef, hnw) 
    
    for i in 1:hnw
        @show i
        input_smpl[i]= SparseIR.valueim(basis.smpl_wn_f.sampling_points[hnw+i], beta)
    end
    
    return input_smpl#, input_gw
end

beta = 100.
input_smpl = generate_input_data(beta)
@show input_smpl
