#
# Project : Gardenia
# Source  : nac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/09/27
#

struct ImagDomainData{T<:Real}
    N_imag::Int64               #The number of points used in Nevanlinna algorithm
    freq  ::Array{Complex{T},1} #The values of Matsubara frequencies
    val   ::Array{Complex{T},1} #The values of negative of Green function
end

struct RealDomainData{T<:Real}
    N_real  ::Int64               #The number of mesh in real axis
    w_max   ::Float64             #The energy cutoff of real axis
    eta     ::Float64             #The paramer. The retarded Green function is evaluated at omega+i*eta
    sum_rule::Float64             #The value of sum of spectral function
    freq    ::Array{Complex{T},1} #The values of frequencies of retarded Green function
    val     ::Array{Complex{T},1} #The values of negative of retarded Green function
end

mutable struct NevanlinnaSolver{T<:Real}
    imags::ImagDomainData{T}          #imaginary domain data
    reals::RealDomainData{T}          #real domain data
    phis::Vector{Complex{T}}          #phis in schur algorithm
    abcd::Array{Complex{T},3}         #continued fractions
    H_max::Int64                      #upper cut off of H
    H_min::Int64                      #lower cut off of H
    H::Int64                          #current value of H
    ab_coeff::Vector{ComplexF64}      #current solution for H
    hardy_matrix::Array{Complex{T},2} #hardy_matrix for H
    iter_tol::Int64                   #upper bound of iteration
    lambda::Float64                   #regularization parameter for second derivative term
    ini_iter_tol::Int64               #upper bound of iteration for H_min
    verbose::Bool                       
end

function NevanlinnaSolver(
                  wn          ::Vector{Complex{T}},
                  gw          ::Vector{Complex{T}},
                  N_real      ::Int64,
                  w_max       ::Float64,
                  eta         ::Float64,
                  sum_rule    ::Float64,
                  H_max       ::Int64,
                  iter_tol    ::Int64,
                  lambda      ::Float64
                  ;
                  verbose     ::Bool=false,
                  pick_check  ::Bool=true,
                  optimization::Bool=true,
                  ini_iter_tol::Int64=500,
                  mesh        ::Symbol=:linear,
                  ham_option  ::Bool=false #option for using in Hamburger moment problem
                  )::NevanlinnaSolver{T} where {T<:Real}

    if N_real%2 == 1
        error("N_real must be even number!")
    end

    @assert length(wn) == length(gw)
    N_imag = length(wn) 
    @show N_imag

    if pick_check
        opt_N_imag =  calc_opt_N_imag(N_imag, wn, gw, verbose=verbose)
    else 
        opt_N_imag = N_imag
    end
    @show opt_N_imag

    imags = ImagDomainData(wn, gw, opt_N_imag)
    reals = RealDomainData(N_real, w_max, eta, sum_rule, T=T, mesh=mesh)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)

    H_min::Int64 = 1
    ab_coeff = zeros(ComplexF64, 2*H_min)
    hardy_matrix = calc_hardy_matrix(reals, H_min)

    sol = NevanlinnaSolver(imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, lambda, ini_iter_tol, verbose)

    if ham_option
        return sol
    end
    
    if optimization
        calc_H_min(sol)
    else
        evaluation!(sol)
    end

    return sol
end

function solve(S::NevanACSolver, rd::RawData)
    N_real    = 1000  #demension of array of output
    omega_max = 10.0  #energy cutoff of real axis
    eta       = 0.001 #broaden parameter 
    sum_rule  = 1.0   #sum rule
    H_max     = 12    #cutoff of Hardy basis
    lambda    = 1e-4  #regularization parameter
    iter_tol  = 1000  #upper bound of iteration

    T = BigFloat
    setprecision(128)

    input_smpl = zeros(Complex{BigFloat},52)
    input_gw = zeros(Complex{BigFloat},52)

    dlm = readdlm("gw.data")
    @show size(dlm)
    @. input_smpl = dlm[:,1] * im
    @. input_gw = dlm[:,2] + dlm[:,3] * im
    wo_sol = NevanlinnaSolver(input_smpl, input_gw, N_real, omega_max, eta, sum_rule, H_max, iter_tol, lambda, verbose=true)

    open("twopeak_wo_opt.dat","w") do f
        for i in 1:wo_sol.reals.N_real
            println(f, "$(Float64(real.(wo_sol.reals.freq[i])))",  "\t", "$(Float64(imag.(wo_sol.reals.val[i]/pi)))")
        end
    end
end
