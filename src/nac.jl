#
# Project : Gardenia
# Source  : nac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/09/29
#

struct ImagDomainData
    Nopt :: I64         # The number of points used in Nevanlinna algorithm
    freq :: Vector{APC} # The values of Matsubara frequencies
    val  :: Vector{APC} # The values of negative of Green function
end

function ImagDomainData(wn::Vector{APC}, gw::Vector{APC}, Nopt::I64)
    freq = calc_mobius(wn[1:Nopt])
    val = calc_mobius(-gw[1:Nopt])

    success = calc_pick(Nopt, val, freq)
    if success
        println("Pick matrix is positive semi-definite.")
    else
        println("Pick matrix is non positive semi-definite matrix in Schur method.")
    end
    
    freq = reverse(wn[1:Nopt])
    val  = reverse(val)
    
    return ImagDomainData(Nopt, freq, val)
end

struct RealDomainData
    N_real  ::I64         # The number of mesh in real axis
    w_max   ::F64         # The energy cutoff of real axis
    eta     ::F64         # The paramer. The retarded Green function is evaluated at omega+i*eta
    sum_rule::F64         # The value of sum of spectral function
    freq    ::Vector{APC} # The values of frequencies of retarded Green function
    val     ::Vector{APC} # The values of negative of retarded Green function
end

function RealDomainData(N_real::I64, w_max::F64, eta::F64, sum_rule::F64)
    val = Array{APC}(collect(LinRange(-w_max, w_max, N_real)))
    freq = val .+ eta * im
    return RealDomainData(N_real, w_max, eta, sum_rule, freq, val)
end

mutable struct NevanlinnaSolver
    imags::ImagDomainData      # imaginary domain data
    reals::RealDomainData      # real domain data
    phis::Vector{APC}          # phis in schur algorithm
    abcd::Array{APC,3}         # continued fractions
    H_max::I64                 # upper cut off of H
    H_min::I64                 # lower cut off of H
    H::I64                     # current value of H
    ab_coeff::Vector{C64}      # current solution for H
    hardy_matrix::Array{APC,2} # hardy_matrix for H
    iter_tol::I64              # upper bound of iteration
    lambda::F64                # regularization parameter for second derivative term
    ini_iter_tol::I64          # upper bound of iteration for H_min
end

function NevanlinnaSolver(
                  wn          ::Vector{APC},
                  gw          ::Vector{APC},
                  N_real      ::I64,
                  w_max       ::F64,
                  eta         ::F64,
                  sum_rule    ::F64,
                  H_max       ::I64,
                  iter_tol    ::I64,
                  lambda      ::F64
                  ;
                  pick_check  ::Bool=true,
                  optimization::Bool=true,
                  ini_iter_tol::I64=500,
                  ham_option  ::Bool=false #option for using in Hamburger moment problem
                  )

    if N_real%2 == 1
        error("N_real must be even number!")
    end

    @assert length(wn) == length(gw)
    N_imag = length(wn) 
    @show N_imag

    if pick_check
        @time opt_N_imag =  calc_Nopt(wn, gw)
    else 
        opt_N_imag = N_imag
    end
    @show opt_N_imag

    imags = ImagDomainData(wn, gw, opt_N_imag)
    reals = RealDomainData(N_real, w_max, eta, sum_rule)

    phis = calc_phis(imags)
    abcd = calc_abcd(imags, reals, phis)

    H_min::Int64 = 1
    ab_coeff = zeros(ComplexF64, 2*H_min)
    hardy_matrix = calc_hardy_matrix(reals, H_min)

    sol = NevanlinnaSolver(imags, reals, phis, abcd, H_max, H_min, H_min, ab_coeff, hardy_matrix, iter_tol, lambda, ini_iter_tol)

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

function calc_mobius(z::Vector{APC})
    _z = similar(z)
    @. _z = (z - im) / (z + im)
    return _z
end

function calc_pick(k::I64, λ::Vector{APC}, ℎ::Vector{APC})
    pick = zeros(APC, k, k)
    for j = 1:k
        for i = 1:k
            num = one(APC) - λ[i] * conj(λ[j])
            den = one(APC) - ℎ[i] * conj(ℎ[j])
            pick[i,j] = num / den
        end
        pick[j,j] += APC(1e-250)
    end
    return issuccess(cholesky(pick, check = false))
end

function calc_Nopt(wn::Vector{APC}, gw::Vector{APC})
    N = length(wn)

    freq = calc_mobius(wn)
    val = calc_mobius(-gw)

    k = 0
    success = true
    while success && k ≤ N
        k += 1
        success = calc_pick(k, val, freq)
    end

    if !success
        println("N_imag is setted as $(k-1)")
        return k-1
    else
        println("N_imag is setted as $(N)")
        return N
    end
end

function calc_phis(imags::ImagDomainData)
    phis  = Array{APC}(undef, imags.Nopt) 
    abcds = Array{APC}(undef, 2, 2, imags.Nopt) 
    phis[1] = imags.val[1]
    
    for i in 1:imags.Nopt
        view(abcds,:,:,i) .= Matrix{APC}(I, 2, 2)
    end
    
    for j in 1:imags.Nopt-1
        for k in j+1:imags.Nopt
            prod = Array{APC}(undef, 2, 2) 
            prod[1,1] = (imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j]) * (imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[2,2] = one(APC)
            view(abcds,:,:,k) .= view(abcds,:,:,k)*prod
        end
        phis[j+1] = (-abcds[2,2,j+1]*imags.val[j+1] + abcds[1,2,j+1]) / (abcds[2,1,j+1]*imags.val[j+1] - abcds[1,1,j+1])
    end
    
    return phis
end

function calc_abcd(imags::ImagDomainData, reals::RealDomainData, phis::Vector{APC})
    abcd = Array{APC}(undef, 2, 2, reals.N_real) 

    for i in 1:reals.N_real
        result = Matrix{APC}(I, 2, 2) 
        z::APC = reals.freq[i]
        for j in 1:imags.Nopt
            prod = Array{APC}(undef, 2, 2)
            prod[1,1] = (z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[2,2] = one(APC)
            result *= prod
        end

        abcd[:,:,i] .= result
    end
    return abcd
end

function check_causality(hardy_matrix::Array{Complex{T},2},
                         ab_coeff::Vector{Complex{S}})::Bool where {S<:Real, T<:Real}

    param = hardy_matrix*ab_coeff

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        println("max_theta=",max_theta)
        println("hardy optimization was success.")
        causality = true
    else
        println("max_theta=",max_theta)
        println("hardy optimization was failure.")
        causality = false
    end
    return causality
end

function evaluation!(sol::NevanlinnaSolver)::Bool

    causality = check_causality(sol.hardy_matrix, sol.ab_coeff)
    if causality
        param = sol.hardy_matrix*sol.ab_coeff
        theta = (sol.abcd[1,1,:].* param .+ sol.abcd[1,2,:]) ./ (sol.abcd[2,1,:].*param .+ sol.abcd[2,2,:])
        sol.reals.val .= im * (one(APC) .+ theta) ./ (one(APC) .- theta)
    end

    return causality
end

function hardy_basis(z::APC, k::I64)
    w = (z-im)/(z+im)
    0.5*im*(w^(k+1)-w^k)/(sqrt(pi))
end

function calc_hardy_matrix(reals::RealDomainData, H::I64)
    hardy_matrix = zeros(APC, reals.N_real, 2*H)
    for k in 1:H
        hardy_matrix[:,2*k-1] .=      hardy_basis.(reals.freq,k-1)
        hardy_matrix[:,2*k]   .= conj(hardy_basis.(reals.freq,k-1))
    end
    return hardy_matrix
end

function calc_H_min(sol::NevanlinnaSolver,)::Nothing
    H_bound::Int64 = 50
    for iH in 1:H_bound
        println("H=$(iH)")
        zero_ab_coeff = zeros(ComplexF64, 2*iH)

        causality, optim = hardy_optim!(sol, iH, zero_ab_coeff, iter_tol=sol.ini_iter_tol)

        #break if we find optimal H in which causality is preserved and optimize is successful
        if causality && optim
            sol.H_min = sol.H
            break
        end

        if isdefined(Main, :IJulia)
            Main.IJulia.stdio_bytes[] = 0
        end

        if iH == H_bound
            error("H_min does not exist")
        end
    end
end

function calc_functional(
                    sol::NevanlinnaSolver,
                    H::Int64, 
                    ab_coeff::Vector{Complex{S}}, 
                    hardy_matrix::Array{Complex{T},2};
                    )::Float64 where {S<:Real, T<:Real}

    param = hardy_matrix*ab_coeff

    theta = (sol.abcd[1,1,:].* param .+ sol.abcd[1,2,:]) ./ (sol.abcd[2,1,:].*param .+ sol.abcd[2,2,:])
    green = im * (one(T) .+ theta) ./ (one(T) .- theta)
    A = Float64.(imag(green)./pi)

    tot_int = trapz(sol.reals.freq, A)
    second_der = integrate_squared_second_deriv(sol.reals.freq, A) 

    max_theta = findmax(abs.(param))[1]
    func = abs(sol.reals.sum_rule-tot_int)^2 + sol.lambda*second_der

    return func
end

function hardy_optim!(
                sol::NevanlinnaSolver,
                H::Int64,
                ab_coeff::Array{ComplexF64,1};
                iter_tol::Int64=sol.iter_tol,
                )::Tuple{Bool, Bool}

    loc_hardy_matrix = calc_hardy_matrix(sol.reals, H)

    function functional(x::Vector{ComplexF64})::Float64
        return calc_functional(sol, H, x, loc_hardy_matrix)
    end

    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = true))
    
    if  !(Optim.converged(res))
        println("Faild to optimize!")
    end
    
    causality = check_causality(loc_hardy_matrix, Optim.minimizer(res))

    if causality && (Optim.converged(res))
        sol.H = H
        sol.ab_coeff = Optim.minimizer(res)
        sol.hardy_matrix = loc_hardy_matrix
        evaluation!(sol)
    end
    
    return causality, (Optim.converged(res))
end

"""
Compute second derivative
If the length of `x` and `y` is `N`, the length of the returned vector is `N-2`.
"""
function second_deriv(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(ArgumentError("x and y must be the same length")) 
    end

    N = length(x)
    dx_backward = view(x, 2:(N-1)) - view(x, 1:(N-2))
    dx_forward = view(x, 3:N) - view(x, 2:(N-1))

    y_forward = view(y, 3:N)
    y_mid = view(y, 2:(N-1))
    y_backward = view(y, 1:(N-2))

    n = dx_backward .* y_forward + dx_forward .* y_backward - (dx_forward + dx_backward) .* y_mid
    d = (dx_forward.^2) .* dx_backward + (dx_backward.^2) .* dx_forward
    return 2 .* n ./ d
end

"""
Integrate the squarre of the abs of the second derivative
"""
function integrate_squared_second_deriv(x::AbstractVector, y::AbstractVector)
    N = length(x)
    sd = second_deriv(x, y)

    x_sd = view(x, 2:(N-1))
    return trapz(x_sd, abs.(sd) .^ 2)
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
    wo_sol = NevanlinnaSolver(input_smpl, input_gw, N_real, omega_max, eta, sum_rule, H_max, iter_tol, lambda)

    open("twopeak_wo_opt.dat","w") do f
        for i in 1:wo_sol.reals.N_real
            println(f, "$(Float64(real.(wo_sol.reals.freq[i])))",  "\t", "$(Float64(imag.(wo_sol.reals.val[i]/pi)))")
        end
    end
end
