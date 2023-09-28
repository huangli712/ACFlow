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

function ImagDomainData(wn     ::Array{Complex{T},1},
                        gw     ::Array{Complex{T},1},
                        N_imag ::Int64;
                        verbose::Bool = false
                        )::ImagDomainData{T} where {T<:Real}

    val  = Array{Complex{T}}(undef, N_imag) 
    freq = Array{Complex{T}}(undef, N_imag) 
    
    for i in 1:N_imag
        freq[i] = wn[i]
        val[i]  = (-gw[i] - im) / (-gw[i] + im) 
    end
    
    Pick = Array{Complex{T}}(undef, N_imag, N_imag)
    
    for j in 1:N_imag
        for i in 1:N_imag
            freq_i = (freq[i] - im) / (freq[i] + im)
            freq_j = (freq[j] - im) / (freq[j] + im)
            nom = one(T) - val[i] * conj(val[j])
            den = one(T) - freq_i * conj(freq_j)
            Pick[i,j] = nom / den
        end
        Pick[j,j] += T(1e-250)
    end
    
    success = issuccess(cholesky(Pick,check = false))
    
    if verbose
        if success
            println("Pick matrix is positive semi-definite.")
        else
            println("Pick matrix is non positive semi-definite matrix in Schur method.")
        end
    end
    
    freq = reverse(freq)
    val  = reverse(val)
    
    return ImagDomainData(N_imag, freq, val)
end

struct RealDomainData{T<:Real}
    N_real  ::Int64               #The number of mesh in real axis
    w_max   ::Float64             #The energy cutoff of real axis
    eta     ::Float64             #The paramer. The retarded Green function is evaluated at omega+i*eta
    sum_rule::Float64             #The value of sum of spectral function
    freq    ::Array{Complex{T},1} #The values of frequencies of retarded Green function
    val     ::Array{Complex{T},1} #The values of negative of retarded Green function
end

function RealDomainData(N_real  ::Int64,
                        w_max   ::Float64,
                        eta     ::Float64,
                        sum_rule::Float64
                        ;
                        T::Type=BigFloat,
                        small_omega::Float64 = 1e-5,
                        mesh::Symbol=:linear
                        )::RealDomainData{T}

    if mesh === :linear
        val = Array{Complex{T}}(collect(LinRange(-w_max, w_max, N_real)))
        freq = val .+ eta * im
        return RealDomainData(N_real, w_max, eta, sum_rule, freq, val)
    elseif mesh === :log
        half_N = N_real ÷ 2
        mesh = exp.(LinRange(log.(small_omega), log.(w_max), half_N))
        val = Array{Complex{T}}([reverse(-mesh); mesh])
        freq = val .+ eta * im
        return RealDomainData(N_real, w_max, eta, sum_rule, freq, val)
    elseif mesh === :test
        val  = Array{Complex{T}}(undef, N_real) 
        freq = Array{Complex{T}}(undef, N_real) 
        inter::T = big(2.0*w_max) / (N_real-1)
        temp ::T = big(-w_max)
        freq[1] = -big(w_max) + big(eta)*im
        for i in 2:N_real
            temp += inter
            freq[i] = temp + big(eta)*im
        end
        return RealDomainData(N_real, w_max, eta, sum_rule, freq, val)
    else
        throw(ArgumentError("Invalid mesh"))
    end
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
        @time opt_N_imag =  calc_Nopt(wn, gw)
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

function calc_Nopt(wn::Vector{APC}, gw::Vector{APC})
    N = length(wn)

    freq = (wn  .- im) ./ (wn  .+ im)
    val  = (-gw .- im) ./ (-gw .+ im)

    k = 0
    success = true

    Pick = Array{APC}(undef, N, N)
    while success
        k += 1

        for j in 1:k
            for i in 1:k
                num = one(APC) - val[i]  * conj(val[j])
                den = one(APC) - freq[i] * conj(freq[j])
                Pick[i,j] = num / den
            end

            Pick[j,j] += APC(1e-250)
        end

        success = issuccess(cholesky(Pick[1:k,1:k],check = false))

        if k == N
            break
        end
    end

    @show "haha"
    if !(success)
        println("N_imag is setted as $(k-1)")
    else
        println("N_imag is setted as $(N)")
    end

    if !(success)
        return (k-1)
    else
        return (N)
    end
end

function calc_phis(imags::ImagDomainData{T})::Vector{Complex{T}} where {T<:Real}
    phis  = Array{Complex{T}}(undef, imags.N_imag) 
    abcds = Array{Complex{T}}(undef, 2, 2, imags.N_imag) 
    phis[1] = imags.val[1]
    
    for i in 1:imags.N_imag
        view(abcds,:,:,i) .= Matrix{Complex{T}}(I, 2, 2)
    end
    
    for j in 1:imags.N_imag-1
        for k in j+1:imags.N_imag
            prod = Array{Complex{T}}(undef, 2, 2) 
            prod[1,1] = (imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j]) * (imags.freq[k] - imags.freq[j]) / (imags.freq[k] - conj(imags.freq[j]))
            prod[2,2] = one(T)
            view(abcds,:,:,k) .= view(abcds,:,:,k)*prod
        end
        phis[j+1] = (-abcds[2,2,j+1]*imags.val[j+1] + abcds[1,2,j+1]) / (abcds[2,1,j+1]*imags.val[j+1] - abcds[1,1,j+1])
    end
    
    return phis
end

function calc_abcd(imags::ImagDomainData{T}, 
                   reals::RealDomainData{T}, 
                   phis::Vector{Complex{T}}
                   )::Array{Complex{T},3} where {T<:Real}
    abcd = Array{Complex{T}}(undef, 2, 2, reals.N_real) 

    for i in 1:reals.N_real
        result = Matrix{Complex{T}}(I, 2, 2) 
        z::Complex{T} = reals.freq[i]
        for j in 1:imags.N_imag
            prod = Array{Complex{T}}(undef, 2, 2)
            prod[1,1] = (z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[1,2] = phis[j]
            prod[2,1] = conj(phis[j])*(z - imags.freq[j]) / (z - conj(imags.freq[j]))
            prod[2,2] = one(T)
            result *= prod
        end

        abcd[:,:,i] .= result
    end
    return abcd
end

function check_causality(hardy_matrix::Array{Complex{T},2},
                         ab_coeff::Vector{Complex{S}};
                         verbose::Bool=false
                         )::Bool where {S<:Real, T<:Real}

    param = hardy_matrix*ab_coeff

    max_theta = findmax(abs.(param))[1]
    if max_theta <= 1.0
        if verbose
           println("max_theta=",max_theta)
           println("hardy optimization was success.")
        end
        causality = true
    else
        if verbose
          println("max_theta=",max_theta)
          println("hardy optimization was failure.")
        end
        causality = false
    end
    return causality
end

function evaluation!(sol::NevanlinnaSolver{T};
                     verbose::Bool=false
                    )::Bool where {T<:Real}

    causality = check_causality(sol.hardy_matrix, sol.ab_coeff, verbose=verbose)
    if causality
        param = sol.hardy_matrix*sol.ab_coeff
        theta = (sol.abcd[1,1,:].* param .+ sol.abcd[1,2,:]) ./ (sol.abcd[2,1,:].*param .+ sol.abcd[2,2,:])
        sol.reals.val .= im * (one(T) .+ theta) ./ (one(T) .- theta)
    end

    return causality
end

function hardy_basis(z::Complex{T}, k::Int64) where {T<:Real}
    w = (z-im)/(z+im)
    0.5*im*(w^(k+1)-w^k)/(sqrt(pi))
end

function calc_hardy_matrix(reals::RealDomainData{T}, 
                           H::Int64
                           )::Array{Complex{T}, 2} where {T<:Real}
    hardy_matrix = Array{Complex{T}}(undef, reals.N_real, 2*H)
    for k in 1:H
        hardy_matrix[:,2*k-1] .=      hardy_basis.(reals.freq,k-1)
        hardy_matrix[:,2*k]   .= conj(hardy_basis.(reals.freq,k-1))
    end
    return hardy_matrix
end

function calc_H_min(sol::NevanlinnaSolver{T},)::Nothing where {T<:Real}
    H_bound::Int64 = 50
    for iH in 1:H_bound
        if sol.verbose
            println("H=$(iH)")
        end
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
                    sol::NevanlinnaSolver{T},
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
                sol::NevanlinnaSolver{T},
                H::Int64,
                ab_coeff::Array{ComplexF64,1};
                iter_tol::Int64=sol.iter_tol,
                )::Tuple{Bool, Bool} where {T<:Real}

    loc_hardy_matrix = calc_hardy_matrix(sol.reals, H)

    function functional(x::Vector{ComplexF64})::Float64
        return calc_functional(sol, H, x, loc_hardy_matrix)
    end

    function jacobian(J::Vector{ComplexF64}, x::Vector{ComplexF64})
        J .= gradient(functional, x)[1] 
    end

    res = optimize(functional, jacobian, ab_coeff, BFGS(), 
                   Optim.Options(iterations = iter_tol,
                                 show_trace = sol.verbose))
    
    if  !(Optim.converged(res)) && sol.verbose
        println("Faild to optimize!")
    end
    
    causality = check_causality(loc_hardy_matrix, Optim.minimizer(res), verbose=sol.verbose)

    if causality && (Optim.converged(res))
        sol.H = H
        sol.ab_coeff = Optim.minimizer(res)
        sol.hardy_matrix = loc_hardy_matrix
        evaluation!(sol, verbose=false)
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
    wo_sol = NevanlinnaSolver(input_smpl, input_gw, N_real, omega_max, eta, sum_rule, H_max, iter_tol, lambda, verbose=true)

    open("twopeak_wo_opt.dat","w") do f
        for i in 1:wo_sol.reals.N_real
            println(f, "$(Float64(real.(wo_sol.reals.freq[i])))",  "\t", "$(Float64(imag.(wo_sol.reals.val[i]/pi)))")
        end
    end
end
