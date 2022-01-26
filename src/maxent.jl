#
# Project : Gardenia
# Source  : maxent.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/26
#

module MaxEnt

using LsqFit
using Einsum
using LinearAlgebra

import ..ACFlow: I64, F64, C64
import ..ACFlow: @cswitch
import ..ACFlow: FermionicMatsubaraGrid, BosonicMatsubaraGrid
import ..ACFlow: AbstractMesh
import ..ACFlow: RawData, GreenData
import ..ACFlow: make_grid, make_mesh, make_model, make_kernel, make_data
import ..ACFlow: make_singular_space
import ..ACFlow: write_spectrum
import ..ACFlow: get_c, get_m
import ..ACFlow: secant, newton, trapz

mutable struct MaxEntContext
    Gdata :: Vector{F64}
    E :: Vector{F64}
    mesh :: AbstractMesh
    model :: Vector{F64}
    kernel :: Array{F64,2}
    V_svd :: Array{F64,2}
    W₂ :: Array{F64,2}
    W₃ :: Array{F64,3}
    Bₘ :: Vector{F64}
    d2chi2 :: Array{F64,2}
end

function solve(rd::RawData)
    mec = maxent_init(rd)
    maxent_run(mec)
end

function maxent_init(rd::RawData)
    G = make_data(rd)
    E = 1.0 ./ G.var
    Gdata = G.value

    grid = make_grid(rd)
    mesh = make_mesh()

    model = make_model(mesh)

    kernel = make_kernel(mesh, grid)
    U_svd, V_svd, S_svd = make_singular_space(kernel)

    @timev W₂, W₃, Bₘ, d2chi2 = precompute(Gdata, E, mesh, model, kernel, U_svd, V_svd, S_svd)

    return MaxEntContext(Gdata, E, mesh, model, kernel, V_svd, W₂, W₃, Bₘ, d2chi2)
end

function maxent_run(mec::MaxEntContext)
    method = get_m("method")

    @cswitch method begin
        @case "historic"
            D, sol = maxent_historic(mec)
            break

        @case "classic"
            D, sol = maxent_classic(mec)
            break
        
        @case "bryan"
            D, sol = maxent_bryan(mec)
            break
        
        @case "chi2kink"
            D, sol = maxent_chi2kink(mec)
            break
    end

    write_spectrum(mec.mesh, sol[:A_opt])
end

function maxent_historic(mec::MaxEntContext)
    println("Using historic algorithm to solve the maximum entropy problem")

    use_bayes = false
    alpha = get_m("alpha")
    n_svd = length(mec.Bₘ)
    ustart = zeros(F64, n_svd)
    optarr = []

    conv = 0.0
    while conv < 1.0
        sol = maxent_optimize(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        alpha = alpha / 10.0
        conv = n_svd / sol[:chi2]
    end

    ustart = optarr[end-1][:u_opt]
    alpha = optarr[end][:alpha]

    function root_fun(_alpha, _u)
        res = maxent_optimize(mec, _alpha, _u, use_bayes)
        @. _u = res[:u_opt]
        return n_svd / res[:chi2] - 1.0
    end
    alpha_opt = secant(root_fun, alpha, ustart)

    sol = maxent_optimize(mec, alpha_opt, ustart, use_bayes)

    return optarr, sol
end

function maxent_classic(mec::MaxEntContext)
    println("Using historic algorithm to solve the maximum entropy problem")

    use_bayes = true
    alpha = get_m("alpha")
    n_svd = length(mec.Bₘ)
    ustart = zeros(F64, n_svd)
    optarr = []

    conv = 0.0
    while conv < 1.0
        sol = maxent_optimize(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        alpha = alpha / 10.0
        @. ustart = sol[:u_opt]
        conv = sol[:conv]
    end

    convarr = [x[:conv] for x in optarr]
    alpharr = [x[:alpha] for x in optarr]
    exp_opt = log10(alpharr[end] / alpharr[end-1])
    exp_opt = exp_opt / log10(convarr[end] / convarr[end-1])
    exp_opt = log10(alpharr[end-1]) - log10(convarr[end-1]) * exp_opt
    alpha = 10.0 ^ exp_opt
    ustart = optarr[end-1][:u_opt]

    function root_fun(_alpha, _u)
        res = maxent_optimize(mec, _alpha, _u, use_bayes)
        @. _u = res[:u_opt]
        return res[:conv] - 1.0
    end
    alpha_opt = secant(root_fun, alpha, ustart)

    sol = maxent_optimize(mec, alpha_opt, ustart, use_bayes)

    return optarr, sol 
end

function maxent_bryan(mec::MaxEntContext)
    println("Using bryan algorithm to solve the maximum entropy problem")

    use_bayes = true
    alpha = get_m("alpha")
    n_svd = length(mec.Bₘ)
    ustart = zeros(F64, n_svd)
    optarr = []

    maxprob = 0.0
    while true
        sol = maxent_optimize(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        alpha = alpha / 1.1
        @. ustart = sol[:u_opt]
        prob = sol[:prob]
        if prob > maxprob
            maxprob = prob
        elseif prob < 0.01 * maxprob
            break
        end
    end

    alpharr = map(x->x[:alpha], optarr)
    probarr = map(x->x[:prob], optarr)
    specarr = map(x->x[:A_opt], optarr)
    probarr = probarr ./ (-trapz(alpharr, probarr))

    nprob = length(probarr)
    nmesh = length(specarr[1])
    A_opt = zeros(F64, nmesh)
    Spectra = zeros(F64, nprob, nmesh)
    for i = 1:nprob
        Spectra[i,:] = specarr[i] * probarr[i]
    end
    for i = 1:nmesh
        A_opt[i] = -trapz(alpharr, Spectra[:,i])
    end

    sol = Dict{Symbol,Any}()
    sol[:A_opt] = A_opt

    return optarr, sol
end

function maxent_chi2kink(mec::MaxEntContext)
    println("Using chi2kink algorithm to solve the maximum entropy problem")
    
    use_bayes = false
    alpha = get_m("alpha")
    n_svd = length(mec.Bₘ)
    ustart = zeros(F64, n_svd)
    optarr = []

    chi_arr = []
    alpharr = []
    while true
        sol = maxent_optimize(mec, alpha, ustart, use_bayes)
        push!(optarr, sol)
        push!(chi_arr, sol[:chi2])
        push!(alpharr, alpha)
        @. ustart = sol[:u_opt]
        alpha = alpha / 10.0
        if alpha < 0.001
            break
        end
    end

    function fitfun(x, p)
        return @. p[1] + p[2] / (1.0 + exp(-p[4] * (x - p[3])))
    end

    good_data = isfinite.(chi_arr)
    fit = curve_fit(fitfun, log10.(alpharr[good_data]), log10.(chi_arr[good_data]), [0.0, 5.0, 2.0, 0.0])
    _, _, c, d = fit.param

    fit_position = 2.5
    a_opt = c - fit_position / d
    closest_idx = argmin( abs.( log10.(alpharr) .- a_opt ) )
    ustart = optarr[closest_idx][:u_opt]
    alpha_opt = 10.0 ^ a_opt

    sol = maxent_optimize(mec, alpha_opt, ustart, use_bayes)

    return optarr, sol
end

function maxent_optimize(mec::MaxEntContext,
                         alpha::F64,
                         ustart::Vector{F64},
                         use_bayes::Bool)
    solution, nfev = newton(f_and_J, ustart, mec, alpha)
    u_opt = copy(solution)
    A_opt = svd_to_real(mec, solution) 
    entr = calc_entropy(mec, A_opt, u_opt)
    chisq = calc_chi2(mec, A_opt)
    norm = trapz(mec.mesh, A_opt)

    result_dict = Dict{Symbol,Any}()
    result_dict[:u_opt] = u_opt
    result_dict[:A_opt] = A_opt

    result_dict[:alpha] = alpha
    result_dict[:entropy] = entr
    result_dict[:chi2] = chisq
    result_dict[:norm] = norm
    result_dict[:Q] = alpha * entr - 0.5 * chisq

    if use_bayes
        ng, tr, conv, prob = calc_bayes(mec, A_opt, entr, chisq, alpha)
        result_dict[:n_good] = ng
        result_dict[:trace] = tr
        result_dict[:conv] = conv
        result_dict[:prob] = prob
    end

    println("log10(alpha) = $(log10(alpha)) chi2 = $chisq S = $entr nfev = $nfev norm = $norm")

    return result_dict
end

function precompute(Gdata::Vector{F64}, E::Vector{F64}, 
                    mesh::AbstractMesh, model::Vector{F64},
                    kernel::Matrix{F64},
                    U::Matrix{F64}, V::Matrix{F64}, S::Vector{F64})
    nmesh = mesh.nmesh
    weight = mesh.weight
    n_svd = length(S)

    W₂ = zeros(F64, n_svd, nmesh)
    W₃ = zeros(F64, n_svd, n_svd, nmesh)
    Bₘ = zeros(F64, n_svd)
    d2chi2 = zeros(F64, nmesh, nmesh)

    @einsum W₂[m,l] = E[k] * U[k,m] * S[m] * U[k,n] * S[n] * V[l,n] * weight[l] * model[l]

    for i = 1:nmesh
        W₃[:,:,i] = W₂[:,i] * (V[i,:])'
    end

    @einsum Bₘ[m] = S[m] * U[k,m] * E[k] * Gdata[k]

    @einsum d2chi2[i,j] = weight[i] * weight[j] * kernel[k,i] * kernel[k,j] * E[k]

    return W₂, W₃, Bₘ, d2chi2
end

function f_and_J(u::Vector{F64}, mec::MaxEntContext, alpha::F64)
    v = mec.V_svd * u
    w = exp.(v)

    n_svd = length(mec.Bₘ)
    term_1 = zeros(F64, n_svd)
    term_2 = zeros(F64, n_svd, n_svd)

    for i = 1:n_svd
        term_1[i] = dot(mec.W₂[i,:], w)
    end

    for j = 1:n_svd
        for i = 1:n_svd
            term_2[i,j] = dot(mec.W₃[i,j,:], w)
        end
    end

    f = alpha * u + term_1 - mec.Bₘ
    J = alpha * diagm(ones(n_svd)) + term_2

    return f, J
end

function svd_to_real(mec::MaxEntContext, u::Vector{F64})
    return mec.model .* exp.(mec.V_svd * u)
end

function calc_entropy(mec::MaxEntContext, A::Vector{F64}, u::Vector{F64})
    f = A - mec.model - A .* (mec.V_svd * u)
    return trapz(mec.mesh, f)
end

function calc_chi2(mec::MaxEntContext, A::Vector{F64})
    ndim, _ = size(mec.kernel)

    T = zeros(F64, ndim)
    for i = 1:ndim
        T[i] = mec.Gdata[i] - trapz(mec.mesh, mec.kernel[i,:] .* A)
    end
    value = sum(mec.E .* T .^ 2.0)

    return value
end

function calc_bayes(mec::MaxEntContext, A::Vector{F64}, ent::F64, chisq::F64, alpha::F64)
    mesh = mec.mesh

    T = sqrt.(A ./ mesh.weight)
    Λ = (T * T') .* mec.d2chi2
    
    lam = eigvals(Hermitian(Λ))
    ng = -2.0 * alpha * ent
    tr = sum(lam ./ (alpha .+ lam))
    conv = tr / ng

    eig_sum = sum(log.(alpha ./ (alpha .+ lam)))
    log_prob = alpha * ent - 0.5 * chisq + log(alpha) + 0.5 * eig_sum

    return ng, tr, conv, exp(log_prob)
end

end