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
import ..ACFlow: secant, trapz

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

    W₂, W₃, Bₘ, d2chi2 = precompute(Gdata, E, mesh, model, kernel, U_svd, V_svd, S_svd)

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
        conv = sol[:convergence]
    end

    bayes_conv = [x[:convergence] for x in optarr]
    alpharr = [x[:alpha] for x in optarr]
    exp_opt = log10(alpharr[end] / alpharr[end-1])
    exp_opt = exp_opt / log10(bayes_conv[end] / bayes_conv[end-1])
    exp_opt = log10(alpharr[end-1]) - log10(bayes_conv[end-1]) * exp_opt
    alpha = 10.0 ^ exp_opt
    ustart = optarr[end-1][:u_opt]

    function root_fun(_alpha, _u)
        res = maxent_optimize(mec, _alpha, _u, use_bayes)
        @. _u = res[:u_opt]
        return res[:convergence] - 1.0
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
        o = maxent_optimize(mec, alpha, ustart, use_bayes)
        ustart = copy(o[:u_opt])
        push!(optarr, o)
        alpha = alpha / 1.1
        prob = o[:probability]
        if prob > maxprob
            maxprob = prob
        elseif prob < 0.01 * maxprob
            break
        end
    end

    alpharr = map(x->x[:alpha], optarr)
    probarr = map(x->x[:probability], optarr)
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
    alpha_start = 1e9
    alpha_end = 1e-3
    alpha_div = 10.0
    fit_position = 2.5

    alpha = alpha_start
    chis = []
    alphas = []
    optarr = []
    ustart = zeros(F64, length(mec.Bₘ))
    mesh = mec.mesh

    use_bayes = false
    while true
        o = maxent_optimize(mec, alpha, ustart, use_bayes)
        push!(optarr, o)
        ustart = copy(o[:u_opt])
        push!(chis, o[:chi2])
        push!(alphas, alpha)

        alpha = alpha / 10.0
        if alpha < alpha_end
            break
        end
    end

    function fitfun(x, p)
        return @. p[1] + p[2] / (1.0 + exp(-p[4] * (x - p[3])))
    end

    good_numbers = isfinite.(chis)
    fit = curve_fit(fitfun, log10.(alphas[good_numbers]), log10.(chis[good_numbers]), [0.0, 5.0, 2.0, 0.0])
    a, b, c, d = fit.param

    a_opt = c - fit_position / d
    alpha_opt = 10.0 ^ a_opt

    closest_idx = argmin( abs.( log10.(alphas) .- a_opt ) )
    ustart = optarr[closest_idx][:u_opt]

    sol = maxent_optimize(mec, alpha_opt, ustart, use_bayes)

    return optarr, sol
    #A_opt = sol[:A_opt]
    #niw = mesh.nmesh
    #open("chi2kink.data", "w") do fout
    #    for i = 1:niw
    #        println(fout, mesh[i], " ", A_opt[i])
    #    end
    #end
end

function maxent_optimize(mec::MaxEntContext, alpha, ustart, use_bayes::Bool)
    max_hist = 1

    #@show alpha, ustart
    solution, nfev = newton(mec, alpha, ustart, max_hist, function_and_jacobian)
    u_opt = copy(solution)
    A_opt = singular_to_realspace_diag(mec, solution) 
    #@show u_opt
    #@show A_opt
    entr = calc_entropy(mec, A_opt, u_opt)
    #@show entr
    #@show size(mec.kernel)
    chisq = real(calc_chi2(mec, A_opt))
    norm = trapz(mec.mesh, A_opt)
    #error()
    #@show chisq
    #@show norm

    result_dict = Dict{Symbol,Any}()
    result_dict[:u_opt] = u_opt
    result_dict[:A_opt] = A_opt

    result_dict[:alpha] = alpha
    result_dict[:entropy] = entr
    result_dict[:chi2] = chisq
    result_dict[:blacktransform] = nothing
    result_dict[:norm] = norm
    result_dict[:Q] = alpha * entr - 0.5 * chisq
    #@show result_dict[:Q]

    if use_bayes
        ng, tr, conv = calc_bayes(mec, A_opt, entr, alpha)
        #@show ng, tr, conv
        #error()
        result_dict[:n_good] = ng
        result_dict[:trace] = tr
        result_dict[:convergence] = conv

        probability = calc_posterior(mec, A_opt, alpha, entr, chisq)
        #@show probability
        result_dict[:probability] = probability
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

    A = reshape(W₂, (n_svd, 1, nmesh))
    B = reshape(V', (1, n_svd, nmesh))
    for i = 1:nmesh
        W₃[:,:,i] = A[:,:,i] * B[:,:,i]
    end

    @einsum Bₘ[m] = S[m] * U[k,m] * E[k] * Gdata[k]

    @einsum d2chi2[i,j] = weight[i] * weight[j] * kernel[k,i] * kernel[k,j] * E[k]

    return W₂, W₃, Bₘ, d2chi2
end

function calc_entropy(mec::MaxEntContext, A, u)
    f = A - mec.model - A .* (mec.V_svd * u)
    trapz(mec.mesh, f)
end

function calc_chi2(mec::MaxEntContext, A)
    ndim, _ = size(mec.kernel)

    T = zeros(C64, ndim)
    for i = 1:ndim
        T[i] = mec.Gdata[i] - trapz(mec.mesh, mec.kernel[i,:] .* A)
    end

    value = sum(mec.E .* T .^ 2.0)
    return value
end

function calc_bayes(mec::MaxEntContext, A, entr, alpha)
    T = sqrt.(A ./ mec.mesh.weight)
    len = length(T)
    Tl = reshape(T, (len, 1))
    Tr = reshape(T, (1, len))
    #@show size(Tl), size(Tr), size(mec.d2chi2)
    Λ = zeros(F64, len, len)
    for i = 1:len
        for j = 1:len
            Λ[i,j] = Tl[i] * mec.d2chi2[i,j] * Tr[j]
        end
    end
    #@show Λ[2,:]

    lam = eigvals(Hermitian(Λ))
    #@show lam
    ng = -2.0 * alpha * entr
    tr = sum(lam ./ (alpha .+ lam))
    conv = tr / ng
    return ng, tr, conv
end

function calc_posterior(mec::MaxEntContext, A, alpha, entr, chisq)
    T = sqrt.(A ./ mec.mesh.weight)
    len = length(T)
    Tl = reshape(T, (len, 1))
    Tr = reshape(T, (1, len))
    #@show size(Tl), size(Tr), size(mec.d2chi2)
    Λ = zeros(F64, len, len)
    for i = 1:len
        for j = 1:len
            Λ[i,j] = Tl[i] * mec.d2chi2[i,j] * Tr[j]
        end
    end
    lam = eigvals(Hermitian(Λ))

    eig_sum = sum(log.(alpha ./ (alpha .+ lam)))
    log_prob = alpha * entr - 0.5 * chisq + log(alpha) + 0.5 * eig_sum
    return exp(log_prob)
end

function singular_to_realspace_diag(mec::MaxEntContext, u)
    return mec.model .* exp.(mec.V_svd * u)
end

function function_and_jacobian(mec::MaxEntContext, u, alpha)
    v = mec.V_svd * u
    w = exp.(v)

    n_svd, niw = size(mec.W₂)
    term_1 = zeros(F64, n_svd)
    term_2 = zeros(F64, n_svd, n_svd)
    for i = 1:n_svd
        term_1[i] = dot(mec.W₂[i,:], w)
        #@show i, term_1[i]
    end

    for i = 1:n_svd
        for j = 1:n_svd
            term_2[i,j] = dot(mec.W₃[i,j,:], w)
            #@show i, j, term_2[i,j]
        end
    end

    f = alpha * u + term_1 - mec.Bₘ
    J = alpha * diagm(ones(n_svd)) + term_2

    return f, J
end

function newton(mec::MaxEntContext, alpha, ustart, max_hist, function_and_jacobian)
    max_hist = 1
    max_iter = 20000
    opt_size = length(ustart)
    mixing = 0.5
    props = []
    res = []
    push!(props, ustart)

    f, J = function_and_jacobian(mec, props[1], alpha)
    initial_result = iteration_function(props[1], f, J)

    push!(res, initial_result)

    counter = 0
    converged = false
    result = nothing
    while !converged
        n_iter = length(props)
        new_proposal = props[n_iter]
        f_i = res[n_iter] - props[n_iter]
        update = mixing * f_i
        prop = new_proposal + update
        #@show prop
        #error()

        f, J = function_and_jacobian(mec, prop, alpha)
        #@show f
        #@show J
        #error()

        result = iteration_function(prop, f, J)
        #@show result

        push!(props, prop)
        push!(res, result)

        converged = counter > max_iter || maximum( abs.(result - prop) ) < 1.e-4
        counter = counter + 1
        #@show counter, maximum( abs.(result - prop) )
        if any(isnan.(result))
            error("have NaN")
        end
        #error()
    end

    if counter > max_iter
        error("max_iter is too small")
    end

    #@show result, counter
    return result, counter
end

function iteration_function(proposal, function_vector, jacobian_matrix)
    increment = nothing
    try
        increment = - pinv(jacobian_matrix) * function_vector
        #@show increment
    catch
        increment = zeros(F64, length(proposal))
    end
    #@show increment
    step_reduction = 1.0
    significance_limit = 1e-4

    #@. proposal = 1.0
    if any(x -> x > significance_limit, abs.(proposal))
        ratio = abs.(increment ./ proposal)
        #@show ratio
        #filter!(x -> abs(x) > significance_limit, ratio)
        #@show ratio
        max_ratio = maximum( ratio[ abs.(proposal) .> significance_limit ] )
        #@show max_ratio
        if max_ratio > 1.0
            step_reduction = 1.0 / max_ratio
        end
        #@show step_reduction
    end

    result = proposal + step_reduction .* increment
    return result
end

end