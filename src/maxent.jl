#
# Project : Gardenia
# Source  : maxent.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/24
#

module MaxEnt

using LsqFit
using Einsum
using LinearAlgebra

import ..ACFlow: I64, F64, C64
import ..ACFlow: FermionicMatsubaraGrid, BosonicMatsubaraGrid
import ..ACFlow: AbstractMesh, UniformMesh
import ..ACFlow: RawData, GreenData
import ..ACFlow: make_grid, make_mesh, make_model, make_kernel, make_data
import ..ACFlow: make_singular_space
import ..ACFlow: get_c
import ..ACFlow: secant, trapz

mutable struct MaxEntContext
    mesh :: AbstractMesh
    model
    imdata :: Vector{F64}
    E :: Vector{F64}
    kernel :: Array{F64,2}
    d2chi2
    W2
    W3
    Bm
    n_sv
    V_svd
end

function solve(rd::RawData)
    mec = maxent_init(rd)
    maxent_run(mec)
end

function maxent_init(rd::RawData)
    G = make_data(rd)
    grid = make_grid(rd)
    mesh = make_mesh()

    model = make_model(mesh)

    kernel = make_kernel(mesh, grid)
    U_svd, V_svd, S_svd = make_singular_space(kernel)
    n_sv = length(S_svd)

    E = 1.0 ./ G.var
    imdata = G.value

    niw = mesh.nmesh
    dw = mesh.weight
    W2 = zeros(F64, n_sv, niw)

    @einsum W2[m,l] = E[k] * U_svd[k,m] * S_svd[m] * U_svd[k,n] * S_svd[n] * V_svd[l,n] * dw[l] * model[l]
    A = reshape(W2, (n_sv, 1, niw))
    B = reshape(V_svd', (1, n_sv, niw))
    W3 = zeros(F64, n_sv, n_sv, niw)
    for i = 1:niw
        W3[:,:,i] = A[:,:,i] * B[:,:,i]
    end

    Bm = zeros(F64, n_sv)

    t0 = time_ns()
    @einsum Bm[m] = S_svd[m] * U_svd[k,m] * E[k] * imdata[k]
    t1 = time_ns()
    println(t1 - t0)

    d2chi2 = zeros(F64, niw, niw)
    @einsum d2chi2[i,j] = dw[i] * dw[j] * kernel[k,i] * kernel[k,j] * E[k]

    return MaxEntContext(mesh, model, imdata, E, kernel, d2chi2, W2, W3, Bm, n_sv, V_svd)
end

function maxent_run(mec::MaxEntContext)
    maxent_bryan(mec)
    #maxent_historic(mec)
    #maxent_classic(mec)
    #maxent_chi2kink(mec)
end

function maxent_historic(mec::MaxEntContext)
    println("Solving")
    alpha = 10 ^ 6.0
    ustart = zeros(F64, mec.n_sv)
    optarr = []
    converged = false
    conv = 0.0
    mesh = mec.mesh

    use_bayes = false
    niw = 20
    while conv < 1.0
        o = maxent_optimize(mec, alpha, ustart, mesh, use_bayes)
        push!(optarr, o)
        alpha = alpha / 10.0
        conv = niw / o[:chi2]
        #@show alpha, conv
    end

    ustart = optarr[end-1][:u_opt]
    alpha = optarr[end][:alpha]
    #@show ustart
    #@show alpha

    function root_fun(_alpha, _u)
        res = maxent_optimize(mec, _alpha, _u, mesh, use_bayes)
        _u[:] = copy(res[:u_opt])
        #@show niw / res[:chi2] - 1.0
        #@show _u
        return niw / res[:chi2] - 1.0
    end

    alpha_opt = secant(root_fun, alpha, ustart)
    #@show ustart, alpha_opt

    sol = maxent_optimize(mec, alpha_opt, ustart, mesh, use_bayes)

    A_opt = sol[:A_opt]
    niw = mesh.nmesh
    open("historic.data", "w") do fout
        for i = 1:niw
            println(fout, mesh[i], " ", A_opt[i])
        end
    end
end

function maxent_classic(mec::MaxEntContext)
    println("Solving")
    optarr = []
    alpha = 10 ^ 6.0
    ustart = zeros(F64, mec.n_sv)
    converged = false
    conv = 0.0
    mesh = mec.mesh

    use_bayes = true
    while conv < 1.0
        o = maxent_optimize(mec, alpha, ustart, mesh, use_bayes)
        push!(optarr, o)
        alpha = alpha / 10.0
        conv = o[:convergence]
        ustart = copy(o[:u_opt])
    end
    #@show ustart, alpha

    bayes_conv = [x[:convergence] for x in optarr]
    alpharr = [x[:alpha] for x in optarr]
    #@show bayes_conv
    #@show alpharr
    
    exp_opt = ( log10(alpharr[end]) - log10(alpharr[end-1]) ) / ( log10(bayes_conv[end]) - log10(bayes_conv[end-1]) )
    exp_opt = log10(alpharr[end-1]) - log10(bayes_conv[end-1]) * exp_opt
    alpha_opt = 10 ^ exp_opt
    #@show exp_opt, alpha_opt

    ustart = optarr[end-1][:u_opt]

    function root_fun(_alpha, _u)
        res = maxent_optimize(mec, _alpha, _u, mesh, use_bayes)
        _u[:] = copy(res[:u_opt])
        #@show _u
        return res[:convergence] - 1.0
    end

    alpha_opt = secant(root_fun, alpha_opt, ustart)
    #@show alpha_opt
    #@show ustart

    sol = maxent_optimize(mec, alpha_opt, ustart, mesh, use_bayes)

    A_opt = sol[:A_opt]
    niw = mesh.nmesh
    open("classic.data", "w") do fout
        for i = 1:niw
            println(fout, mesh[i], " ", A_opt[i])
        end
    end
end

function maxent_bryan(mec::MaxEntContext)
    println("hehe")
    alpha = 500
    ustart = zeros(F64, mec.n_sv)
    #@show ustart
    mesh = mec.mesh

    maxprob = 0.0
    optarr = []
    use_bayes = true
    while true
        o = maxent_optimize(mec, alpha, ustart, mesh, use_bayes)
        ustart = copy(o[:u_opt])
        push!(optarr, o)
        alpha = alpha / 1.1
        probability = o[:probability]
        if probability > maxprob
            maxprob = probability
        elseif probability < 0.01 * maxprob
            break
        end
    end

    alpharr = []
    probarr = []
    specarr = []

    for i in eachindex(optarr)
        push!(alpharr, optarr[i][:alpha])
        push!(probarr, optarr[i][:probability])
        push!(specarr, optarr[i][:A_opt])
    end
    
    probarr = probarr ./ (-trapz(alpharr, probarr))
    len = length(probarr)
    niw = mesh.nmesh

    Spectra = zeros(F64, len, niw)
    for i = 1:len
        Spectra[i,:] = specarr[i] * probarr[i]
    end

    A_opt = zeros(F64, niw)
    for i = 1:niw
        A_opt[i] = -trapz(alpharr, Spectra[:,i])
    end

    open("mem.data", "w") do fout
        for i = 1:niw
            println(fout, mesh[i], " ", A_opt[i])
        end
    end
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
    ustart = zeros(F64, mec.n_sv)
    mesh = mec.mesh

    use_bayes = false
    while true
        o = maxent_optimize(mec, alpha, ustart, mesh, use_bayes)
        push!(optarr, o)
        ustart = copy(o[:u_opt])
        push!(chis, o[:chi2])
        push!(alphas, alpha)

        alpha = alpha / 10.0
        if alpha < alpha_end
            break
        end
    end

    #@show alphas
    #@show chis

    function fitfun(x, p)
        #@show coefs
        return @. p[1] + p[2] / (1.0 + exp(-p[4] * (x - p[3])))
    end

    good_numbers = isfinite.(chis)
    #@show log10.(alphas[good_numbers])
    fit = curve_fit(fitfun, log10.(alphas[good_numbers]), log10.(chis[good_numbers]), [0.0, 5.0, 2.0, 0.0])
    a, b, c, d = fit.param

    a_opt = c - fit_position / d
    alpha_opt = 10.0 ^ a_opt

    closest_idx = argmin( abs.( log10.(alphas) .- a_opt ) )
    #@show alpha_opt, closest_idx
    ustart = optarr[closest_idx][:u_opt]
    #@show ustart

    sol = maxent_optimize(mec, alpha_opt, ustart, mesh, use_bayes)

    A_opt = sol[:A_opt]
    niw = mesh.nmesh
    open("chi2kink.data", "w") do fout
        for i = 1:niw
            println(fout, mesh[i], " ", A_opt[i])
        end
    end
end

function maxent_optimize(mec::MaxEntContext, alpha, ustart, mesh::UniformMesh, use_bayes::Bool)
    max_hist = 1

    #@show mec.n_sv, alpha, ustart
    solution, nfev = newton(mec, alpha, ustart, max_hist, function_and_jacobian)
    u_opt = copy(solution)
    A_opt = singular_to_realspace_diag(mec, solution) 
    #@show u_opt
    #@show A_opt
    entr = calc_entropy(mec, A_opt, u_opt, mesh)
    #@show entr
    #@show size(mec.kernel)
    chisq = real(calc_chi2(mec, A_opt, mesh))
    norm = trapz(mesh, A_opt)
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
        ng, tr, conv = calc_bayes(mec, A_opt, entr, alpha, mesh)
        #@show ng, tr, conv
        #error()
        result_dict[:n_good] = ng
        result_dict[:trace] = tr
        result_dict[:convergence] = conv

        probability = calc_posterior(mec, A_opt, alpha, entr, chisq, mesh)
        #@show probability
        result_dict[:probability] = probability
    end

    println("log10(alpha) = $(log10(alpha)) chi2 = $chisq S = $entr nfev = $nfev norm = $norm")

    return result_dict
end

function calc_entropy(mec::MaxEntContext, A, u, mesh::UniformMesh)
    f = A - mec.model - A .* (mec.V_svd * u)
    #@show f
    trapz(mesh, f)
end

function calc_chi2(mec::MaxEntContext, A, mesh::UniformMesh)
    ndim, _ = size(mec.kernel)

    T = zeros(C64, ndim)
    for i = 1:ndim
        T[i] = mec.imdata[i] - trapz(mesh, mec.kernel[i,:] .* A)
        #@show i, T[i]
    end

    value = sum(mec.E .* T .^ 2.0)
    return value
    #mec.kernel * reshape(A, (1,length(A)))
end

function calc_bayes(mec::MaxEntContext, A, entr, alpha, mesh::UniformMesh)
    T = sqrt.(A ./ mesh.weight)
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

function calc_posterior(mec::MaxEntContext, A, alpha, entr, chisq, mesh::UniformMesh)
    T = sqrt.(A ./ mesh.weight)
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

    n_sv, niw = size(mec.W2)
    term_1 = zeros(F64, n_sv)
    term_2 = zeros(F64, n_sv, n_sv)
    for i = 1:n_sv
        term_1[i] = dot(mec.W2[i,:], w)
        #@show i, term_1[i]
    end

    for i = 1:n_sv
        for j = 1:n_sv
            term_2[i,j] = dot(mec.W3[i,j,:], w)
            #@show i, j, term_2[i,j]
        end
    end

    f = alpha * u + term_1 - mec.Bm
    J = alpha * diagm(ones(n_sv)) + term_2
    #@show f
    #@show J
    return f, J
end

function newton(mec::MaxEntContext, alpha, ustart, max_hist, function_and_jacobian)
    max_hist = 1
    max_iter = 20000
    opt_size = mec.n_sv
    mixing = 0.5
    props = []
    res = []
    push!(props, ustart)

    f, J = function_and_jacobian(mec, props[1], alpha)
    initial_result = iteration_function(props[1], f, J)
    #@show initial_result
    #@show f, J
    #error()

    push!(res, initial_result)
    #@show res

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