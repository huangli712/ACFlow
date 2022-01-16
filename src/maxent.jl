#
# Project : Gardenia
# Source  : maxent.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/14
#

using Einsum
using Statistics
using LinearAlgebra

@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

const I64 = Int64
const F64 = Float64
const C64 = ComplexF64

abstract type AbstractData end
abstract type AbstractGrid end

mutable struct GreenData <: AbstractData
    value :: Vector{C64}
    error :: Vector{F64}
    imdata :: Vector{F64}
    var  :: Vector{F64}
    cov  :: Array{F64,2}
    ucov :: Array{F64,2}
end

struct MaxEntGrid
    wmesh :: Vector{F64}
    dw :: Vector{F64}
end

mutable struct MaxEntContext
    E :: Vector{F64}
    kernel :: Array{C64,2}
    d2chi2
    W2
    W3
    n_sv
    V_svd
    Evi
    model
    imdata
end

struct FermionicMatsubaraGrid <: AbstractGrid
    grid :: Vector{F64}
end

function read_data!(::Type{FermionicMatsubaraGrid})
    grid  = F64[] 
    value = C64[]
    error = F64[]
    imdata = F64[]
    var   = F64[]

    niw = 10
    #
    open("green.data", "r") do fin
        for i = 1:niw
            arr = parse.(F64, line_to_array(fin))
            push!(grid, arr[1])
            push!(value, arr[2] + arr[3] * im)
            push!(error, 0.0001)
            push!(var, 0.0001^2)
        end
    end
    cov = diagm(var)
    ucov = diagm(ones(niw))
    value = ucov' * value

    return FermionicMatsubaraGrid(grid), GreenData(value, error, imdata, var, cov, ucov)
end

function maxent_mesh()
    wmesh = collect(LinRange(-5.0, 5.0, 501))

    test = (wmesh[2:end] + wmesh[1:end-1]) / 2.0
    pushfirst!(test, wmesh[1])
    push!(test, wmesh[end])
    dw = diff(test)

    return MaxEntGrid(wmesh, dw)
end

function maxent_model(g::MaxEntGrid)
    len = length(g.wmesh)
    model = ones(F64, len) / 10.0
    return model
end

function maxent_kernel(mesh::MaxEntGrid, ω::FermionicMatsubaraGrid)
    niw = length(ω.grid)
    nw = length(mesh.wmesh)

    kernel = zeros(C64, niw, nw)
    for i = 1:nw
        for j = 1:niw
            kernel[j,i] = 1.0 / (im * ω.grid[j] - mesh.wmesh[i])
        end
    end

    return kernel
end

function maxent_init(G::GreenData, mesh::MaxEntGrid, ω::FermionicMatsubaraGrid)
    E = 1.0 ./ G.var
    kernel = maxent_kernel(mesh, ω)
    kernel = G.ucov' * kernel

    # Only for fermionic frequency
    G.var = vcat(G.var, G.var)
    G.imdata = vcat(real(G.value), imag(G.value))
    kernel = vcat(real(kernel), imag(kernel))
    #@show size(kernel)
    #@show kernel[:,1]
    #@show kernel[:,33]
    #@show kernel[:,end]
    E = vcat(E, E)

    F = svd(kernel)
    U, S, V = F
    #@show U[:,1]
    #@show U[:,11]
    #@show U[:,end]
    #@show S
    n_sv = count( x -> x ≥ 1e-10, S)
    U_svd = U[:, 1:n_sv]
    V_svd = V[:, 1:n_sv]
    Xi_svd = S[1:n_sv]
    #@show size(V)
    #@show n_sv
    #@show size(U_svd), size(V_svd), size(Xi_svd)

    open("svd.data", "r") do fin
        for i = 1:20
            for j = 1:20
                U_svd[i,j] = parse(F64, line_to_array(fin)[3])
            end
        end
        readline(fin)

        for i = 1:501
            for j = 1:20
                V_svd[i,j] = parse(F64, line_to_array(fin)[3])
            end
        end
        readline(fin)
        for i = 1:20
            Xi_svd[i] = parse(F64, line_to_array(fin)[2])
        end
    end

    #@show length(mesh.wmesh)
    #error()
    
    niw = length(mesh.wmesh)
    W2 = zeros(F64, n_sv, niw)
    dw = mesh.dw
    model = maxent_model(mesh)
    @einsum W2[m,l] = E[k] * U_svd[k,m] * Xi_svd[m] * U_svd[k,n] * Xi_svd[n] * V_svd[l,n] * dw[l] * model[l]
    A = reshape(W2, (n_sv, 1, niw))
    B = reshape(V_svd', (1, n_sv, niw))
    W3 = zeros(F64, n_sv, n_sv, niw)
    for i = 1:niw
        W3[:,:,i] = A[:,:,i] * B[:,:,i]
    end
    #@show W3[:,:,end]

    Evi = zeros(F64, n_sv)
    imdata = G.imdata
    @einsum Evi[m] = Xi_svd[m] * U_svd[k,m] * E[k] * imdata[k]
    #@show Evi

    d2chi2 = zeros(F64, niw, niw)
    @einsum d2chi2[i,j] = dw[i] * dw[j] * kernel[k,i] * kernel[k,j] * E[k]
    #@show d2chi2[:,end]

    #@show size(kernel)
    return MaxEntContext(E, kernel, d2chi2, W2, W3, n_sv, V_svd, Evi, model, imdata)
end

function maxent_run_bryan(mec::MaxEntContext, mesh::MaxEntGrid)
    println("hehe")
    alpha = 500
    ustart = zeros(F64, mec.n_sv)
    #@show ustart

    maxprob = 0.0
    optarr = []
    while true
        o = maxent_optimize(mec, alpha, ustart, mesh)
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
    
    #probarr = probarr ./ trapz(alpharr, probarr)
    @show probarr 

end

function function_and_jacobian(mec::MaxEntContext, u, alpha)
    #@show length(u), alpha
    v = mec.V_svd * u
    w = exp.(v)
    #@show w
    #@show size(mec.W2), size(mec.W3)

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

    f = alpha * u + term_1 - mec.Evi
    J = alpha * diagm(ones(n_sv)) + term_2
    #@show f
    #@show J
    return f, J
end

function singular_to_realspace_diag(mec::MaxEntContext, u)
    return mec.model .* exp.(mec.V_svd * u)
end

function entropy_pos(mec::MaxEntContext, A, u, mesh::MaxEntGrid)
    f = A - mec.model - A .* (mec.V_svd * u)
    #@show f
    trapz(mesh.wmesh, f)
end

function trapz(x, y)
    #@show x
    #@show y
    h = x[2] - x[1]
    sum = 0.0
    for i = 2:length(x)-1
        sum = sum + y[i]
    end
    value = (h / 2.0) * (y[1] + y[end] + 2.0 * sum)
        
    return value
end

function chi2(mec::MaxEntContext, A, mesh::MaxEntGrid)
    ndim, _ = size(mec.kernel)

    T = zeros(C64, ndim)
    for i = 1:ndim
        T[i] = mec.imdata[i] - trapz(mesh.wmesh, mec.kernel[i,:] .* A)
        #@show i, T[i]
    end

    value = sum(mec.E .* T .^ 2.0)
    return value
    #mec.kernel * reshape(A, (1,length(A)))
end

function bayes_conv(mec::MaxEntContext, A, entr, alpha, mesh::MaxEntGrid)
    T = sqrt.(A ./ mesh.dw)
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

function posterior_probability(mec::MaxEntContext, A, alpha, entr, chisq, mesh::MaxEntGrid)
    T = sqrt.(A ./ mesh.dw)
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

function maxent_optimize(mec::MaxEntContext, alpha, ustart, mesh::MaxEntGrid)
    use_bayes = true
    max_hist = 1

    solution, nfev = newton(mec, alpha, ustart, max_hist, function_and_jacobian)
    u_opt = copy(solution)
    A_opt = singular_to_realspace_diag(mec, solution) 
    #@show u_opt
    #@show A_opt
    entr = entropy_pos(mec, A_opt, u_opt, mesh)
    #@show entr
    #@show size(mec.kernel)
    chisq = real(chi2(mec, A_opt, mesh))
    norm = trapz(mesh.wmesh, A_opt)
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
        ng, tr, conv = bayes_conv(mec, A_opt, entr, alpha, mesh)
        #@show ng, tr, conv
        #error()
        result_dict[:n_good] = ng
        result_dict[:trace] = tr
        result_dict[:convergence] = conv

        probability = posterior_probability(mec, A_opt, alpha, entr, chisq, mesh)
        #@show probability
        result_dict[:probability] = probability
    end

    println("log10(alpha) = $(log10(alpha)) chi2 = $chisq S = $entr nfev = $nfev norm = $norm")

    return result_dict
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

println("hello")
ω, G = read_data!(FermionicMatsubaraGrid)
mesh = maxent_mesh()
mec = maxent_init(G, mesh, ω)
maxent_run_bryan(mec, mesh)


