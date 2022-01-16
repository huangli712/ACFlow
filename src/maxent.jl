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

    @show length(mesh.wmesh)
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

    return MaxEntContext(E, kernel, d2chi2, W2, W3, n_sv, V_svd, Evi)
end

function maxent_run_bryan(mec::MaxEntContext)
    println("hehe")
    alpha = 500
    ustart = zeros(F64, mec.n_sv)
    #@show ustart

    maxent_optimize(mec, alpha, ustart)
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

function maxent_optimize(mec::MaxEntContext, alpha, ustart)
    use_bayes = false
    max_hist = 1

    newton(mec, alpha, ustart, max_hist, function_and_jacobian)

end

function newton(mec::MaxEntContext, alpha, ustart, max_hist, function_and_jacobian)
    max_hist = 1
    max_iter = 20000
    opt_size = mec.n_sv
    props = []
    res = []
    push!(props, ustart)

    f, J = function_and_jacobian(mec, props[1], alpha)
    initial_result = iteration_function(props[1], f, J)
    #@show f, J
    #error()
end

function iteration_function(proposal, function_vector, jacobian_matrix)
    increment = - pinv(jacobian_matrix) * function_vector
    @show increment
end

println("hello")
ω, G = read_data!(FermionicMatsubaraGrid)
mesh = maxent_mesh()
mec = maxent_init(G, mesh, ω)
maxent_run_bryan(mec)


