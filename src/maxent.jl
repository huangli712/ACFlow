#
# Project : Gardenia
# Source  : maxent.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/13
#

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

struct GreenData <: AbstractData
    value :: Vector{C64}
    error :: Vector{F64}
    var  :: Vector{F64}
    cov  :: Array{F64,2}
    ucov :: Array{F64,2}
end

struct MaxEntGrid
    wmesh :: Vector{F64}
end

mutable struct MaxEntContext
    E :: Vector{F64}
    kernel :: Array{C64,2}
end

struct FermionicMatsubaraGrid <: AbstractGrid
    grid :: Vector{F64}
end

function read_data!(::Type{FermionicMatsubaraGrid})
    grid  = F64[] 
    value = C64[]
    error = F64[]
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

    return FermionicMatsubaraGrid(grid), GreenData(value, error, var, cov, ucov)
end

function maxent_mesh()
    wmesh = collect(LinRange(-5.0, 5.0, 501))

    test = (wmesh[2:end] + wmesh[1:end-1]) / 2.0
    pushfirst!(test, wmesh[1])
    push!(test, wmesh[end])
    dw = diff(test)

    return MaxEntGrid(wmesh)
end

function maxent_model(g::MaxEntGrid)
    len = length(g.wmesh)
    model = ones(F64, len) / 10.0
    return model
end

function maxent_kernel(mesh::MaxEntGrid, ω::FermionicMatsubaraGrid)
    niw = length(ω.grid)
    nw = length(mesh.wmesh)
    #@show niw, nw
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
    #@show E
    #@show mesh, ω
    kernel = maxent_kernel(mesh, ω)
    kernel = G.ucov' * kernel

    return MaxEntContext(E, kernel)
end

println("hello")
ω, G = read_data!(FermionicMatsubaraGrid)
mesh = maxent_mesh()
model = maxent_model(mesh)
maxent_init(G, mesh, ω)


