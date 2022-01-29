#
# Project : Gardenia
# Source  : mesh.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/01/29
#

abstract type AbstractMesh end

mutable struct UniformMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

function UniformMesh(nmesh::I64, wmin::F64, wmax::F64)
    @assert nmesh ≥ 1
    @assert wmax > wmin

    mesh = collect(LinRange(wmin, wmax, nmesh))
    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return UniformMesh(nmesh, wmax, wmin, mesh, weight)
end

function Base.eachindex(um::UniformMesh)
    eachindex(um.mesh)
end

function Base.firstindex(um::UniformMesh)
    firstindex(um.mesh)
end

function Base.lastindex(um::UniformMesh)
    lastindex(um.mesh)
end

function Base.getindex(um::UniformMesh, ind::I64)
    @assert 1 ≤ ind ≤ um.nmesh
    return um.mesh[ind]
end

function Base.getindex(um::UniformMesh, I::UnitRange{I64})
    @assert checkbounds(Bool, um.mesh, I)
    lI = length(I)
    X = similar(um.mesh, lI)
    if lI > 0
        unsafe_copyto!(X, 1, um.mesh, first(I), lI)
    end
    return X
end

mutable struct NonUniformMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

function Base.getindex(num::NonUniformMesh, ind::I64)
    @assert 1 ≤ ind ≤ num.nmesh
    return num.mesh[ind]
end


function make_mesh()
    mesh = get_c("mesh")
    nmesh = get_c("nmesh")
    wmax = get_c("wmax")
    wmin = get_c("wmin")

    if mesh == "uniform"
        return UniformMesh(nmesh, wmin, wmax)
    else
        return NonUniformMesh(nmesh, wmin, wmax)
    end
end






function trapz(x::Vector, y::Vector, uniform::Bool = false)
    if uniform
        h = x[2] - x[1]
        _sum = sum(y[2:end-1])
        value = (h / 2.0) * (y[1] + y[end] + 2.0 * _sum)
    else
        len = length(x)
        value = 0.0
        for i = 1:len-1
            value = value + (y[i] + y[i+1]) * (x[i+1] - x[i])
        end
        value = value / 2.0    
    end

    return value
end

function trapz(x::UniformMesh, y::Vector)
    value = dot(x.weight, y)
    return value
end