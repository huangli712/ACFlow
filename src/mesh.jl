#
# Project : Gardenia
# Source  : mesh.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/04/26
#

#=
### *Struct : LinearMesh*
=#

"""
    LinearMesh(nmesh::I64, wmin::F64, wmax::F64)

A constructor for the LinearMesh struct.

See also: [`LinearMesh`](@ref).
"""
function LinearMesh(nmesh::I64, wmin::F64, wmax::F64)
    @assert nmesh ≥ 1
    @assert wmax > wmin

    mesh = collect(LinRange(wmin, wmax, nmesh))
    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return LinearMesh(nmesh, wmax, wmin, mesh, weight)
end

"""
    LinearMesh(mesh::Vector{F64})

Create a LinearMesh struct from a standard Vector. Be careful, `mesh`
must be equidistant.

See also: [`LinearMesh`](@ref).
"""
function LinearMesh(mesh::Vector{F64})
    nmesh = length(mesh)

    wmin = mesh[1]
    wmax = mesh[end]
    @assert wmax > wmin

    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return LinearMesh(nmesh, wmax, wmin, mesh, weight)
end

#=
### *Struct : TangentMesh*
=#

"""
    TangentMesh(nmesh::I64, wmin::F64, wmax::F64, f1::F64 = 2.1)

A constructor for the TangentMesh struct.

See also: [`TangentMesh`](@ref).
"""
function TangentMesh(nmesh::I64, wmin::F64, wmax::F64, f1::F64 = 2.1)
    @assert nmesh ≥ 1
    @assert wmax > 0.0 > wmin
    @assert wmax == abs(wmin)
    @assert f1 > 0.0

    mesh = collect(LinRange(-π / f1, π / f1, nmesh))
    mesh = wmax * tan.(mesh) / tan(π / f1)
    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return TangentMesh(nmesh, wmax, wmin, mesh, weight)
end

#=
### *Struct : LorentzMesh*
=#

"""
    LorentzMesh(nmesh::I64, wmin::F64, wmax::F64, f1::F64 = 2.1)

A constructor for the LorentzMesh struct.

See also: [`LorentzMesh`](@ref).
"""
function LorentzMesh()
end

#=
### *Common Interface*
=#

#=
*Remarks* :

Here we overload some essential functions in the Base module to support
the basic operations for the following meshes:

* LinearMesh
* TangentMesh
* LorentzMesh
* HalfLorentzMesh

With the help of these functions, you can easily visit the elements in
the mesh objective.
=#

"""
    Base.length(am::AbstractMesh)

Return number of mesh points in a Mesh-like struct.
"""
function Base.length(am::AbstractMesh)
    am.nmesh
end

"""
    Base.iterate(am::AbstractMesh)

Advance the iterator of a Mesh-like struct to obtain the next mesh point.
"""
function Base.iterate(am::AbstractMesh)
    iterate(am.mesh)
end

"""
    Base.iterate(am::AbstractMesh, i::I64)

This is the key method that allows a Mesh-like struct to be iterated,
yielding a sequences of mesh points.
"""
function Base.iterate(am::AbstractMesh, i::I64)
    iterate(am.mesh, i::I64)
end

"""
    Base.eachindex(am::AbstractMesh)

Create an iterable object for visiting each index of a Mesh-like struct.
"""
function Base.eachindex(am::AbstractMesh)
    eachindex(am.mesh)
end

"""
    Base.firstindex(am::AbstractMesh)

Return the first index of a Mesh-like struct.
"""
function Base.firstindex(am::AbstractMesh)
    firstindex(am.mesh)
end

"""
    Base.lastindex(am::AbstractMesh)

Return the last index of a Mesh-like struct.
"""
function Base.lastindex(am::AbstractMesh)
    lastindex(am.mesh)
end

"""
    Base.getindex(am::AbstractMesh, ind::I64)

Retrieve the value(s) stored at the given key or index within a
Mesh-like struct.
"""
function Base.getindex(am::AbstractMesh, ind::I64)
    @assert 1 ≤ ind ≤ am.nmesh
    return am.mesh[ind]
end

"""
    Base.getindex(am::AbstractMesh, I::UnitRange{I64})

Return a subset of a Mesh-like struct as specified by `I`.
"""
function Base.getindex(am::AbstractMesh, I::UnitRange{I64})
    @assert checkbounds(Bool, am.mesh, I)
    lI = length(I)
    X = similar(am.mesh, lI)
    if lI > 0
        unsafe_copyto!(X, 1, am.mesh, first(I), lI)
    end
    return X
end
