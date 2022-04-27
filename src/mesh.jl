#
# Project : Gardenia
# Source  : mesh.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/04/27
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
    #
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
    TangentMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 2.1)

A constructor for the TangentMesh struct.

See also: [`TangentMesh`](@ref).
"""
function TangentMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 2.1)
    @assert nmesh ≥ 1
    @assert wmax > 0.0 > wmin
    @assert wmax == abs(wmin)
    @assert 𝑝 > 0.0

    mesh = collect(LinRange(-π / 𝑝, π / 𝑝, nmesh))
    mesh = wmax * tan.(mesh) / tan(π / 𝑝)
    #
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
    LorentzMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 0.01)

A constructor for the LorentzMesh struct.

See also: [`LorentzMesh`](@ref).
"""
function LorentzMesh(nmesh::I64, wmin::F64, wmax::F64, 𝑝::F64 = 0.01)
    @assert nmesh ≥ 1
    @assert wmax > 0.0 > wmin
    @assert wmax == abs(wmin)
    @assert 1.0 > 𝑝 > 0.0

    temp = zeros(F64, nmesh)
    mesh = zeros(F64, nmesh)

    for i in eachindex(temp)
        f = ( (i - 1) / (nmesh - 1) * (1.0 - 2.0 * 𝑝) + 𝑝 - 0.5 )
        temp[i] = tan(π * f)
    end

    for i in eachindex(mesh)
        # mesh[i] is in [0,1]
        mesh[i] = (temp[i] - temp[1]) / (temp[end] - temp[1])
        #
        # Convert mesh[i] from [0,1] to [wmin,wmax]
        mesh[i] = mesh[i] * 2.0 * wmax - wmax
    end

    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return LorentzMesh(nmesh, wmax, wmin, mesh, weight)
end

#=
### *Struct : HalfLorentzMesh*
=#

"""
    HalfLorentzMesh(nmesh::I64, wmax::F64, 𝑝::F64 = 0.01)

A constructor for the HalfLorentzMesh struct.

See also: [`HalfLorentzMesh`](@ref).
"""
function HalfLorentzMesh(nmesh::I64, wmax::F64, 𝑝::F64 = 0.01)
    @assert nmesh ≥ 1
    @assert wmax > 0.0
    @assert 1.0 > 𝑝 > 0.0

    wmin = 0.0
    temp = zeros(F64, nmesh)
    mesh = zeros(F64, nmesh)

    for i in eachindex(temp)
        f = (i - 2 + nmesh) / ( 2 * nmesh - 3 ) * (1.0 - 2.0 * 𝑝) + 𝑝 - 0.5
        temp[i] = tan(π * f)
    end

    for i in eachindex(mesh)
        # mesh[i] is in [0,1]
        mesh[i] = (temp[i] - temp[1]) / (temp[end] - temp[1])
        #
        # Convert mesh[i] from [0,1] to [0,wmax]
        mesh[i] = mesh[i] * wmax
    end

    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return HalfLorentzMesh(nmesh, wmax, wmin, mesh, weight)
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

With the help of these functions, we can easily visit the elements in
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
