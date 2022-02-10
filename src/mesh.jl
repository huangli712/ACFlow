#
# Project : Gardenia
# Source  : mesh.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/10
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

"""
    Base.length(lm::LinearMesh)
"""
function Base.length(lm::LinearMesh)
    lm.nmesh
end

"""
    Base.iterate(lm::LinearMesh)
"""
function Base.iterate(lm::LinearMesh)
    iterate(lm.mesh)
end

"""
    Base.iterate(lm::LinearMesh, i::I64)
"""
function Base.iterate(lm::LinearMesh, i::I64)
    iterate(lm.mesh, i::I64)
end

"""
    Base.eachindex(lm::LinearMesh)
"""
function Base.eachindex(lm::LinearMesh)
    eachindex(lm.mesh)
end

"""
    Base.firstindex(lm::LinearMesh)
"""
function Base.firstindex(lm::LinearMesh)
    firstindex(lm.mesh)
end

"""
    Base.lastindex(lm::LinearMesh)
"""
function Base.lastindex(lm::LinearMesh)
    lastindex(lm.mesh)
end

"""
    Base.getindex(lm::LinearMesh, ind::I64)
"""
function Base.getindex(lm::LinearMesh, ind::I64)
    @assert 1 ≤ ind ≤ lm.nmesh
    return lm.mesh[ind]
end

"""
    Base.getindex(lm::LinearMesh, I::UnitRange{I64})
"""
function Base.getindex(lm::LinearMesh, I::UnitRange{I64})
    @assert checkbounds(Bool, lm.mesh, I)
    lI = length(I)
    X = similar(lm.mesh, lI)
    if lI > 0
        unsafe_copyto!(X, 1, lm.mesh, first(I), lI)
    end
    return X
end

#=
### *Struct : TangentMesh*
=#

"""
    TangentMesh(nmesh::I64, wmin::F64, wmax::F64, f1::F64 = 2.1)
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

"""
    Base.length(tm::TangentMesh)
"""
function Base.length(tm::TangentMesh)
    tm.nmesh
end

"""
    Base.iterate(tm::TangentMesh)
"""
function Base.iterate(tm::TangentMesh)
    iterate(tm.mesh)
end

"""
    Base.iterate(tm::TangentMesh, i::I64)
"""
function Base.iterate(tm::TangentMesh, i::I64)
    iterate(tm.mesh, i)
end

"""
    Base.eachindex(tm::TangentMesh)
"""
function Base.eachindex(tm::TangentMesh)
    eachindex(tm.mesh)
end

"""
    Base.firstindex(tm::TangentMesh)
"""
function Base.firstindex(tm::TangentMesh)
    firstindex(tm.mesh)
end

"""
    Base.lastindex(tm::TangentMesh)
"""
function Base.lastindex(tm::TangentMesh)
    lastindex(tm.mesh)
end

"""
    Base.getindex(tm::TangentMesh, ind::I64)
"""
function Base.getindex(tm::TangentMesh, ind::I64)
    @assert 1 ≤ ind ≤ tm.nmesh
    return tm.mesh[ind]
end

"""
    Base.getindex(tm::TangentMesh, I::UnitRange{I64})
"""
function Base.getindex(tm::TangentMesh, I::UnitRange{I64})
    @assert checkbounds(Bool, tm.mesh, I)
    lI = length(I)
    X = similar(tm.mesh, lI)
    if lI > 0
        unsafe_copyto!(X, 1, tm.mesh, first(I), lI)
    end
    return X
end
