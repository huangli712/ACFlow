#
# Project : Gardenia
# Source  : mesh.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/08/26
#

#=
### *Struct : LinearMesh*
=#

"""
    LinearMesh(nmesh::I64, wmin::T, wmax::T) where {T}

A constructor for the LinearMesh struct, which is announced in
`src/types.jl`.

### Arguments
* nmesh -> Number of mesh points.
* wmin  -> Left boundary of the mesh.
* wmax  -> Right boundary of the mesh.

### Returns
* lm -> A LinearMesh struct.

See also: [`LinearMesh`](@ref).
"""
function LinearMesh(nmesh::I64, wmin::T, wmax::T) where {T}
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

#=
### *Struct : TangentMesh*
=#

"""
    TangentMesh(nmesh::I64, wmin::T, wmax::T, 𝑝::T = 2.1) where {T}

A constructor for the TangentMesh struct, which is announced in
`src/types.jl`.

### Arguments
* nmesh -> Number of mesh points.
* wmin  -> Left boundary of the mesh.
* wmax  -> Right boundary of the mesh.
* 𝑝     -> A customized parameter.

### Returns
* tm -> A TangentMesh struct.

See also: [`TangentMesh`](@ref).
"""
function TangentMesh(nmesh::I64, wmin::T, wmax::T, 𝑝::T = 2.1) where {T}
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
    LorentzMesh(nmesh::I64, wmin::T, wmax::T, 𝑝::T = 0.01) where {T}

A constructor for the LorentzMesh struct, which is announced in
`src/types.jl`. The algorithm for generating a lorentzian mesh
is taken from:

* https://github.com/CQMP/Maxent.

### Arguments
* nmesh -> Number of mesh points.
* wmin  -> Left boundary of the mesh.
* wmax  -> Right boundary of the mesh.
* 𝑝     -> A customized parameter.

### Returns
* lm -> A LorentzMesh struct.

See also: [`LorentzMesh`](@ref).
"""
function LorentzMesh(nmesh::I64, wmin::T, wmax::T, 𝑝::T = 0.01) where {T}
    @assert nmesh ≥ 1
    @assert wmax > 0.0 > wmin
    @assert wmax == abs(wmin)
    @assert 1.0 > 𝑝 > 0.0

    temp = zeros(T, nmesh)
    mesh = zeros(T, nmesh)

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
    HalfLorentzMesh(nmesh::I64, wmax::T, 𝑝::T = 0.01) where {T}

A constructor for the HalfLorentzMesh struct, which is announced
in `src/types.jl`. The algorithm for generating a half-lorentzian
mesh is taken from:

* https://github.com/CQMP/Maxent.

### Arguments
* nmesh -> Number of mesh points.
* wmax  -> Right boundary of the mesh (wmin ≡ 0.0).
* 𝑝     -> A customized parameter.

### Returns
* hm -> A HalfLorentzMesh struct.

See also: [`HalfLorentzMesh`](@ref).
"""
function HalfLorentzMesh(nmesh::I64, wmax::T, 𝑝::T = 0.01) where {T}
    @assert nmesh ≥ 1
    @assert wmax > 0.0
    @assert 1.0 > 𝑝 > 0.0

    wmin::T = 0.0
    temp = zeros(T, nmesh)
    mesh = zeros(T, nmesh)

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
### *Struct : DynamicMesh*
=#

"""
    DynamicMesh(mesh::Vector{T}) where {T}

A constructor for the DynamicMesh struct, which is announced in
`src/types.jl`. The δ peaks in stochastic analytic continuation methods
or poles in stochastic pole expansion method could be placed in this mesh.
This mesh should not be used to define the spectrum.

### Arguments
* mesh -> Usually a mesh from file `fmesh.inp`. See util/gmesh.jl.

### Returns
* dm -> A DynamicMesh struct.

See also: [`DynamicMesh`](@ref).
"""
function DynamicMesh(mesh::Vector{T}) where {T}
    nmesh = length(mesh)

    wmin = mesh[1]
    wmax = mesh[end]
    @assert wmax > wmin

    weight = (mesh[2:end] + mesh[1:end-1]) / 2.0
    pushfirst!(weight, mesh[1])
    push!(weight, mesh[end])
    weight = diff(weight)

    return DynamicMesh(nmesh, wmax, wmin, mesh, weight)
end

#=
### *Common Interface*
=#

#=
*Remarks* :

Here we overload some essential functions in the Base module to support
the basic operations for the following types of meshes:

* LinearMesh
* TangentMesh
* LorentzMesh
* HalfLorentzMesh
* DynamicMesh

With the help of these functions, we can easily access the mesh's elements.
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

"""
    nearest(am::AbstractMesh, r::F64)

Given a position `r` (0.0 ≤ r ≤ 1.0), and return the index of the nearest
point in the mesh `am`.

### Arguments
See above explanations.

### Returns
See above explanations.

### Examples
```julia
am = LinearMesh(1001, -10.0, 10.0)
pos = nearest(am, 0.2) # pos = 201
println(am[pos]) # -6.0
```

See also: [`AbstractMesh`](@ref).
"""
function nearest(am::AbstractMesh, r::F64)
    # Check r and evaluate the corresponding value
    @assert 0.0 ≤ r ≤ 1.0
    val = am.wmin + (am.wmax - am.wmin) * r

    # Try to locate val in the mesh by using the bisection algorithm
    left = 1
    right = length(am)
    @assert am[left] ≤ val ≤ am[right]

    while right - left ≥ 2
        mid = round(I64, (left + right) / 2)
        if val < am[mid]
            right = mid
        else
            left = mid
        end
    end

    # Well, now we have the left and right boundaries. We should return
    # the closer one.
    if am[right] - val > val - am[left]
        return left
    else
        return right
    end
end
