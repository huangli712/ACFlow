#
# Project : Gardenia
# Source  : grid.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/08/26
#

#=
### *Struct : FermionicImaginaryTimeGrid*
=#

"""
    FermionicImaginaryTimeGrid(ntime::I64, β::T) where {T}

A constructor for the FermionicImaginaryTimeGrid struct, which is defined
in `src/types.jl`.

### Arguments
* ntime -> Number of time slices in imaginary axis.
* β     -> Inverse temperature.

### Returns
* grid -> A FermionicImaginaryTimeGrid struct.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function FermionicImaginaryTimeGrid(ntime::I64, β::T) where {T}
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    τ = collect(LinRange(0.0, β, ntime))
    return FermionicImaginaryTimeGrid(ntime, β, τ)
end

"""
    Base.length(fg::FermionicImaginaryTimeGrid)

Return number of grid points in a FermionicImaginaryTimeGrid struct.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.length(fg::FermionicImaginaryTimeGrid)
    fg.ntime
end

"""
    Base.iterate(fg::FermionicImaginaryTimeGrid)

Advance the iterator of a FermionicImaginaryTimeGrid struct to obtain
the next grid point.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.iterate(fg::FermionicImaginaryTimeGrid)
    iterate(fg.τ)
end

"""
    Base.iterate(fg::FermionicImaginaryTimeGrid, i::I64)

This is the key method that allows a FermionicImaginaryTimeGrid struct
to be iterated, yielding a sequences of grid points.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.iterate(fg::FermionicImaginaryTimeGrid, i::I64)
    iterate(fg.τ, i)
end

"""
    Base.eachindex(fg::FermionicImaginaryTimeGrid)

Create an iterable object for visiting each index of a
FermionicImaginaryTimeGrid struct.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.eachindex(fg::FermionicImaginaryTimeGrid)
    eachindex(fg.τ)
end

"""
    Base.firstindex(fg::FermionicImaginaryTimeGrid)

Return the first index of a FermionicImaginaryTimeGrid struct.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.firstindex(fg::FermionicImaginaryTimeGrid)
    firstindex(fg.τ)
end

"""
    Base.lastindex(fg::FermionicImaginaryTimeGrid)

Return the last index of a FermionicImaginaryTimeGrid struct.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.lastindex(fg::FermionicImaginaryTimeGrid)
    lastindex(fg.τ)
end

"""
    Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
FermionicImaginaryTimeGrid struct.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ntime
    return fg.τ[ind]
end

"""
    Base.getindex(fg::FermionicImaginaryTimeGrid, I::UnitRange{I64})

Return a subset of a FermionicImaginaryTimeGrid struct as specified by `I`.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function Base.getindex(fg::FermionicImaginaryTimeGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, fg.τ, I)
    lI = length(I)
    X = similar(fg.τ, lI)
    if lI > 0
        unsafe_copyto!(X, 1, fg.τ, first(I), lI)
    end
    return X
end

"""
    rebuild!(fg::FermionicImaginaryTimeGrid, ntime::I64, β::T) where {T}

Rebuild the FermionicImaginaryTimeGrid struct via new `ntime` and `β`
parameters.

### Arguments
* fg -> A FermionicImaginaryTimeGrid struct.
* ntime -> Number of time slices.
* β -> Inverse temperature.

### Returns
N/A

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function rebuild!(fg::FermionicImaginaryTimeGrid, ntime::I64, β::T) where {T}
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    fg.ntime = ntime
    fg.β = β
    fg.τ = collect(LinRange(0.0, fg.β, fg.ntime))
end

#=
### *Struct : FermionicFragmentTimeGrid*
=#

"""
    FermionicFragmentTimeGrid(β::T, τ::Vector{T}) where {T}

A constructor for the FermionicFragmentTimeGrid struct, which is defined
in `src/types.jl`.

### Arguments
* β -> Inverse temperature.
* τ -> Given imaginary time points.

### Returns
* grid -> A FermionicFragmentTimeGrid struct.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function FermionicFragmentTimeGrid(β::T, τ::Vector{T}) where {T}
    ntime = length(τ)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    @assert all(x -> (0.0 ≤ x ≤ β), τ)
    return FermionicFragmentTimeGrid(ntime, β, τ)
end

"""
    Base.length(fg::FermionicFragmentTimeGrid)

Return number of grid points in a FermionicFragmentTimeGrid struct.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.length(fg::FermionicFragmentTimeGrid)
    fg.ntime
end

"""
    Base.iterate(fg::FermionicFragmentTimeGrid)

Advance the iterator of a FermionicFragmentTimeGrid struct to obtain
the next grid point.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.iterate(fg::FermionicFragmentTimeGrid)
    iterate(fg.τ)
end

"""
    Base.iterate(fg::FermionicFragmentTimeGrid, i::I64)

This is the key method that allows a FermionicFragmentTimeGrid struct
to be iterated, yielding a sequences of grid points.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.iterate(fg::FermionicFragmentTimeGrid, i::I64)
    iterate(fg.τ, i)
end

"""
    Base.eachindex(fg::FermionicFragmentTimeGrid)

Create an iterable object for visiting each index of a
FermionicFragmentTimeGrid struct.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.eachindex(fg::FermionicFragmentTimeGrid)
    eachindex(fg.τ)
end

"""
    Base.firstindex(fg::FermionicFragmentTimeGrid)

Return the first index of a FermionicFragmentTimeGrid struct.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.firstindex(fg::FermionicFragmentTimeGrid)
    firstindex(fg.τ)
end

"""
    Base.lastindex(fg::FermionicFragmentTimeGrid)

Return the last index of a FermionicFragmentTimeGrid struct.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.lastindex(fg::FermionicFragmentTimeGrid)
    lastindex(fg.τ)
end

"""
    Base.getindex(fg::FermionicFragmentTimeGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
FermionicFragmentTimeGrid struct.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.getindex(fg::FermionicFragmentTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ntime
    return fg.τ[ind]
end

"""
    Base.getindex(fg::FermionicFragmentTimeGrid, I::UnitRange{I64})

Return a subset of a FermionicFragmentTimeGrid struct as specified by `I`.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function Base.getindex(fg::FermionicFragmentTimeGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, fg.τ, I)
    lI = length(I)
    X = similar(fg.τ, lI)
    if lI > 0
        unsafe_copyto!(X, 1, fg.τ, first(I), lI)
    end
    return X
end

"""
    rebuild!(fg::FermionicFragmentTimeGrid, ntime::I64, β::T) where {T}

Rebuild the FermionicFragmentTimeGrid struct via new `ntime` and `β`
parameters. Now its imaginary time points are continuous and complete.

### Arguments
* fg -> A FermionicFragmentTimeGrid struct.
* ntime -> Number of time slices.
* β -> Inverse temperature.

### Returns
N/A

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function rebuild!(fg::FermionicFragmentTimeGrid, ntime::I64, β::T) where {T}
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    fg.ntime = ntime
    fg.β = β
    fg.τ = collect(LinRange(0.0, fg.β, fg.ntime))
end

#=
### *Struct : FermionicMatsubaraGrid*
=#

"""
    FermionicMatsubaraGrid(nfreq::I64, β::T) where {T}

A constructor for the FermionicMatsubaraGrid struct, which is defined in
`src/types.jl`. The Matsubara grid is evaluated as ωₙ = (2n - 1) π / β.

### Arguments
* nfreq -> Number of Matsubara frequencies.
* β     -> Inverse temperature.

### Returns
* grid -> A FermionicMatsubaraGrid struct.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function FermionicMatsubaraGrid(nfreq::I64, β::T) where {T}
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    wmin = π / β
    wmax = (2 * nfreq - 1) * π / β
    ω = collect(LinRange(wmin, wmax, nfreq))
    return FermionicMatsubaraGrid(nfreq, β, ω)
end

"""
    Base.length(fg::FermionicMatsubaraGrid)

Return number of grid points in a FermionicMatsubaraGrid struct.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.length(fg::FermionicMatsubaraGrid)
    fg.nfreq
end

"""
    Base.iterate(fg::FermionicMatsubaraGrid)

Advance the iterator of a FermionicMatsubaraGrid struct to obtain
the next grid point.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.iterate(fg::FermionicMatsubaraGrid)
    iterate(fg.ω)
end

"""
    Base.iterate(fg::FermionicMatsubaraGrid, i::I64)

Create an iterable object for visiting each index of a
FermionicMatsubaraGrid struct.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.iterate(fg::FermionicMatsubaraGrid, i::I64)
    iterate(fg.ω, i)
end

"""
    Base.eachindex(fg::FermionicMatsubaraGrid)

Create an iterable object for visiting each index of a
FermionicMatsubaraGrid struct.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.eachindex(fg::FermionicMatsubaraGrid)
    eachindex(fg.ω)
end

"""
    Base.firstindex(fg::FermionicMatsubaraGrid)

Return the first index of a FermionicMatsubaraGrid struct.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.firstindex(fg::FermionicMatsubaraGrid)
    firstindex(fg.ω)
end

"""
    Base.lastindex(fg::FermionicMatsubaraGrid)

Return the last index of a FermionicMatsubaraGrid struct.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.lastindex(fg::FermionicMatsubaraGrid)
    lastindex(fg.ω)
end

"""
    Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
FermionicMatsubaraGrid struct.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.nfreq
    return fg.ω[ind]
end

"""
    Base.getindex(fg::FermionicMatsubaraGrid, I::UnitRange{I64})

Return a subset of a FermionicMatsubaraGrid struct as specified by `I`.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.getindex(fg::FermionicMatsubaraGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, fg.ω, I)
    lI = length(I)
    X = similar(fg.ω, lI)
    if lI > 0
        unsafe_copyto!(X, 1, fg.ω, first(I), lI)
    end
    return X
end

"""
    rebuild!(fg::FermionicMatsubaraGrid, nfreq::I64, β::T) where {T}

Rebuild the FermionicMatsubaraGrid struct via new `nfreq` and `β`
parameters.

### Arguments
* fg -> A FermionicMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.
* β -> Inverse temperature.

### Returns
N/A

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function rebuild!(fg::FermionicMatsubaraGrid, nfreq::I64, β::T) where {T}
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    fg.nfreq = nfreq
    fg.β = β
    resize!(fg.ω, nfreq)
    for n = 1:nfreq
        fg.ω[n] = (2 * n - 1) * π / fg.β
    end
end

"""
    Base.resize!(fg::FermionicMatsubaraGrid, nfreq::I64)

Reduce the size of the fermionic Matsubara grid. Note that `nfreq` should
be smaller than or equal to `fg.nfreq`. This function is called by the
NevanAC solver only.

### Arguments
* fg -> A FermionicMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.

### Returns
N/A

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.resize!(fg::FermionicMatsubaraGrid, nfreq::I64)
    @assert fg.nfreq ≥ nfreq
    fg.nfreq = nfreq
    resize!(fg.ω, nfreq)
end

"""
    Base.reverse!(fg::FermionicMatsubaraGrid)

Reverse the fermionic Matsubara grid. This function is called by the
`NevanAC` solver only.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function Base.reverse!(fg::FermionicMatsubaraGrid)
    reverse!(fg.ω)
end

#=
### *Struct : FermionicFragmentMatsubaraGrid*
=#

"""
    FermionicFragmentMatsubaraGrid(β::T, ω::Vector{T}) where {T}

A constructor for the FermionicFragmentMatsubaraGrid struct, which is
defined in `src/types.jl`. The Matsubara grid is from input.

### Arguments
* β -> Inverse temperature.
* ω -> Given Matsubara frequency points.

### Returns
* grid -> A FermionicFragmentMatsubaraGrid struct.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function FermionicFragmentMatsubaraGrid(β::T, ω::Vector{T}) where {T}
    nfreq = length(ω)
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    wmin = π / β
    wmax = (2 * nfreq - 1) * π / β
    @assert all(x -> (wmin ≤ x ≤ wmax * 2.0), ω)
    return FermionicFragmentMatsubaraGrid(nfreq, β, ω)
end

"""
    Base.length(fg::FermionicFragmentMatsubaraGrid)

Return number of grid points in a FermionicFragmentMatsubaraGrid struct.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.length(fg::FermionicFragmentMatsubaraGrid)
    fg.nfreq
end

"""
    Base.iterate(fg::FermionicFragmentMatsubaraGrid)

Advance the iterator of a FermionicFragmentMatsubaraGrid struct to obtain
the next grid point.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.iterate(fg::FermionicFragmentMatsubaraGrid)
    iterate(fg.ω)
end

"""
    Base.iterate(fg::FermionicFragmentMatsubaraGrid, i::I64)

Create an iterable object for visiting each index of a
FermionicFragmentMatsubaraGrid struct.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.iterate(fg::FermionicFragmentMatsubaraGrid, i::I64)
    iterate(fg.ω, i)
end

"""
    Base.eachindex(fg::FermionicFragmentMatsubaraGrid)

Create an iterable object for visiting each index of a
FermionicFragmentMatsubaraGrid struct.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.eachindex(fg::FermionicFragmentMatsubaraGrid)
    eachindex(fg.ω)
end

"""
    Base.firstindex(fg::FermionicFragmentMatsubaraGrid)

Return the first index of a FermionicFragmentMatsubaraGrid struct.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.firstindex(fg::FermionicFragmentMatsubaraGrid)
    firstindex(fg.ω)
end

"""
    Base.lastindex(fg::FermionicFragmentMatsubaraGrid)

Return the last index of a FermionicFragmentMatsubaraGrid struct.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.lastindex(fg::FermionicFragmentMatsubaraGrid)
    lastindex(fg.ω)
end

"""
    Base.getindex(fg::FermionicFragmentMatsubaraGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
FermionicFragmentMatsubaraGrid struct.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.getindex(fg::FermionicFragmentMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.nfreq
    return fg.ω[ind]
end

"""
    Base.getindex(fg::FermionicFragmentMatsubaraGrid, I::UnitRange{I64})

Return a subset of a FermionicFragmentMatsubaraGrid struct as specified by `I`.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.getindex(fg::FermionicFragmentMatsubaraGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, fg.ω, I)
    lI = length(I)
    X = similar(fg.ω, lI)
    if lI > 0
        unsafe_copyto!(X, 1, fg.ω, first(I), lI)
    end
    return X
end

"""
    rebuild!(fg::FermionicFragmentMatsubaraGrid, nfreq::I64, β::T) where {T}

Rebuild the FermionicFragmentMatsubaraGrid struct via new `nfreq` and `β`
parameters. Now its Matsubara frequency points are continuous and complete.

### Arguments
* fg -> A FermionicFragmentMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.
* β -> Inverse temperature.

### Returns
N/A

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function rebuild!(fg::FermionicFragmentMatsubaraGrid, nfreq::I64, β::T) where {T}
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    fg.nfreq = nfreq
    fg.β = β
    resize!(fg.ω, nfreq)
    for n = 1:nfreq
        fg.ω[n] = (2 * n - 1) * π / fg.β
    end
end

"""
    Base.resize!(fg::FermionicFragmentMatsubaraGrid, nfreq::I64)

Reduce the size of the fermionic fragment Matsubara grid. Note that
`nfreq` should be smaller than or equal to `fg.nfreq`. This function
is called by the NevanAC solver only.

### Arguments
* fg -> A FermionicFragmentMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.

### Returns
N/A

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.resize!(fg::FermionicFragmentMatsubaraGrid, nfreq::I64)
    @assert fg.nfreq ≥ nfreq
    fg.nfreq = nfreq
    resize!(fg.ω, nfreq)
end

"""
    Base.reverse!(fg::FermionicFragmentMatsubaraGrid)

Reverse the fermionic fragment Matsubara grid. This function is called
by the `NevanAC` solver only.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function Base.reverse!(fg::FermionicFragmentMatsubaraGrid)
    reverse!(fg.ω)
end

#=
### *Struct : BosonicImaginaryTimeGrid*
=#

"""
    BosonicImaginaryTimeGrid(ntime::I64, β::T)

A constructor for the BosonicImaginaryTimeGrid struct, which is defined
in `src/types.jl`.

### Arguments
* ntime -> Number of time slices in imaginary axis.
* β     -> Inverse temperature.

### Returns
* grid -> A BosonicImaginaryTimeGrid struct.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function BosonicImaginaryTimeGrid(ntime::I64, β::T) where {T}
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    τ = collect(LinRange(0.0, β, ntime))
    return BosonicImaginaryTimeGrid(ntime, β, τ)
end

"""
    Base.length(bg::BosonicImaginaryTimeGrid)

Return number of grid points in a BosonicImaginaryTimeGrid struct.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.length(bg::BosonicImaginaryTimeGrid)
    bg.ntime
end

"""
    Base.iterate(bg::BosonicImaginaryTimeGrid)

Advance the iterator of a BosonicImaginaryTimeGrid struct to obtain
the next grid point.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.iterate(bg::BosonicImaginaryTimeGrid)
    iterate(bg.τ)
end

"""
    Base.iterate(bg::BosonicImaginaryTimeGrid, i::I64)

Create an iterable object for visiting each index of a
BosonicImaginaryTimeGrid struct.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.iterate(bg::BosonicImaginaryTimeGrid, i::I64)
    iterate(bg.τ, i)
end

"""
    Base.eachindex(bg::BosonicImaginaryTimeGrid)

Create an iterable object for visiting each index of a
BosonicImaginaryTimeGrid struct.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.eachindex(bg::BosonicImaginaryTimeGrid)
    eachindex(bg.τ)
end

"""
    Base.firstindex(bg::BosonicImaginaryTimeGrid)

Return the first index of a BosonicImaginaryTimeGrid struct.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.firstindex(bg::BosonicImaginaryTimeGrid)
    firstindex(bg.τ)
end

"""
    Base.lastindex(bg::BosonicImaginaryTimeGrid)

Return the last index of a BosonicImaginaryTimeGrid struct.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.lastindex(bg::BosonicImaginaryTimeGrid)
    lastindex(bg.τ)
end

"""
    Base.getindex(bg::BosonicImaginaryTimeGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
BosonicImaginaryTimeGrid struct.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.getindex(bg::BosonicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.ntime
    return bg.τ[ind]
end

"""
    Base.getindex(bg::BosonicImaginaryTimeGrid, I::UnitRange{I64})

Return a subset of a BosonicImaginaryTimeGrid struct as specified by `I`.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function Base.getindex(bg::BosonicImaginaryTimeGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, bg.τ, I)
    lI = length(I)
    X = similar(bg.τ, lI)
    if lI > 0
        unsafe_copyto!(X, 1, bg.τ, first(I), lI)
    end
    return X
end

"""
    rebuild!(bg::BosonicImaginaryTimeGrid, ntime::I64, β::T) where {T}

Rebuild the BosonicImaginaryTimeGrid struct via new `ntime` and `β`
parameters.

### Arguments
* bg -> A BosonicImaginaryTimeGrid struct.
* ntime -> Number of time slices.
* β -> Inverse temperature.

### Returns
N/A

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function rebuild!(bg::BosonicImaginaryTimeGrid, ntime::I64, β::T) where {T}
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    bg.ntime = ntime
    bg.β = β
    bg.τ = collect(LinRange(0.0, bg.β, bg.ntime))
end

#=
### *Struct : BosonicFragmentTimeGrid*
=#

"""
    BosonicFragmentTimeGrid(β::T, τ::Vector{T}) where {T}

A constructor for the BosonicFragmentTimeGrid struct, which is defined
in `src/types.jl`.

### Arguments
* β -> Inverse temperature.
* τ -> Given imaginary time points.

### Returns
* grid -> A BosonicFragmentTimeGrid struct.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function BosonicFragmentTimeGrid(β::T, τ::Vector{T}) where {T}
    ntime = length(τ)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    @assert all(x -> (0.0 ≤ x ≤ β), τ)
    return BosonicFragmentTimeGrid(ntime, β, τ)
end

"""
    Base.length(bg::BosonicFragmentTimeGrid)

Return number of grid points in a BosonicFragmentTimeGrid struct.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.length(bg::BosonicFragmentTimeGrid)
    bg.ntime
end

"""
    Base.iterate(bg::BosonicFragmentTimeGrid)

Advance the iterator of a BosonicFragmentTimeGrid struct to obtain
the next grid point.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.iterate(bg::BosonicFragmentTimeGrid)
    iterate(bg.τ)
end

"""
    Base.iterate(bg::BosonicFragmentTimeGrid, i::I64)

Create an iterable object for visiting each index of a
BosonicFragmentTimeGrid struct.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.iterate(bg::BosonicFragmentTimeGrid, i::I64)
    iterate(bg.τ, i)
end

"""
    Base.eachindex(bg::BosonicFragmentTimeGrid)

Create an iterable object for visiting each index of a
BosonicFragmentTimeGrid struct.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.eachindex(bg::BosonicFragmentTimeGrid)
    eachindex(bg.τ)
end

"""
    Base.firstindex(bg::BosonicFragmentTimeGrid)

Return the first index of a BosonicFragmentTimeGrid struct.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.firstindex(bg::BosonicFragmentTimeGrid)
    firstindex(bg.τ)
end

"""
    Base.lastindex(bg::BosonicFragmentTimeGrid)

Return the last index of a BosonicFragmentTimeGrid struct.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.lastindex(bg::BosonicFragmentTimeGrid)
    lastindex(bg.τ)
end

"""
    Base.getindex(bg::BosonicFragmentTimeGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
BosonicFragmentTimeGrid struct.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.getindex(bg::BosonicFragmentTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.ntime
    return bg.τ[ind]
end

"""
    Base.getindex(bg::BosonicFragmentTimeGrid, I::UnitRange{I64})

Return a subset of a BosonicFragmentTimeGrid struct as specified by `I`.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function Base.getindex(bg::BosonicFragmentTimeGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, bg.τ, I)
    lI = length(I)
    X = similar(bg.τ, lI)
    if lI > 0
        unsafe_copyto!(X, 1, bg.τ, first(I), lI)
    end
    return X
end

"""
    rebuild!(bg::BosonicFragmentTimeGrid, ntime::I64, β::T) where {T}

Rebuild the BosonicFragmentTimeGrid struct via new `ntime` and `β`
parameters. Now its imaginary time points are continuous and complete.

### Arguments
* bg -> A BosonicFragmentTimeGrid struct.
* ntime -> Number of time slices.
* β -> Inverse temperature.

### Returns
N/A

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function rebuild!(bg::BosonicFragmentTimeGrid, ntime::I64, β::T) where {T}
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    bg.ntime = ntime
    bg.β = β
    bg.τ = collect(LinRange(0.0, bg.β, bg.ntime))
end

#=
### *Struct : BosonicMatsubaraGrid*
=#

"""
    BosonicMatsubaraGrid(nfreq::I64, β::T) where {T}

A constructor for the BosonicMatsubaraGrid struct, which is defined in
`src/types.jl`. The Matsubara grid is evaluated as ωₙ = (2n - 2) π / β.

### Arguments
* nfreq -> Number of Matsubara frequencies.
* β     -> Inverse temperature.

### Returns
* grid -> A BosonicMatsubaraGrid struct.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function BosonicMatsubaraGrid(nfreq::I64, β::T) where {T}
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    wmin = 0.0
    wmax = (2 * nfreq - 2) * π / β
    ω = collect(LinRange(wmin, wmax, nfreq))
    return BosonicMatsubaraGrid(nfreq, β, ω)
end

"""
    Base.length(bg::BosonicMatsubaraGrid)

Return number of grid points in a BosonicMatsubaraGrid struct.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.length(bg::BosonicMatsubaraGrid)
    bg.nfreq
end

"""
    Base.iterate(bg::BosonicMatsubaraGrid)

Advance the iterator of a BosonicMatsubaraGrid struct to obtain
the next grid point.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.iterate(bg::BosonicMatsubaraGrid)
    iterate(bg.ω)
end

"""
    Base.iterate(bg::BosonicMatsubaraGrid, i::I64)

Create an iterable object for visiting each index of a
BosonicMatsubaraGrid struct.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.iterate(bg::BosonicMatsubaraGrid, i::I64)
    iterate(bg.ω, i)
end

"""
    Base.eachindex(bg::BosonicMatsubaraGrid)

Create an iterable object for visiting each index of a
BosonicMatsubaraGrid struct.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.eachindex(bg::BosonicMatsubaraGrid)
    eachindex(bg.ω)
end

"""
    Base.firstindex(bg::BosonicMatsubaraGrid)

Return the first index of a BosonicMatsubaraGrid struct.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.firstindex(bg::BosonicMatsubaraGrid)
    firstindex(bg.ω)
end

"""
    Base.lastindex(bg::BosonicMatsubaraGrid)

Return the last index of a BosonicMatsubaraGrid struct.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.lastindex(bg::BosonicMatsubaraGrid)
    lastindex(bg.ω)
end

"""
    Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
BosonicMatsubaraGrid struct.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.nfreq
    return bg.ω[ind]
end

"""
    Base.getindex(bg::BosonicMatsubaraGrid, I::UnitRange{I64})

Return a subset of a BosonicMatsubaraGrid struct as specified by `I`.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.getindex(bg::BosonicMatsubaraGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, bg.ω, I)
    lI = length(I)
    X = similar(bg.ω, lI)
    if lI > 0
        unsafe_copyto!(X, 1, bg.ω, first(I), lI)
    end
    return X
end

"""
    rebuild!(bg::BosonicMatsubaraGrid, nfreq::I64, β::T) where {T}

Rebuild the BosonicMatsubaraGrid struct via new `nfreq` and `β`
parameters.

### Arguments
* bg -> A BosonicMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.
* β -> Inverse temperature.

### Returns
N/A

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function rebuild!(bg::BosonicMatsubaraGrid, nfreq::I64, β::T) where {T}
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    bg.nfreq = nfreq
    bg.β = β
    resize!(bg.ω, nfreq)
    for n = 1:nfreq
        bg.ω[n] = (2 * n - 2) * π / bg.β
    end
end

"""
    Base.resize!(bg::BosonicMatsubaraGrid, nfreq::I64)

Reduce the size of the bosonic Matsubara grid. Note that `nfreq` should
be smaller than or equal to `bg.nfreq`. This function is called by the
NevanAC solver only.

### Arguments
* bg -> A BosonicMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.

### Returns
N/A

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.resize!(bg::BosonicMatsubaraGrid, nfreq::I64)
    @assert bg.nfreq ≥ nfreq
    bg.nfreq = nfreq
    resize!(bg.ω, nfreq)
end

"""
    Base.reverse!(bg::BosonicMatsubaraGrid)

Reverse the bosonic Matsubara grid. This function is called by the
`NevanAC` solver only.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function Base.reverse!(bg::BosonicMatsubaraGrid)
    reverse!(bg.ω)
end

#=
### *Struct : BosonicFragmentMatsubaraGrid*
=#

"""
    BosonicFragmentMatsubaraGrid(β::T, ω::Vector{T}) where {T}

A constructor for the BosonicFragmentMatsubaraGrid struct, which is
defined in `src/types.jl`. The Matsubara grid is from input.

### Arguments
* β -> Inverse temperature.
* ω -> Given Matsubara frequency points.

### Returns
* grid -> A BosonicFragmentMatsubaraGrid struct.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function BosonicFragmentMatsubaraGrid(β::T, ω::Vector{T}) where {T}
    nfreq = length(ω)
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    wmin = 0.0
    wmax = (2 * nfreq - 2) * π / β
    @assert all(x -> (wmin ≤ x ≤ wmax * 2.0), ω)
    return BosonicFragmentMatsubaraGrid(nfreq, β, ω)
end

"""
    Base.length(bg::BosonicFragmentMatsubaraGrid)

Return number of grid points in a BosonicFragmentMatsubaraGrid struct.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.length(bg::BosonicFragmentMatsubaraGrid)
    bg.nfreq
end

"""
    Base.iterate(bg::BosonicFragmentMatsubaraGrid)

Advance the iterator of a BosonicFragmentMatsubaraGrid struct to obtain
the next grid point.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.iterate(bg::BosonicFragmentMatsubaraGrid)
    iterate(bg.ω)
end

"""
    Base.iterate(bg::BosonicFragmentMatsubaraGrid, i::I64)

Create an iterable object for visiting each index of a
BosonicFragmentMatsubaraGrid struct.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.iterate(bg::BosonicFragmentMatsubaraGrid, i::I64)
    iterate(bg.ω, i)
end

"""
    Base.eachindex(bg::BosonicFragmentMatsubaraGrid)

Create an iterable object for visiting each index of a
BosonicFragmentMatsubaraGrid struct.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.eachindex(bg::BosonicFragmentMatsubaraGrid)
    eachindex(bg.ω)
end

"""
    Base.firstindex(bg::BosonicFragmentMatsubaraGrid)

Return the first index of a BosonicFragmentMatsubaraGrid struct.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.firstindex(bg::BosonicFragmentMatsubaraGrid)
    firstindex(bg.ω)
end

"""
    Base.lastindex(bg::BosonicFragmentMatsubaraGrid)

Return the last index of a BosonicFragmentMatsubaraGrid struct.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.lastindex(bg::BosonicFragmentMatsubaraGrid)
    lastindex(bg.ω)
end

"""
    Base.getindex(bg::BosonicFragmentMatsubaraGrid, ind::I64)

Retrieve the value(s) stored at the given key or index within a
BosonicFragmentMatsubaraGrid struct.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.getindex(bg::BosonicFragmentMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.nfreq
    return bg.ω[ind]
end

"""
    Base.getindex(bg::BosonicFragmentMatsubaraGrid, I::UnitRange{I64})

Return a subset of a BosonicFragmentMatsubaraGrid struct as specified by `I`.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.getindex(bg::BosonicFragmentMatsubaraGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, bg.ω, I)
    lI = length(I)
    X = similar(bg.ω, lI)
    if lI > 0
        unsafe_copyto!(X, 1, bg.ω, first(I), lI)
    end
    return X
end

"""
    rebuild!(bg::BosonicFragmentMatsubaraGrid, nfreq::I64, β::T) where {T}

Rebuild the BosonicFragmentMatsubaraGrid struct via new `nfreq` and `β`
parameters. Now its Matsubara frequency points are continuous and complete.

### Arguments
* bg -> A BosonicFragmentMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.
* β -> Inverse temperature.

### Returns
N/A

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function rebuild!(bg::BosonicFragmentMatsubaraGrid, nfreq::I64, β::T) where {T}
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    bg.nfreq = nfreq
    bg.β = β
    resize!(bg.ω, nfreq)
    for n = 1:nfreq
        bg.ω[n] = (2 * n - 2) * π / bg.β
    end
end

"""
    Base.resize!(bg::BosonicFragmentMatsubaraGrid, nfreq::I64)

Reduce the size of the bosonic fragment Matsubara grid. Note that
`nfreq` should be smaller than or equal to `bg.nfreq`. This function
is called by the NevanAC solver only.

### Arguments
* bg -> A BosonicFragmentMatsubaraGrid struct.
* nfreq -> Number of Matsubara frequencies.

### Returns
N/A

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.resize!(bg::BosonicFragmentMatsubaraGrid, nfreq::I64)
    @assert bg.nfreq ≥ nfreq
    bg.nfreq = nfreq
    resize!(bg.ω, nfreq)
end

"""
    Base.reverse!(bg::BosonicFragmentMatsubaraGrid)

Reverse the bosonic fragment Matsubara grid. This function is called
by the `NevanAC` solver only.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function Base.reverse!(bg::BosonicFragmentMatsubaraGrid)
    reverse!(bg.ω)
end
