#
# Project : Gardenia
# Source  : grid.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/10/13
#

#=
### *Struct : FermionicImaginaryTimeGrid*
=#

"""
    FermionicImaginaryTimeGrid(ntime::I64, β::F64)

A constructor for the FermionicImaginaryTimeGrid struct, which is defined
in `src/types.jl`.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function FermionicImaginaryTimeGrid(ntime::I64, β::F64)
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
    rebuild(fg::FermionicImaginaryTimeGrid, ntime::I64, β::F64)

Rebuild the FermionicImaginaryTimeGrid struct via new `ntime` and `β`
parameters.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function rebuild(fg::FermionicImaginaryTimeGrid, ntime::I64, β::F64)
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
    FermionicMatsubaraGrid(nfreq::I64, β::F64)

A constructor for the FermionicMatsubaraGrid struct, which is defined in
`src/types.jl`. The Matsubara grid is evaluated as ωₙ = (2n - 1) π / β.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function FermionicMatsubaraGrid(nfreq::I64, β::F64)
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
    rebuild(fg::FermionicMatsubaraGrid, nfreq::I64, β::F64)

Rebuild the FermionicMatsubaraGrid struct via new `nfreq` and `β`
parameters.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function rebuild(fg::FermionicMatsubaraGrid, nfreq::I64, β::F64)
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    fg.nfreq = nfreq
    fg.β = β
    resize!(fg.ω, nfreq)
    for n = 1:nfreq
        fg.ω[n] = (2 * n - 1) * π / fg.β
    end
end

#=
### *Struct : BosonicImaginaryTimeGrid*
=#

"""
    BosonicImaginaryTimeGrid(ntime::I64, β::F64)

A constructor for the BosonicImaginaryTimeGrid struct, which is defined
in `src/types.jl`.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function BosonicImaginaryTimeGrid(ntime::I64, β::F64)
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
    rebuild(bg::BosonicImaginaryTimeGrid, ntime::I64, β::F64)

Rebuild the BosonicImaginaryTimeGrid struct via new `ntime` and `β`
parameters.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function rebuild(bg::BosonicImaginaryTimeGrid, ntime::I64, β::F64)
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
    BosonicMatsubaraGrid(nfreq::I64, β::F64)

A constructor for the BosonicMatsubaraGrid struct, which is defined in
`src/types.jl`. The Matsubara grid is evaluated as ωₙ = (2n - 2) π / β.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function BosonicMatsubaraGrid(nfreq::I64, β::F64)
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
    rebuild(bg::BosonicMatsubaraGrid, nfreq::I64, β::F64)

Rebuild the BosonicMatsubaraGrid struct via new `nfreq` and `β`
parameters.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function rebuild(bg::BosonicMatsubaraGrid, nfreq::I64, β::F64)
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    bg.nfreq = nfreq
    bg.β = β
    resize!(bg.ω, nfreq)
    for n = 1:nfreq
        bg.ω[n] = (2 * n - 2) * π / bg.β
    end
end
