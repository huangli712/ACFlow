#
# Project : Gardenia
# Source  : grid.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
#

abstract type AbstractGrid end

mutable struct FermionicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

function FermionicImaginaryTimeGrid(ntime::I64, β::F64)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    τ = collect(LinRange(0.0, β, ntime))
    return FermionicImaginaryTimeGrid(ntime, β, τ)
end

function Base.eachindex(fg::FermionicImaginaryTimeGrid)
    eachindex(fg.τ)
end

function Base.firstindex(fg::FermionicImaginaryTimeGrid)
    firstindex(fg.τ)
end

function Base.lastindex(fg::FermionicImaginaryTimeGrid)
    lastindex(fg.τ)
end

function Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ntime
    return fg.τ[ind]
end

function Base.getindex(fg::FermionicImaginaryTimeGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, fg.τ, I)
    lI = length(I)
    X = similar(fg.τ, lI)
    if lI > 0
        unsafe_copyto!(X, 1, fg.τ, first(I), lI)
    end
    return X
end

function rebuild_grid(fg::FermionicImaginaryTimeGrid, ntime::I64, β::F64)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    fg.ntime = ntime
    fg.β = β
    fg.τ = collect(LinRange(0.0, fg.β, fg.ntime))
end

mutable struct FermionicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

function FermionicMatsubaraGrid(nfreq::I64, β::F64)
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    wmin = π / β
    wmax = (2 * nfreq - 1) * π / β
    ω = collect(LinRange(wmin, wmax, nfreq))
    return FermionicMatsubaraGrid(nfreq, β, ω)
end

function Base.eachindex(fg::FermionicMatsubaraGrid)
    eachindex(fg.ω)
end

function Base.firstindex(fg::FermionicMatsubaraGrid)
    firstindex(fg.ω)
end

function Base.lastindex(fg::FermionicMatsubaraGrid)
    lastindex(fg.ω)
end

function Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.nfreq
    return fg.ω[ind]
end

function Base.getindex(fg::FermionicMatsubaraGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, fg.ω, I)
    lI = length(I)
    X = similar(fg.ω, lI)
    if lI > 0
        unsafe_copyto!(X, 1, fg.ω, first(I), lI)
    end
    return X
end

function rebuild_grid(fg::FermionicMatsubaraGrid, nfreq::I64, β::F64)
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    fg.nfreq = nfreq
    fg.β = β
    resize!(fg.ω, nfreq)
    for n = 1:nfreq
        fg.ω[n] = (2 * n - 1) * π / fg.β
    end
end

mutable struct BosonicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

function BosonicImaginaryTimeGrid(ntime::I64, β::F64)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    τ = collect(LinRange(0.0, β, ntime))
    return BosonicImaginaryTimeGrid(ntime, β, τ)
end

function Base.eachindex(bg::BosonicImaginaryTimeGrid)
    eachindex(bg.τ)
end

function Base.firstindex(bg::BosonicImaginaryTimeGrid)
    firstindex(bg.τ)
end

function Base.lastindex(bg::BosonicImaginaryTimeGrid)
    lastindex(bg.τ)
end

function Base.getindex(bg::BosonicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.ntime
    return bg.τ[ind]
end

function Base.getindex(bg::BosonicImaginaryTimeGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, bg.τ, I)
    lI = length(I)
    X = similar(bg.τ, lI)
    if lI > 0
        unsafe_copyto!(X, 1, bg.τ, first(I), lI)
    end
    return X
end

function rebuild_grid(bg::BosonicImaginaryTimeGrid)
    bg.τ = collect(LinRange(0.0, bg.β, bg.ntime))
end

mutable struct BosonicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

function BosonicMatsubaraGrid(nfreq::I64, β::F64)
    @assert nfreq ≥ 1
    @assert β ≥ 0.0
    wmin = 0.0
    wmax = (2 * nfreq - 2) * π / β
    ω = collect(LinRange(wmin, wmax, nfreq))
    return BosonicMatsubaraGrid(nfreq, β, ω)
end

function Base.eachindex(bg::BosonicMatsubaraGrid)
    eachindex(bg.ω)
end

function Base.firstindex(bg::BosonicMatsubaraGrid)
    firstindex(bg.ω)
end

function Base.lastindex(bg::BosonicMatsubaraGrid)
    lastindex(bg.ω)
end

function Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.nfreq
    return bg.ω[ind]
end

function Base.getindex(bg::BosonicMatsubaraGrid, I::UnitRange{I64})
    @assert checkbounds(Bool, bg.ω, I)
    lI = length(I)
    X = similar(bg.ω, lI)
    if lI > 0
        unsafe_copyto!(X, 1, bg.ω, first(I), lI)
    end
    return X
end

function rebuild_grid(bg::BosonicMatsubaraGrid)
    for n in eachindex(bg)
        bg.ω[n] = (2 * n - 2) * π / bg.β
    end
end

function make_grid(rd::RawData)
    return make_grid(rd.mesh)
end

function make_grid(v::Vector{F64})
    grid = get_c("grid")
    ngrid = get_c("ngrid")
    kernel = get_c("kernel")
    @assert ngrid == length(v)

    if grid == "matsubara"
        β = 2.0 * π / (v[2] - v[1])
        @assert abs(β - get_c("beta")) ≤ 1e-10
        if kernel == "fermionic"
            _grid = FermionicMatsubaraGrid(ngrid, β, v)
        else
            _grid = BosonicMatsubaraGrid(ngrid, β, v)
        end
        return _grid
    else
        β = v[end]
        @assert β == get_c("beta")
        if kernel == "fermionic"
            _grid = FermionicImaginaryTimeGrid(ngrid, β, v)
        else
            _grid = BosonicImaginaryTimeGrid(ngrid, β, v)
        end
        return _grid
    end
end