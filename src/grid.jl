#
# Project : Gardenia
# Source  : grid.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/06
#

"""
    FermionicImaginaryTimeGrid
"""
function FermionicImaginaryTimeGrid(ntime::I64, β::F64)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    τ = collect(LinRange(0.0, β, ntime))
    return FermionicImaginaryTimeGrid(ntime, β, τ)
end

"""
    Base.length
"""
function Base.length(fg::FermionicImaginaryTimeGrid)
    fg.ntime
end

"""
    Base.iterate
"""
function Base.iterate(fg::FermionicImaginaryTimeGrid)
    iterate(fg.τ)
end

"""
    Base.iterate
"""
function Base.iterate(fg::FermionicImaginaryTimeGrid, i::I64)
    iterate(fg.τ, i)
end

"""
    Base.eachindex
"""
function Base.eachindex(fg::FermionicImaginaryTimeGrid)
    eachindex(fg.τ)
end

"""
    Base.firstindex
"""
function Base.firstindex(fg::FermionicImaginaryTimeGrid)
    firstindex(fg.τ)
end

"""
    Base.lastindex
"""
function Base.lastindex(fg::FermionicImaginaryTimeGrid)
    lastindex(fg.τ)
end

"""
    Base.getindex
"""
function Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ntime
    return fg.τ[ind]
end

"""
    Base.getindex
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
    rebuild
"""
function rebuild(fg::FermionicImaginaryTimeGrid, ntime::I64, β::F64)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    fg.ntime = ntime
    fg.β = β
    fg.τ = collect(LinRange(0.0, fg.β, fg.ntime))
end

"""
    FermionicMatsubaraGrid
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
    Base.length
"""
function Base.length(fg::FermionicMatsubaraGrid)
    fg.nfreq
end

"""
    Base.iterate
"""
function Base.iterate(fg::FermionicMatsubaraGrid)
    iterate(fg.ω)
end

"""
    Base.iterate
"""
function Base.iterate(fg::FermionicMatsubaraGrid, i::I64)
    iterate(fg.ω, i)
end

"""
    Base.eachindex
"""
function Base.eachindex(fg::FermionicMatsubaraGrid)
    eachindex(fg.ω)
end

"""
    Base.firstindex
"""
function Base.firstindex(fg::FermionicMatsubaraGrid)
    firstindex(fg.ω)
end

"""
    Base.lastindex
"""
function Base.lastindex(fg::FermionicMatsubaraGrid)
    lastindex(fg.ω)
end

"""
    Base.getindex
"""
function Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.nfreq
    return fg.ω[ind]
end

"""
    Base.getindex
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
    rebuild
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

"""
    BosonicImaginaryTimeGrid
"""
function BosonicImaginaryTimeGrid(ntime::I64, β::F64)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    τ = collect(LinRange(0.0, β, ntime))
    return BosonicImaginaryTimeGrid(ntime, β, τ)
end

"""
    Base.length
"""
function Base.length(bg::BosonicImaginaryTimeGrid)
    bg.ntime
end

"""
    Base.iterate
"""
function Base.iterate(bg::BosonicImaginaryTimeGrid)
    iterate(bg.τ)
end

"""
    Base.iterate
"""
function Base.iterate(bg::BosonicImaginaryTimeGrid, i::I64)
    iterate(bg.τ, i)
end

"""
    Base.eachindex
"""
function Base.eachindex(bg::BosonicImaginaryTimeGrid)
    eachindex(bg.τ)
end

"""
    Base.firstindex
"""
function Base.firstindex(bg::BosonicImaginaryTimeGrid)
    firstindex(bg.τ)
end

"""
    Base.lastindex
"""
function Base.lastindex(bg::BosonicImaginaryTimeGrid)
    lastindex(bg.τ)
end

"""
    Base.getindex
"""
function Base.getindex(bg::BosonicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.ntime
    return bg.τ[ind]
end

"""
    Base.getindex
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
    rebuild
"""
function rebuild(bg::BosonicImaginaryTimeGrid, ntime::I64, β::F64)
    @assert ntime ≥ 1
    @assert β ≥ 0.0
    bg.ntime = ntime
    bg.β = β
    bg.τ = collect(LinRange(0.0, bg.β, bg.ntime))
end

"""
    BosonicMatsubaraGrid
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
    Base.length
"""
function Base.length(bg::BosonicMatsubaraGrid)
    bg.nfreq
end

"""
    Base.iterate
"""
function Base.iterate(bg::BosonicMatsubaraGrid)
    iterate(bg.ω)
end

"""
    Base.iterate
"""
function Base.iterate(bg::BosonicMatsubaraGrid, i::I64)
    iterate(bg.ω, i)
end

"""
    Base.eachindex
"""
function Base.eachindex(bg::BosonicMatsubaraGrid)
    eachindex(bg.ω)
end

"""
    Base.firstindex
"""
function Base.firstindex(bg::BosonicMatsubaraGrid)
    firstindex(bg.ω)
end

"""
    Base.lastindex
"""
function Base.lastindex(bg::BosonicMatsubaraGrid)
    lastindex(bg.ω)
end

"""
    Base.getindex
"""
function Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.nfreq
    return bg.ω[ind]
end

"""
    Base.getindex
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
    rebuild
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
