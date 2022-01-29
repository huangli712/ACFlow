#
# Project : Gardenia
# Source  : grid.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/01/29
#

function Base.getindex(fg::FermionicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.ntime
    return fg.τ[ind]
end

function Base.getindex(fg::FermionicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ fg.nfreq
    return fg.ω[ind]
end

function Base.getindex(bg::BosonicImaginaryTimeGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.ntime
    return bg.τ[ind]
end

function Base.getindex(bg::BosonicMatsubaraGrid, ind::I64)
    @assert 1 ≤ ind ≤ bg.nfreq
    return bg.ω[ind]
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
        @assert β == get_c("beta")
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

function rebuild_grid(fg::FermionicImaginaryTimeGrid)
    fg.τ = collect(LinRange(0.0, fg.β, fg.ntime))
end

function rebuild_grid(fg::FermionicMatsubaraGrid)
    for n in eachindex(fg.ω)
        fg.ω[n] = (2 * n - 1) * π / fg.β
    end
end

function rebuild_grid(bg::BosonicImaginaryTimeGrid)
    bg.τ = collect(LinRange(0.0, bg.β, bg.ntime))
end

function rebuild_grid(bg::BosonicMatsubaraGrid)
    for n in eachindex(bg.ω)
        bg.ω[n] = (2 * n - 2) * π / bg.β
    end
end