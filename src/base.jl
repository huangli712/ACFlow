#
# Project : Gardenia
# Source  : flow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/12
#

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T})

Solve the analytical continuation problem. The arguments `grid`, `Gval`,
and `Gerr` are the grid, value, and error bar, respectively.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}, Gerr::Vector{T}) where {T}
    solve(RawData(grid, Gval, Gerr))
end

"""
    solve(grid::Vector{F64}, Gval::Vector{T}, err::T)

Solve the analytical continuation problem. The arguments `grid`, `Gval`,
and `err` are the grid, value, and error bar, respectively.
"""
function solve(grid::Vector, Gval::Vector{T}, err::T) where {T}
    Gerr = similar(Gval)
    fill!(Gerr, err)
    solve(RawData(grid, Gval, Gerr))
end

"""
    solve(grid::Vector{F64}, Gval::Vector{T})

Solve the analytical continuation problem. The arguments `grid` and `Gval`
are the grid and value, respectively. The error bar is set to 1.0e-4.
"""
function solve(grid::Vector{F64}, Gval::Vector{T}) where {T}
    Gerr = similar(Gval)
    fill!(Gerr, 1.0e-4)
    solve(RawData(grid, Gval, Gerr))
end

"""
    solve(rd::RawData)

See also: [`RawData`](@ref).
"""
function solve(rd::RawData)
    solver = get_c("solver")

    @cswitch solver begin
        @case "MaxEnt"
            solve(MaxEntSolver(), rd)
            break

        @case "StochAC"
            solve(StochACSolver(), rd)
            break

        @case "StochOM"
            solve(StochOMSolver(), rd)
            break

        @default
            sorry()
            break
    end
end

"""
    reprod
"""
function reprod(kernel::Matrix{F64}, am::AbstractMesh, A::Vector{F64})
    ndim, nmesh = size(kernel)
    @assert nmesh == length(am) == length(A)

    Ac = reshape(A, (1, nmesh))
    KA = kernel .* Ac

    G = zeros(F64, ndim)
    for i = 1:ndim
        G[i] = trapz(am, KA[i,:])
    end

    return G
end

"""
    setup_param
"""
function setup_param()
end

"""
    read_param
"""
function read_param()
    cfg = inp_toml(query_args(), true)
    fil_dict(cfg)
    chk_dict()
end

"""
    read_data
"""
function read_data()
    finput = get_c("finput")
    ngrid = get_c("ngrid")

    @cswitch get_c("grid") begin
        @case "ftime"
            return read_real_data(finput, ngrid)
            break

        @case "btime"
            return read_real_data(finput, ngrid)
            break

        @case "ffreq"
            return read_complex_data(finput, ngrid)
            break

        @case "bfreq"
            return read_complex_data(finput, ngrid, true)
            break

        @default
            sorry()
            break
    end
end

"""
    make_data
"""
function make_data(rd::RawData)
    grid = get_c("grid")
    val = rd.value
    err = rd.error

    if grid == "ffreq"
        value = vcat(real(val), imag(val))
        error = vcat(real(err), imag(err))
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    else
        value = real(val)
        error = real(err)
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    end
end

"""
    make_grid
"""
function make_grid(rd::RawData)
    grid = get_c("grid")
    ngrid = get_c("ngrid")
    v = rd._grid
    @assert ngrid == length(v)

    _grid = nothing
    @cswitch grid begin
        @case "ftime"
            β = v[end]
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = FermionicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "btime"
            β = v[end]
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = BosonicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "ffreq"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = FermionicMatsubaraGrid(ngrid, β, v)
            break

        @case "bfreq"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = BosonicMatsubaraGrid(ngrid, β, v)
            break

        @default
            sorry()
            break
    end

    return _grid
end

"""
    make_mesh
"""
function make_mesh()
    nmesh = get_c("nmesh")
    mesh = get_c("mesh")
    wmax = get_c("wmax")
    wmin = get_c("wmin")

    if mesh == "linear"
        return LinearMesh(nmesh, wmin, wmax)
    else
        return TangentMesh(nmesh, wmin, wmax)
    end
end

"""
    make_model
"""
function make_model(am::AbstractMesh)
    mtype = get_c("mtype")

    @cswitch mtype begin
        @case "flat"
            return build_flat_model(am)
            break

        @case "gauss"
            return build_gaussian_model(am)
            break

        @case "file"
            return build_file_model(am)
            break

        @case "func"
            sorry()
            break

        @default
            sorry()
            break
    end
end

"""
    make_kernel
"""
function make_kernel(am::AbstractMesh, ag::AbstractGrid)
    ktype = get_c("ktype")
    grid = get_c("grid")

    if ktype == "fermi" || ktype == "boson"
        return build_kernel(am, ag)
    else
        @assert grid in ("btime", "bfreq")
        return build_kernel_symm(am, ag)
    end
end
