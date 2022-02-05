#
# Project : Gardenia
# Source  : flow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/05
#

abstract type AbstractData end

mutable struct RawData{T} <: AbstractData
    _grid :: Vector{F64}
    value :: Vector{T}
    error :: Vector{T}
end

mutable struct GreenData <: AbstractData
    value :: Vector{F64}
    error :: Vector{F64}
    covar :: Vector{F64}
end

function solve(rd::RawData)
    solver = get_c("solver")
    
    @cswitch solver begin
        @case "MaxEnt"
            MaxEnt.solve(rd)
            break

        @case "StochOM"
            break

        @case "StochAC"
            StochAC.solve(rd)
            break

        @default
            sorry()
            break
    end
end

function setup_param()
end

function read_param()
    cfg = inp_toml(query_args(), true)
    fil_dict(cfg)
    chk_dict()
end

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

function make_grid(rd::RawData)
    grid = get_c("grid")
    ngrid = get_c("ngrid")
    v = rd._grid
    @assert ngrid == length(v)

    @cswitch grid begin
        @case "ftime"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = FermionicMatsubaraGrid(ngrid, β, v)
            break

        @case "btime"
            β = 2.0 * π / (v[2] - v[1])
            @assert abs(β - get_c("beta")) ≤ 1e-10
            _grid = BosonicMatsubaraGrid(ngrid, β, v)
            break

        @case "ffreq"
            β = v[end]
            @assert β == get_c("beta")
            _grid = FermionicImaginaryTimeGrid(ngrid, β, v)
            break

        @case "bfreq"
            β = v[end]
            @assert β == get_c("beta")
            _grid = BosonicImaginaryTimeGrid(ngrid, β, v)
            break

        @default
            sorry()
            break
    end

    return _grid
end

function make_mesh()
    nmesh = get_c("nmesh")
    mesh = get_c("mesh")
    wmax = get_c("wmax")
    wmin = get_c("wmin")

    if mesh == "uniform"
        return LinearMesh(nmesh, wmin, wmax)
    else
        return TangentMesh(nmesh, wmin, wmax)
    end
end

function make_model(m::AbstractMesh)
    model = get_c("model")
    if model == "flat"
        return make_flat_model(m)
    elseif model == "gauss"
        return make_gaussian_model(m)
    elseif model == "file"
        return make_file_model(m)
    end
end
