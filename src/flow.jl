#
# Project : Gardenia
# Source  : flow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
#

abstract type AbstractData end

struct RawData{T} <: AbstractData
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

function make_data(rd::RawData)
    return make_data(rd.value, rd.error)
end

function make_data(val::Vector{T}, err::Vector{T}) where {T}
    grid = get_c("grid")
    kernel = get_c("kernel")

    if grid == "matsubara" && kernel == "fermionic"
        value = vcat(real(val), imag(val))
        error = vcat(real(err), imag(err))
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    end

    if grid == "matsubara" && kernel == "bosonic"
        value = real(val)
        error = real(err)
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    end

    if grid == "time" && kernel == "fermionic"
        value = real(val)
        error = real(err)
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    end

    if grid == "time" && kernel == "bosonic"
        value = real(val)
        error = real(err)
        covar = error .^ 2.0
        _data = GreenData(value, error, covar)
        return _data
    end
end

function read_data()
    finput = get_c("finput")
    ngrid = get_c("ngrid")
    grid = get_c("grid")

    @cswitch grid begin
        @case "matsubara"
            return read_freq_data(finput, ngrid)
            break

        @case "time"
            return read_time_data(finput, ngrid)
            break

        @default
            sorry()
            break
    end
end
