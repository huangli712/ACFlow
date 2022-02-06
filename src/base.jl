#
# Project : Gardenia
# Source  : flow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/06
#

function solve(rd::RawData)
    solver = get_c("solver")

    @cswitch solver begin
        @case "MaxEnt"
            solve(MaxEntSolver(), rd)
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

function make_model(am::AbstractMesh)
    mtype = get_c("mtype")

    @cswitch mtype begin
        @case "flat"
            return build_flat_model(am)
            break

        @case "gauss"
            return build_gaussian_model(am)
            break

        @case "func"
            sorry()
            break

        @default
            sorry()
            break
    end
end

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
