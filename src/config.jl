#
# Project : Gardenia
# Source  : config.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/05
#

#=
### *Customized Types*
=#

"Customized types. It is used to define the following dicts."
const DType = Any

"Customized types. It is used to define the following dicts."
const ADT = Array{DType,1}

const PCOMM    = Dict{String,ADT}(
    "finput"  => [missing, 1, :String, "Filename for input correlation function"],
    "solver"  => [missing, 1, :String, "Solve for the analytical continuation problem"],
    "ktype"   => [missing, 1, :String, "Type of kernel function"],
    "mtype"   => [missing, 1, :String, "Type of default model"],
    "grid"    => [missing, 1, :String, "Grid for input correlation function"],
    "mesh"    => [missing, 1, :String, "Mesh for output spectrum"],
    "ngrid"   => [missing, 1, :I64   , "Number of grid points"],
    "nmesh"   => [missing, 1, :I64   , "Number of mesh points"],
    "wmax"    => [missing, 1, :F64   , "Maximum value of mesh"],
    "wmin"    => [missing, 1, :F64   , "Minimum value of mesh"],
    "beta"    => [missing, 1, :F64   , "Inverse temperature"],
    "offdiag" => [missing, 1, :Bool  , "Is it an offdiagonal element in matrix-valued function"],
)

const PMaxEnt  = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to determine the optimized α parameter"],
    "nalph"   => [missing, 1, :I64   , "Number of the α parameters"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
    "blur"    => [missing, 1, :F64   , "Shall we blur the kernel and spectrum"],
)

const PStochOM = Dict{String,ADT}(
)

const PStochAC = Dict{String,ADT}(
    "nfine"   => [missing, 1, :I64   , "Number of points for very fine mesh"],
    "ngamm"   => [missing, 1, :I64   , "Number of δ functions"],
    "nalph"   => [missing, 1, :I64   , "Number of the α parameters"],
    "nwarm"   => [missing, 1, :I64   , "Number of Monte Carlo warming-up steps"],
    "nstep"   => [missing, 1, :I64   , "Number of Monte Carlo sweeping steps"],
    "ndump"   => [missing, 1, :I64   , "Intervals for monitoring Monte Carlo sweeps"],
    "alpha"   => [missing, 1, :F64   , "Starting value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
)

"""
    inp_toml(f::String, key::String, necessary::Bool)

Parse the configuration file (in toml format). It reads only parts of
the configuration file, which depends on the value of `key`.
"""
function inp_toml(f::String, key::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)

        if haskey(dict, key)
            dict[key]
        else
            error("Do not have this key: $key in file: $f")
        end
    else
        if necessary
            error("Please make sure that the file $f really exists")
        else
            nothing
        end
    end
end

"""
    inp_toml(f::String, necessary::Bool)

Parse the configuration file (in toml format). It reads the whole file.
"""
function inp_toml(f::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)
        return dict
    else
        if necessary
            error("Please make sure that the file $f really exists")
        else
            nothing
        end
    end
end

"""
    fil_dict(cfg::Dict{String,Any})

Transfer configurations from dict `cfg` to internal dicts (including
`PCOMM`, `PMaxEnt`, `PStochOM`, and `PStochAC` etc).
"""
function fil_dict(cfg::Dict{String,Any})
    # For COMM block
    COMM = cfg["COMM"]
    for key in keys(COMM)
        if haskey(PCOMM, key)
            PCOMM[key][1] = COMM[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # For MaxEnt block
    if haskey(cfg, "MaxEnt")
        MaxEnt = cfg["MaxEnt"]
        for key in keys(MaxEnt)
            if haskey(PMaxEnt, key)
                PMaxEnt[key][1] = MaxEnt[key]
            else
                error("Sorry, $key is not supported currently")
            end
        end
    end

    # For StochOM block
    if haskey(cfg, "StochOM")
        StochOM = cfg["StochOM"]
        for key in keys(StochOM)
            if haskey(PStochOM, key)
                PStochOM[key][1] = StochOM[key]
            else
                error("Sorry, $key is not supported currently")
            end
        end
    end

    # For StochAC block
    if haskey(cfg, "StochAC")
        StochAC = cfg["StochAC"]
        for key in keys(StochAC)
            if haskey(PStochAC, key)
                PStochAC[key][1] = StochAC[key]
            else
                error("Sorry, $key is not supported currently")
            end
        end
    end
end

"""
    chk_dict()

Validate the correctness and consistency of configurations.

See also: [`fil_dict`](@ref), [`_v`](@ref).
"""
function chk_dict()
    @assert get_c("solver") in ("MaxEnt", "StochOM", "StochAC")
    @assert get_c("ktype") in ("fermi", "boson", "bsymm")
    @assert get_c("mtype") in ("flat", "gauss", "func")
    @assert get_c("grid") in ("ftime", "btime", "ffreq", "bfreq")
    @assert get_c("mesh") in ("linear", "tangent")
    @assert get_c("ngrid") ≥ 1
    @assert get_c("nmesh") ≥ 1
    @assert get_c("wmax") > get_c("wmin")
    @assert get_c("beta") ≥ 0.0

    PA = [PCOMM]
    @cswitch get_c("solver") begin
        @case "MaxEnt"
            push!(PA, PMaxEnt)
            @assert get_m("method") in ("historic", "classic", "bryan", "chi2kink")
            @assert get_m("nalph") ≥ 1
            @assert get_m("alpha") > 0.0
            @assert get_m("ratio") > 0.0
            break

        @case "StochOM"
            push!(PA, PStochOM)
            break

        @case "StochAC"
            push!(PA, PStochAC)
            @assert get_a("nfine") ≥ 1000
            @assert get_a("ngamm") ≥ 100
            @assert get_a("nalph") ≥ 1
            @assert get_a("nwarm") ≥ 100
            @assert get_a("nstep") ≥ 1000
            @assert get_a("ndump") ≥ 100
            @assert get_a("alpha") > 0.0
            @assert get_a("ratio") > 0.0
            break
    end

    for P in PA
        foreach(x -> _v(x.first, x.second), P)
    end
end

"""
    _v(key::String, val::Array{Any,1})

Verify the value array. Called by chk_dict() function only.

See also: [`chk_dict`](@ref).
"""
@inline function _v(key::String, val::Array{Any,1})
    # To check if the value is updated
    if isa(val[1], Missing) && val[2] > 0
        error("Sorry, key ($key) shoule be set")
    end

    # To check if the type of value is correct
    if !isa(val[1], Missing) && !isa(val[1], eval(val[3]))
        error("Sorry, type of key ($key) is wrong")
    end
end

"""
    get_c(key::String)

Extract configurations from dict: PCOMM.
"""
@inline function get_c(key::String)
    if haskey(PCOMM, key)
        PCOMM[key][1]
    else
        error("Sorry, PCOMM does not contain key: $key")
    end
end

"""
    get_m(key::String)

Extract configurations from dict: PMaxEnt.
"""
@inline function get_m(key::String)
    if haskey(PMaxEnt, key)
        PMaxEnt[key][1]
    else
        error("Sorry, PMaxEnt does not contain key: $key")
    end
end

"""
    get_s(key::String)

Extract configurations from dict: PStochOM.
"""
@inline function get_s(key::String)
    if haskey(PStochOM, key)
        PStochOM[key][1]
    else
        error("Sorry, PStochOM does not contain key: $key")
    end
end

"""
    get_a(key::String)

Extract configurations from dict: PStochAC.
"""
@inline function get_a(key::String)
    if haskey(PStochAC, key)
        PStochAC[key][1]
    else
        error("Sorry, PStochAC does not contain key: $key")
    end
end
