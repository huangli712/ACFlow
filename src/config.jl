#
# Project : Gardenia
# Source  : config.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/03
#

#=
### *Customized Types*
=#

"Customized types. It is used to define the following dicts."
const DType = Any

"Customized types. It is used to define the following dicts."
const ADT = Array{DType,1}

const PCOMM    = Dict{String,ADT}(
    "finput"  => [missing, 1, :String, ""],
    "solver"  => [missing, 1, :String, ""],
    "kernel"  => [missing, 1, :String, ""],
    "model"   => [missing, 1, :String, ""],
    "grid"    => [missing, 1, :String, ""],
    "mesh"    => [missing, 1, :String, ""],
    "ngrid"   => [missing, 1, :I64, ""],
    "nmesh"   => [missing, 1, :I64, ""],
    "wmax"    => [missing, 1, :F64, ""],
    "wmin"    => [missing, 1, :F64, ""],
    "beta"    => [missing, 1, :F64, ""],
    "offdiag" => [missing, 1, :Bool, ""],
)

const PMaxEnt  = Dict{String,ADT}(
    "method" => [missing, 1, :String, ""],
    "alpha"  => [missing, 1, :F64, ""],
    "blur"   => [missing, 1, :F64, ""],
)

const PStochOM = Dict{String,ADT}(
)

const PStochAC = Dict{String,ADT}(
    "nfine" => [missing, 1, :I64, ""],
    "ngamm" => [missing, 1, :I64, ""],
    "nalph" => [missing, 1, :I64, ""],
    "nwarm" => [missing, 1, :I64, ""],
    "nstep" => [missing, 1, :I64, ""],
    "ndump" => [missing, 1, :I64, ""],
    "alpha" => [missing, 1, :F64, ""],
    "ratio" => [missing, 1, :F64, ""],
)

function read_param()
    cfg = inp_toml(query_args(), true)
    fil_dict(cfg)
end

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
`PCOMM`, `PMaxEnt`, `PStochOM`, and `PStochAC`).
"""
function fil_dict(cfg::Dict{String,Any})
    # For dft block
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

function chk_dict()
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