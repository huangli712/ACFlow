#
# Project : Gardenia
# Source  : config.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/04
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
    "ktype"   => [missing, 1, :String, "Type of kernel"],
    "mtype"   => [missing, 1, :String, "Type of default model"],
    "grid"    => [missing, 1, :String, "Grid of input correlation function"],
    "mesh"    => [missing, 1, :String, "Mesh of output spectrum"],
    "ngrid"   => [missing, 1, :I64   , "Number of grid points"],
    "nmesh"   => [missing, 1, :I64   , "Number of mesh points"],
    "wmax"    => [missing, 1, :F64   , "Maximum value of mesh"],
    "wmin"    => [missing, 1, :F64   , "Minimum value of mesh"],
    "beta"    => [missing, 1, :F64   , "Inverse temperature"],
    "fermi"   => [missing, 1, :Bool  , "Statistics of input correlation function"],
    "offdiag" => [missing, 1, :Bool  , "Whether it is an offdiagonal element in matrix-valued function"],
)

const PMaxEnt  = Dict{String,ADT}(
    "method"  => [missing, 1, :String, "How to determine the optimized α parameter"],
    "alpha0"  => [missing, 1, :F64   , "Starting value for the α parameter"],
    "alpha1"  => [missing, 1, :F64   , "Ending value for the α parameter"],
    "ratio"   => [missing, 1, :F64   , "Scaling factor for the α parameter"],
    "blur"    => [missing, 1, :F64   , "Shall we blur the kernel and spectrum"],
)

const PStochOM = Dict{String,ADT}(
)

const PStochAC = Dict{String,ADT}(
    "nfine"   => [missing, 1, :I64   , "Number of points for very fine mesh"],
    "ngamm"   => [missing, 1, :I64   , "Number of δ functions"],
    "nalph"   => [missing, 1, :I64   , "Number of α parameters"],
    "nwarm"   => [missing, 1, :I64   , ""],
    "nstep"   => [missing, 1, :I64   , ""],
    "ndump"   => [missing, 1, :I64   , ""],
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