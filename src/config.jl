#
# Project : Gardenia
# Source  : config.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/21
#

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
`PCOMM`, `PMaxEnt`, `PStochAC`, and `PStochOM` etc).
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
end

"""
    rev_dict(COMM::Dict{String,Any})
    rev_dict(COMM::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PCOMM`.

See also: [`PCOMM`](@ref).
"""
function rev_dict(COMM::Dict{String,Any})
    for key in keys(COMM)
        if haskey(PCOMM, key)
            PCOMM[key][1] = COMM[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PCOMM)
end

function rev_dict(COMM::Dict{String,Vector{Any}})
    for key in keys(COMM)
        if haskey(PCOMM, key)
            PCOMM[key][1] = COMM[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PCOMM)
end

"""
    rev_dict(S::MaxEntSolver, MaxEnt::Dict{String,Any})
    rev_dict(S::MaxEntSolver, MaxEnt::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PMaxEnt`.

See also: [`PMaxEnt`](@ref).
"""
function rev_dict(S::MaxEntSolver, MaxEnt::Dict{String,Any})
    for key in keys(MaxEnt)
        if haskey(PMaxEnt, key)
            PMaxEnt[key][1] = MaxEnt[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PMaxEnt)
end

function rev_dict(S::MaxEntSolver, MaxEnt::Dict{String,Vector{Any}})
    for key in keys(MaxEnt)
        if haskey(PMaxEnt, key)
            PMaxEnt[key][1] = MaxEnt[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PMaxEnt)
end

"""
    rev_dict(S::StochACSolver, StochAC::Dict{String,Any})
    rev_dict(S::StochACSolver, StochAC::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochAC`.

See also: [`PStochAC`](@ref).
"""
function rev_dict(S::StochACSolver, StochAC::Dict{String,Any})
    for key in keys(StochAC)
        if haskey(PStochAC, key)
            PStochAC[key][1] = StochAC[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochAC)
end

function rev_dict(S::StochACSolver, StochAC::Dict{String,Vector{Any}})
    for key in keys(StochAC)
        if haskey(PStochAC, key)
            PStochAC[key][1] = StochAC[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochAC)
end

"""
    rev_dict(S::StochOMSolver, StochOM::Dict{String,Any})
    rev_dict(S::StochOMSolver, StochOM::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochOM`.

See also: [`PStochOM`](@ref).
"""
function rev_dict(S::StochOMSolver, StochOM::Dict{String,Any})
    for key in keys(StochOM)
        if haskey(PStochOM, key)
            PStochOM[key][1] = StochOM[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochOM)
end

function rev_dict(S::StochOMSolver, StochOM::Dict{String,Vector{Any}})
    for key in keys(StochOM)
        if haskey(PStochOM, key)
            PStochOM[key][1] = StochOM[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochOM)
end

"""
    chk_dict()

Validate the correctness and consistency of configurations.

See also: [`fil_dict`](@ref), [`_v`](@ref).
"""
function chk_dict()
    @assert get_c("solver") in ("MaxEnt", "StochOM", "StochAC")
    @assert get_c("ktype") in ("fermi", "boson", "bsymm")
    @assert get_c("mtype") in ("flat", "gauss", "file")
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

        @case "StochAC"
            push!(PA, PStochAC)
            @assert get_a("nfine") ≥ 1000
            @assert get_a("ngamm") ≥ 100
            @assert get_a("nwarm") ≥ 100
            @assert get_a("nstep") ≥ 1000
            @assert get_a("ndump") ≥ 100
            @assert get_a("nalph") ≥ 1
            @assert get_a("alpha") > 0.0
            @assert get_a("ratio") > 0.0
            break

        @case "StochOM"
            push!(PA, PStochOM)
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
