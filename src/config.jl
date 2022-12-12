#
# Project : Gardenia
# Source  : config.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/12/12
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
`PBASE`, `PMaxEnt`, `PStochAC`, `PStochSK`, and `PStochOM` etc).
"""
function fil_dict(cfg::Dict{String,Any})
    # For BASE block
    BASE = cfg["BASE"]
    for key in keys(BASE)
        if haskey(PBASE, key)
            PBASE[key][1] = BASE[key]
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

    # For StochSK block
    if haskey(cfg, "StochSK")
        StochSK = cfg["StochSK"]
        for key in keys(StochSK)
            if haskey(PStochSK, key)
                PStochSK[key][1] = StochSK[key]
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

    # For StochPX block
    if haskey(cfg, "StochPX")
        StochPX = cfg["StochPX"]
        for key in keys(StochPX)
            if haskey(PStochPX, key)
                PStochPX[key][1] = StochPX[key]
            else
                error("Sorry, $key is not supported currently")
            end
        end
    end
end

"""
    rev_dict(BASE::Dict{String,Any})

Setup the configuration dictionary: `PBASE`.

See also: [`PBASE`](@ref).
"""
function rev_dict(BASE::Dict{String,Any})
    for key in keys(BASE)
        if haskey(PBASE, key)
            PBASE[key][1] = BASE[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PBASE)
end

"""
    rev_dict(BASE::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PBASE`.

See also: [`PBASE`](@ref).
"""
function rev_dict(BASE::Dict{String,Vector{Any}})
    for key in keys(BASE)
        if haskey(PBASE, key)
            PBASE[key][1] = BASE[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PBASE)
end

"""
    rev_dict(S::MaxEntSolver, MaxEnt::Dict{String,Any})

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

"""
    rev_dict(S::MaxEntSolver, MaxEnt::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PMaxEnt`.

See also: [`PMaxEnt`](@ref).
"""
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

"""
    rev_dict(S::StochACSolver, StochAC::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochAC`.

See also: [`PStochAC`](@ref).
"""
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
    rev_dict(S::StochSKSolver, StochSK::Dict{String,Any})

Setup the configuration dictionary: `PStochSK`.

See also: [`PStochSK`](@ref).
"""
function rev_dict(S::StochSKSolver, StochSK::Dict{String,Any})
    for key in keys(StochSK)
        if haskey(PStochSK, key)
            PStochSK[key][1] = StochSK[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochSK)
end

"""
    rev_dict(S::StochSKSolver, StochSK::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochSK`.

See also: [`PStochSK`](@ref).
"""
function rev_dict(S::StochSKSolver, StochSK::Dict{String,Vector{Any}})
    for key in keys(StochSK)
        if haskey(PStochSK, key)
            PStochSK[key][1] = StochSK[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochSK)
end

"""
    rev_dict(S::StochOMSolver, StochOM::Dict{String,Any})

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

"""
    rev_dict(S::StochOMSolver, StochOM::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochOM`.

See also: [`PStochOM`](@ref).
"""
function rev_dict(S::StochOMSolver, StochOM::Dict{String,Vector{Any}})
    for key in keys(StochOM)
        if haskey(PStochOM, key)
            PStochOM[key][1] = StochOM[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochOM)
end

"""
    rev_dict(S::StochPXSolver, StochPX::Dict{String,Any})

Setup the configuration dictionary: `PStochPX`.

See also: [`PStochPX`](@ref).
"""
function rev_dict(S::StochPXSolver, StochPX::Dict{String,Any})
    for key in keys(StochPX)
        if haskey(PStochPX, key)
            PStochPX[key][1] = StochPX[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochPX)
end

"""
    rev_dict(S::StochPXSolver, StochPX::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochPX`.

See also: [`PStochPX`](@ref).
"""
function rev_dict(S::StochPXSolver, StochPX::Dict{String,Vector{Any}})
    for key in keys(StochPX)
        if haskey(PStochPX, key)
            PStochPX[key][1] = StochPX[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PStochPX)
end

"""
    chk_dict()

Validate the correctness and consistency of configurations.

See also: [`fil_dict`](@ref), [`_v`](@ref).
"""
function chk_dict()
    @assert get_b("solver") in ("MaxEnt", "StochAC", "StochSK", "StochOM", "StochPX")
    @assert get_b("ktype") in ("fermi", "boson", "bsymm")
    @assert get_b("mtype") in ("flat", "gauss", "1gauss", "2gauss", "lorentz", "1lorentz", "2lorentz", "risedecay", "file")
    @assert get_b("grid") in ("ftime", "btime", "ffreq", "bfreq")
    @assert get_b("mesh") in ("linear", "tangent", "lorentz", "halflorentz")
    @assert get_b("ngrid") ≥ 1
    @assert get_b("nmesh") ≥ 1
    @assert get_b("wmax") > get_b("wmin")
    @assert get_b("beta") ≥ 0.0

    PA = [PBASE]
    #
    @cswitch get_b("solver") begin
        # For MaxEnt solver
        @case "MaxEnt"
            push!(PA, PMaxEnt)
            #
            @assert get_m("method") in ("historic", "classic", "bryan", "chi2kink")
            @assert get_m("nalph") ≥ 1
            @assert get_m("alpha") > 0.0
            @assert get_m("ratio") > 0.0
            break

        # For StochAC solver
        @case "StochAC"
            push!(PA, PStochAC)
            @assert get_b("mtype") == "flat"
            #
            @assert get_a("nfine") ≥ 1000
            @assert get_a("ngamm") ≥ 100
            @assert get_a("nwarm") ≥ 100
            @assert get_a("nstep") ≥ 1000
            @assert get_a("ndump") ≥ 100
            @assert get_a("nalph") ≥ 1
            @assert get_a("alpha") > 0.0
            @assert get_a("ratio") > 0.0
            break

        # For StochSK solver
        @case "StochSK"
            push!(PA, PStochSK)
            #
            @assert get_k("method") in ("chi2min", "chi2kink")
            @assert get_k("nfine") ≥ 10000
            @assert get_k("ngamm") ≥ 1000
            @assert get_k("nwarm") ≥ 1000
            @assert get_k("nstep") ≥ 10000
            @assert get_k("ndump") ≥ 100
            @assert get_k("retry") ≥ 10
            @assert get_k("theta") > 1e+4
            @assert get_k("ratio") > 0.0
            break

        # For StochOM solver
        @case "StochOM"
            push!(PA, PStochOM)
            @assert get_b("grid") in ("btime", "ffreq", "bfreq")
            #
            @assert get_s("ntry")  ≥ 40
            @assert get_s("nstep") ≥ 1000
            @assert get_s("nbox")  ≥ 100
            @assert get_s("sbox")  > 0.0
            @assert get_s("wbox")  > 0.0
            break

        # For StochPX solver
        @case "StochPX"
            push!(PA, PStochPX)
            @assert get_b("grid") in ("ffreq", "bfreq")
            #
            @assert get_x("method") in ("best", "mean")
            @assert get_x("nfine") ≥ 10000
            @assert get_x("npole") ≥ 1
            @assert get_x("ntry")  ≥ 10
            @assert get_x("nstep") ≥ 100
            @assert get_x("theta") ≥ 0.00
            @assert get_x("eta")   ≥ 1e-8
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
    get_b(key::String)

Extract configurations from dict: PBASE.

See also: [`PBASE`](@ref).
"""
@inline function get_b(key::String)
    if haskey(PBASE, key)
        PBASE[key][1]
    else
        error("Sorry, PBASE does not contain key: $key")
    end
end

"""
    get_m(key::String)

Extract configurations from dict: PMaxEnt.

See also: [`PMaxEnt`](@ref).
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

See also: [`PStochAC`](@ref).
"""
@inline function get_a(key::String)
    if haskey(PStochAC, key)
        PStochAC[key][1]
    else
        error("Sorry, PStochAC does not contain key: $key")
    end
end

"""
    get_k(key::String)

Extract configurations from dict: PStochSK.

See also: [`PStochSK`](@ref).
"""
@inline function get_k(key::String)
    if haskey(PStochSK, key)
        PStochSK[key][1]
    else
        error("Sorry, PStochSK does not contain key: $key")
    end
end

"""
    get_s(key::String)

Extract configurations from dict: PStochOM.

See also: [`PStochOM`](@ref).
"""
@inline function get_s(key::String)
    if haskey(PStochOM, key)
        PStochOM[key][1]
    else
        error("Sorry, PStochOM does not contain key: $key")
    end
end

"""
    get_x(key::String)

Extract configurations from dict: PStochPX.

See also: [`PStochPX`](@ref).
"""
@inline function get_x(key::String)
    if haskey(PStochPX, key)
        PStochPX[key][1]
    else
        error("Sorry, PStochPX does not contain key: $key")
    end
end
