#
# Project : Gardenia
# Source  : config.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/01
#

"""
    inp_toml(f::String, key::String, necessary::Bool)

Parse the configuration file (in toml format). It reads only parts of
the configuration file, which depends on the value of `key`.

### Arguments
* f -> Filename of configuration.
* key -> Parameter's name.
* necessary -> If it is true and configuration is absent, raise an error.

### Returns
* value -> Parameter's value.
"""
function inp_toml(f::String, key::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)
        #
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

### Arguments
* f -> Filename of configuration.
* necessary -> If it is true and configuration is absent, raise an error.

### Returns
* dict -> A Dictionary struct that contains all the key-value pairs.
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
`PBASE`, `PMaxEnt`, `PBarRat`, `PNevanAC`, `PStochAC`, `PStochSK`,
`PStochOM` and `PStochPX` etc). In other words, all the relevant internal
dicts should be filled / updated in this function.

### Arguments
* cfg -> A dict struct that contains all the configurations (from ac.toml).

### Returns
N/A
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

    # For BarRat block
    if haskey(cfg, "BarRat")
        BarRat = cfg["BarRat"]
        for key in keys(BarRat)
            if haskey(PBarRat, key)
                PBarRat[key][1] = BarRat[key]
            else
                error("Sorry, $key is not supported currently")
            end
        end
    end

    # For NevanAC block
    if haskey(cfg, "NevanAC")
        NevanAC = cfg["NevanAC"]
        for key in keys(NevanAC)
            if haskey(PNevanAC, key)
                PNevanAC[key][1] = NevanAC[key]
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
    see_dict()

Display all of the relevant configuration parameters to the terminal.

### Arguments
N/A

### Returns
N/A

See also: [`fil_dict`](@ref).
"""
function see_dict()
    println("[ Param: base ]")
    #
    println("finput  : ", get_b("finput") )
    println("solver  : ", get_b("solver") )
    println("ktype   : ", get_b("ktype")  )
    println("mtype   : ", get_b("mtype")  )
    println("grid    : ", get_b("grid")   )
    println("mesh    : ", get_b("mesh")   )
    println("ngrid   : ", get_b("ngrid")  )
    println("nmesh   : ", get_b("nmesh")  )
    println("wmax    : ", get_b("wmax")   )
    println("wmin    : ", get_b("wmin")   )
    println("beta    : ", get_b("beta")   )
    println("offdiag : ", get_b("offdiag"))
    println("fwrite  : ", get_b("fwrite") )
    println("pmodel  : ", get_b("pmodel") )
    println("pmesh   : ", get_b("pmesh")  )
    println("exclude : ", get_b("exclude"))
    #
    println()
    #
    println("[ Param: solver ]")
    #
    @cswitch get_b("solver") begin
        # For MaxEnt solver
        @case "MaxEnt"
            println("method : ", get_m("method"))
            println("stype  : ", get_m("stype") )
            println("nalph  : ", get_m("nalph") )
            println("alpha  : ", get_m("alpha") )
            println("ratio  : ", get_m("ratio") )
            println("blur   : ", get_m("blur")  )
            break

        # For BarRat solver
        @case "BarRat"
            println("atype   : ", get_r("atype")  )
            println("denoise : ", get_r("denoise"))
            println("epsilon : ", get_r("epsilon"))
            println("pcut    : ", get_r("pcut")   )
            println("eta     : ", get_r("eta")    )
            break

        # For NevanAC solver
        @case "NevanAC"
            println("pick  : ", get_n("pick") )
            println("hardy : ", get_n("hardy"))
            println("hmax  : ", get_n("hmax") )
            println("alpha : ", get_n("alpha"))
            println("eta   : ", get_n("eta")  )
            break

        # For StochAC solver
        @case "StochAC"
            println("nfine : ", get_a("nfine"))
            println("ngamm : ", get_a("ngamm"))
            println("nwarm : ", get_a("nwarm"))
            println("nstep : ", get_a("nstep"))
            println("ndump : ", get_a("ndump"))
            println("nalph : ", get_a("nalph"))
            println("alpha : ", get_a("alpha"))
            println("ratio : ", get_a("ratio"))
            break

        # For StochSK solver
        @case "StochSK"
            println("method : ", get_k("method"))
            println("nfine  : ", get_k("nfine") )
            println("ngamm  : ", get_k("ngamm") )
            println("nwarm  : ", get_k("nwarm") )
            println("nstep  : ", get_k("nstep") )
            println("ndump  : ", get_k("ndump") )
            println("retry  : ", get_k("retry") )
            println("theta  : ", get_k("theta") )
            println("ratio  : ", get_k("ratio") )
            break

        # For StochOM solver
        @case "StochOM"
            println("ntry  : ", get_s("ntry") )
            println("nstep : ", get_s("nstep"))
            println("nbox  : ", get_s("nbox") )
            println("sbox  : ", get_s("sbox") )
            println("wbox  : ", get_s("wbox") )
            println("norm  : ", get_s("norm") )
            break

        # For StochPX solver
        @case "StochPX"
            println("method : ", get_x("method"))
            println("nfine  : ", get_x("nfine") )
            println("npole  : ", get_x("npole") )
            println("ntry   : ", get_x("ntry")  )
            println("nstep  : ", get_x("nstep") )
            println("theta  : ", get_x("theta") )
            println("eta    : ", get_x("eta")   )
            break
    end
    #
    println()
    #
    flush(stdout)
end

"""
    rev_dict_b(BASE::Dict{String,Any})

Setup the configuration dictionary: `PBASE`.

### Arguments
* BASE -> A dict struct that contains configurations from the [BASE] block.

### Returns
N/A

See also: [`PBASE`](@ref).
"""
function rev_dict_b(BASE::Dict{String,Any})
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
    rev_dict_b(BASE::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PBASE`.

### Arguments
* BASE -> A dict struct that contains configurations from the [BASE] block.

### Returns
N/A

See also: [`PBASE`](@ref).
"""
function rev_dict_b(BASE::Dict{String,Vector{Any}})
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
    rev_dict_m(S::MaxEntSolver, MaxEnt::Dict{String,Any})

Setup the configuration dictionary: `PMaxEnt`.

### Arguments
* S -> A MaxEntSolver struct.
* MaxEnt -> A dict struct that contains configurations from the [MaxEnt] block.

### Returns
N/A

See also: [`PMaxEnt`](@ref).
"""
function rev_dict_m(S::MaxEntSolver, MaxEnt::Dict{String,Any})
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
    rev_dict_m(S::MaxEntSolver, MaxEnt::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PMaxEnt`.

### Arguments
* S -> A MaxEntSolver struct.
* MaxEnt -> A dict struct that contains configurations from the [MaxEnt] block.

### Returns
N/A

See also: [`PMaxEnt`](@ref).
"""
function rev_dict_m(S::MaxEntSolver, MaxEnt::Dict{String,Vector{Any}})
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
    rev_dict_r(S::BarRatSolver, BarRat::Dict{String,Any})

Setup the configuration dictionary: `PBarRat`.

### Arguments
* S -> A BarRatSolver struct.
* BarRat -> A dict struct that contains configurations from the [BarRat] block.

### Returns
N/A

See also: [`PBarRat`](@ref).
"""
function rev_dict_r(S::BarRatSolver, BarRat::Dict{String,Any})
    for key in keys(BarRat)
        if haskey(PBarRat, key)
            PBarRat[key][1] = BarRat[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PBarRat)
end

"""
    rev_dict_r(S::BarRatSolver, BarRat::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PBarRat`.

### Arguments
* S -> A BarRatSolver struct.
* BarRat -> A dict struct that contains configurations from the [BarRat] block.

### Returns
N/A

See also: [`PBarRat`](@ref).
"""
function rev_dict_r(S::BarRatSolver, BarRat::Dict{String,Vector{Any}})
    for key in keys(BarRat)
        if haskey(PBarRat, key)
            PBarRat[key][1] = BarRat[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PBarRat)
end

"""
    rev_dict_n(S::NevanACSolver, NevanAC::Dict{String,Any})

Setup the configuration dictionary: `PNevanAC`.

### Arguments
* S -> A NevanACSolver struct.
* NevanAC -> A dict struct that contains configurations from the [NevanAC] block.

### Returns
N/A

See also: [`PNevanAC`](@ref).
"""
function rev_dict_n(S::NevanACSolver, NevanAC::Dict{String,Any})
    for key in keys(NevanAC)
        if haskey(PNevanAC, key)
            PNevanAC[key][1] = NevanAC[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PNevanAC)
end

"""
    rev_dict_n(S::NevanACSolver, NevanAC::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PNevanAC`.

### Arguments
* S -> A NevanACSolver struct.
* NevanAC -> A dict struct that contains configurations from the [NevanAC] block.

### Returns
N/A

See also: [`PNevanAC`](@ref).
"""
function rev_dict_n(S::NevanACSolver, NevanAC::Dict{String,Vector{Any}})
    for key in keys(NevanAC)
        if haskey(PNevanAC, key)
            PNevanAC[key][1] = NevanAC[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PNevanAC)
end

"""
    rev_dict_a(S::StochACSolver, StochAC::Dict{String,Any})

Setup the configuration dictionary: `PStochAC`.

### Arguments
* S -> A StochACSolver struct.
* StochAC -> A dict struct that contains configurations from the [StochAC] block.

### Returns
N/A

See also: [`PStochAC`](@ref).
"""
function rev_dict_a(S::StochACSolver, StochAC::Dict{String,Any})
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
    rev_dict_a(S::StochACSolver, StochAC::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochAC`.

### Arguments
* S -> A StochACSolver struct.
* StochAC -> A dict struct that contains configurations from the [StochAC] block.

### Returns
N/A

See also: [`PStochAC`](@ref).
"""
function rev_dict_a(S::StochACSolver, StochAC::Dict{String,Vector{Any}})
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
    rev_dict_k(S::StochSKSolver, StochSK::Dict{String,Any})

Setup the configuration dictionary: `PStochSK`.

### Arguments
* S -> A StochSKSolver struct.
* StochSK -> A dict struct that contains configurations from the [StochSK] block.

### Returns
N/A

See also: [`PStochSK`](@ref).
"""
function rev_dict_k(S::StochSKSolver, StochSK::Dict{String,Any})
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
    rev_dict_k(S::StochSKSolver, StochSK::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochSK`.

### Arguments
* S -> A StochSKSolver struct.
* StochSK -> A dict struct that contains configurations from the [StochSK] block.

### Returns
N/A

See also: [`PStochSK`](@ref).
"""
function rev_dict_k(S::StochSKSolver, StochSK::Dict{String,Vector{Any}})
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
    rev_dict_s(S::StochOMSolver, StochOM::Dict{String,Any})

Setup the configuration dictionary: `PStochOM`.

### Arguments
* S -> A StochOMSolver struct.
* StochOM -> A dict struct that contains configurations from the [StochOM] block.

### Returns
N/A

See also: [`PStochOM`](@ref).
"""
function rev_dict_s(S::StochOMSolver, StochOM::Dict{String,Any})
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
    rev_dict_s(S::StochOMSolver, StochOM::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochOM`.

### Arguments
* S -> A StochOMSolver struct.
* StochOM -> A dict struct that contains configurations from the [StochOM] block.

### Returns
N/A

See also: [`PStochOM`](@ref).
"""
function rev_dict_s(S::StochOMSolver, StochOM::Dict{String,Vector{Any}})
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
    rev_dict_x(S::StochPXSolver, StochPX::Dict{String,Any})

Setup the configuration dictionary: `PStochPX`.

### Arguments
* S -> A StochPXSolver struct.
* StochPX -> A dict struct that contains configurations from the [StochPX] block.

### Returns
N/A

See also: [`PStochPX`](@ref).
"""
function rev_dict_x(S::StochPXSolver, StochPX::Dict{String,Any})
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
    rev_dict_x(S::StochPXSolver, StochPX::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PStochPX`.

### Arguments
* S -> A StochPXSolver struct.
* StochPX -> A dict struct that contains configurations from the [StochPX] block.

### Returns
N/A

See also: [`PStochPX`](@ref).
"""
function rev_dict_x(S::StochPXSolver, StochPX::Dict{String,Vector{Any}})
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

### Arguments
N/A

### Returns
N/A

See also: [`fil_dict`](@ref), [`_v`](@ref).
"""
function chk_dict()
    @assert get_b("solver") in ("MaxEnt", "BarRat", "NevanAC", "StochAC", "StochSK", "StochOM", "StochPX")
    @assert get_b("ktype") in ("fermi", "boson", "bsymm")
    @assert get_b("mtype") in ("flat", "gauss", "1gauss", "2gauss", "lorentz", "1lorentz", "2lorentz", "risedecay", "file")
    @assert get_b("grid") in ("ftime", "fpart", "btime", "bpart", "ffreq", "ffrag", "bfreq", "bfrag")
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
            @assert get_m("stype") in ("sj", "br")
            @assert get_m("nalph") ≥ 1
            @assert get_m("alpha") > 0.0
            @assert get_m("ratio") > 0.0
            break

        # For BarRat solver
        @case "BarRat"
            push!(PA, PBarRat)
            # It does not support imaginary time data.
            # The Prony approximation doesn't support broken data.
            @assert get_b("grid") in ("ffreq", "ffrag", "bfreq", "bfrag")
            #
            @assert get_r("atype") in ("cont", "delta")
            @assert get_r("denoise") in ("none", "prony_s", "prony_o")
            @assert get_r("epsilon") ≥ 0.0
            @assert get_r("pcut")    ≥ 1e-6
            @assert get_r("eta")     ≥ 1e-8
            break

        # For NevanAC solver
        @case "NevanAC"
            push!(PA, PNevanAC)
            # It does not support imaginary time data.
            # It does not support bosonic frequency data directly.
            @assert get_b("grid") in ("ffreq", "ffrag")
            #
            @assert get_n("hmax")  ≥ 10
            @assert get_n("alpha") ≥ 1e-8
            @assert get_n("eta")   ≥ 1e-8
            break

        # For StochAC solver
        @case "StochAC"
            push!(PA, PStochAC)
            @assert get_b("mtype") == "flat"
            #
            @assert get_a("nfine") ≥ 1000
            @assert get_a("ngamm") ≥ 1
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
            @assert get_k("ngamm") ≥ 1
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
            #
            @assert get_s("ntry")  ≥ 200
            @assert get_s("nstep") ≥ 1000
            @assert get_s("nbox")  ≥ 2
            @assert get_s("sbox")  > 0.0
            @assert get_s("wbox")  > 0.0
            break

        # For StochPX solver
        @case "StochPX"
            push!(PA, PStochPX)
            # It does not support imaginary time data.
            @assert get_b("grid") in ("ffreq", "ffrag", "bfreq", "bfrag")
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

### Arguments
* key -> Key of parameter.
* val -> Value of parameter.

### Returns
N/A

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

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

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

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

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
    get_r(key::String)

Extract configurations from dict: PBarRat.

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

See also: [`PBarRat`](@ref).
"""
@inline function get_r(key::String)
    if haskey(PBarRat, key)
        PBarRat[key][1]
    else
        error("Sorry, PBarRat does not contain key: $key")
    end
end

"""
    get_n(key::String)

Extract configurations from dict: PNevanAC.

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

See also: [`PNevanAC`](@ref).
"""
@inline function get_n(key::String)
    if haskey(PNevanAC, key)
        PNevanAC[key][1]
    else
        error("Sorry, PNevanAC does not contain key: $key")
    end
end

"""
    get_a(key::String)

Extract configurations from dict: PStochAC.

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

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

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

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

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

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

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

See also: [`PStochPX`](@ref).
"""
@inline function get_x(key::String)
    if haskey(PStochPX, key)
        PStochPX[key][1]
    else
        error("Sorry, PStochPX does not contain key: $key")
    end
end
