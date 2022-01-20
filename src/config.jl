
#=
### *Customized Types*
=#

"Customized types. It is used to define the following dicts."
const DType = Any

"Customized types. It is used to define the following dicts."
const ADT = Array{DType,1}

const PCOMM    = Dict{String,ADT}(
    "solver" => [missing, 1, :String, ""],
    "kernel" => [missing, 1, :String, ""],
    "model"  => [missing, 1, :String, ""],
    "mesh"   => [missing, 1, :String, ""],
    "wmax"   => [missing, 1, :F64, ""],
    "wmin"   => [missing, 1, :F64, ""],
    "beta"   => [missing, 1, :F64, ""],
)

const PMaxEnt  = Dict{String,ADT}(
    "method" => [missing, 1, :String, ""],
    "alpha"  => [missing, 1, :F64, ""],
)

const PStochOM = Dict{String,ADT}(
)

const PStochAC = Dict{String,ADT}(
)


"""
    inp_toml(f::String, key::String, necessary::Bool)

Parse the configuration file (in toml format). It reads only parts of
the configuration file, which depends on the value of `key`.

See also: [`setup`](@ref).
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

See also: [`setup`](@ref).
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
`PCASE`, `PDFT`, `PDMFT`, `PIMP`, and `PSOLVER`).

See also: [`chk_dict`](@ref).
"""
function fil_dict(cfg::Dict{String,Any})
    # For case block
    #
    # Pay attention to that the case block only includes one child element
    PCASE["case"][1] = cfg["case"]

    # For dft block
    dft = cfg["dft"]
    for key in keys(dft)
        if haskey(PDFT, key)
            PDFT[key][1] = dft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # For dmft block
    dmft = cfg["dmft"]
    for key in keys(dmft)
        if haskey(PDMFT, key)
            PDMFT[key][1] = dmft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # For impurity block
    impurity = cfg["impurity"]
    for key in keys(impurity)
        if haskey(PIMP, key)
            PIMP[key][1] = impurity[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # For solver block
    solver = cfg["solver"]
    for key in keys(solver)
        if haskey(PSOLVER, key)
            PSOLVER[key][1] = solver[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
end
