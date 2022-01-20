
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
