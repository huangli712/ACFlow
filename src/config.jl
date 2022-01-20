
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