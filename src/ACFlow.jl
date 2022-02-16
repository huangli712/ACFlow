#
# Project : Gardenia
# Source  : ACFlow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/16
#

"""
    ACFlow

ACFlow is a modern software package for solving the many-body analytical
continuation problem. It is usually used to convert the single-particle
or two-particle correlation functions from imaginary axis to real axis.
Now this package is under heavy development. **PLEASE USE IT AT YOUR OWN
RISK**.

ACFlow supports the following algorithms:

* Maximum Entropy Method (MEM)
* Stochastic Analytical Continuation (SAC)
* Stochastic Optimization Method (SOM)


"""
module ACFlow

#=
### *Using Standard Libraries*
=#

using Distributed
using LinearAlgebra
using Printf
using Dates
using Random
using TOML

#=
### *Using Third-party Libraries*
=#

using Einsum
using LsqFit

#=
### *Includes*
=#

include("global.jl")
include("types.jl")
include("util.jl")
include("grid.jl")
include("mesh.jl")
include("config.jl")
include("inout.jl")
include("kernel.jl")
include("model.jl")
include("maxent.jl")
include("sac.jl")
include("base.jl")

#=
### *Exports*
=#

export welcome
export overview
export newton
export secant
export read_param
export read_data
export setup_args
export precompute
export make_singular_space
export make_kernel
export make_model
export make_grid
export make_mesh
export solve

#=
### *PreCompile*
=#

export _precompile

"""
    _precompile()

Here, we would like to precompile the whole `ACFlow` package to reduce
the runtime latency and speed up the successive calculations.
"""
function _precompile()
    prompt("Loading...")

    # Get an array of the names exported by the `ACFlow` module
    nl = names(ACFlow)

    # Go through each name
    cf = 0 # Counter
    for i in eachindex(nl)
        # Please pay attention to that nl[i] is a Symbol, we need to
        # convert it into string and function, respectively.
        str = string(nl[i])
        fun = eval(nl[i])

        # For methods only (macros must be excluded)
        if fun isa Function && !startswith(str, "@")
            # Increase the counter
            cf = cf + 1

            # Extract the signature of the function
            # Actually, `types` is a Core.SimpleVector.
            types = typeof(fun).name.mt.defs.sig.types

            # Convert `types` from SimpleVector into Tuple
            # If length(types) is 1, the method is without arguments.
            T = ()
            if length(types) > 1
                T = tuple(types[2:end]...)
            end

            # Precompile them one by one
            # println(i, " -> ", str, " -> ", length(types), " -> ", T)
            precompile(fun, T)
            @printf("Function %15s (#%3i) is compiled.\r", str, cf)
        end
    end

    prompt("Well, ACFlow is compiled and loaded ($cf functions).")
    prompt("We are ready to go!")
    println()
    flush(stdout)
end

"""
    __init__()

This function would be executed immediately after the module is loaded
at runtime for the first time. It works at the REPL mode only.
"""
__init__() = begin
    isinteractive() && _precompile()
end

end # END OF MODULE
