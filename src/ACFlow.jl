#
# Project : Gardenia
# Source  : ACFlow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/17
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
### *Includes And Exports* : *global.jl*
=#

#=
*Summary* :

Define some type aliases and string constants for the ACFlow package.

*Members* :

```text
I32, I64    -> Numerical types (Integer).
F32, F64    -> Numerical types (Float).
C32, C64    -> Numerical types (Complex).
R32, R64    -> Numerical types (Union of Integer and Float).
N32, N64    -> Numerical types (Union of Integer, Float, and Complex).
#
__LIBNAME__ -> Name of this julia package.
__VERSION__ -> Version of this julia package.
__RELEASE__ -> Released date of this julia package.
__AUTHORS__ -> Authors of this julia package.
#
authors     -> Print the authors of ACFlow to screen.
```
=#

#
include("global.jl")
#
export I32, I64
export F32, F64
export C32, C64
export R32, R64
export N32, N64
#
export __LIBNAME__
export __VERSION__
export __RELEASE__
export __AUTHORS__
#
export authors

#=
### *Includes And Exports* : *types.jl*
=#

#=
*Summary* :

Define some dicts and structs, which are used to store the config
parameters or represent some essential data structures.

*Members* :

```text
DType           -> Customized type.
ADT             -> Customized type.
PCOMM           ->
PMaxEnt         ->
PStochAC        ->
PStochOM        ->
#
AbstractSolver  ->
MaxEntSolver    ->
StochACSolver   -> 
StochOMSolver   -> 
#
AbstractData
RawData
GreenData
#
AbstractGrid
FermionicImaginaryTimeGrid
FermionicMatsubaraGrid
BosonicImaginaryTimeGrid
BosonicMatsubaraGrid
#
AbstractMesh
LinearMesh
TangentMesh
```
=#

#
include("types.jl")
#
export DType
export ADT
export PCOMM
export PMaxEnt
export PStochAC
export PStochOM
#
export AbstractSolver
export MaxEntSolver
export StochACSolver
export StochOMSolver
#
export AbstractData
export RawData
export GreenData
#
export AbstractGrid
export FermionicImaginaryTimeGrid
export FermionicMatsubaraGrid
export BosonicImaginaryTimeGrid
export BosonicMatsubaraGrid
#
export AbstractMesh
export LinearMesh
export TangentMesh

#=
### *Includes And Exports* : *util.jl*
=#

#=
*Summary* :

To provide some useful utility macros and functions. They can be used
to colorize the output strings, query the environments, and parse the
input strings, etc.

*Members* :

```text
@cswitch      -> C-style switch.
@time_call    -> Evaluate a function call and print the elapsed time.
@pcs          -> Print colorful strings.
require       -> Check julia envirnoment.
setup_args    -> Setup ARGS manually.
query_args    -> Query program's arguments.
query_case    -> Query case (job's name).
query_inps    -> Query input files.
query_stop    -> Query case.stop file.
query_test    -> Query case.test file.
query_home    -> Query home directory of Zen framework.
query_core    -> Query home directory of ZenCore (where is ZenCore.jl).
query_dft     -> Query home directory of DFT engine.
query_dmft    -> Query home directory of DMFT engine.
query_solver  -> Query home directory of quantum impurity solvers.
is_vasp       -> Test whether the DFT backend is the vasp code.
is_qe         -> Test whether the DFT backend is the quantum espresso code.
is_plo        -> Test whether the projector is the projected local orbitals.
is_wannier    -> Test whether the projector is the wannier functions.
welcome       -> Print welcome message.
overview      -> Print runtime information of ZenCore.
goodbye       -> Say goodbye.
sorry         -> Say sorry.
prompt        -> Print some messages or logs to the output devices.
line_to_array -> Convert a line to a string array.
line_to_cmplx -> Convert a line to a cmplx number.
erf           -> Gauss error function.
subscript     -> Convert a number to subscript.
str_to_struct -> Convert a string to an instance of specified struct.
```
=#

#
include("util.jl")
#
export @cswitch
export @time_call
export @pcs
export require
export setup_args
export query_args
export welcome
export overview
export welcome
export overview
export goodbye
export sorry
export prompt
export line_to_array
export secant
export newton
export trapz
export simpson
export AbstractInterpolation
export LinearInterpolation
export QuadraticInterpolation
export CubicSplineInterpolation
export munge_data

#=
### *Includes And Exports* : *grid.jl*
=#

#
include("grid.jl")
#
export rebuild

#=
### *Includes And Exports* : *mesh.jl*
=#

#
include("mesh.jl")
#

#=
### *Includes And Exports* : *config.jl*
=#

#=
*Summary* :

To extract, parse, verify, and print the configuration parameters.
They are stored in external files (case.toml) or dictionaries.

*Members* :

```text
setup    -> Setup parameters.
renew    -> Renew some parameters dynamically.
inp_toml -> Parse case.toml, return raw configuration information.
fil_dict -> Fill dicts for configuration parameters.
rev_dict -> Update dicts for configuration parameters.
chk_dict -> Check dicts for configuration parameters.
exhibit  -> Display parameters for reference.
_v       -> Verify dict's values.
cat_c    -> Print dict (PCASE dict).
cat_d    -> Print dict (PDFT dict).
cat_m    -> Print dict (PDMFT dict).
cat_i    -> Print dict (PIMP dict).
cat_s    -> Print dict (PSOLVER dict).
get_c    -> Extract value from dict (PCASE dict), return raw value.
get_d    -> Extract value from dict (PDFT dict), return raw value.
get_m    -> Extract value from dict (PDMFT dict), return raw value.
get_i    -> Extract value from dict (PIMP dict), return raw value.
get_s    -> Extract value from dict (PSOLVER dict), return raw value.
str_c    -> Extract value from dict (PCASE dict), return string.
str_d    -> Extract value from dict (PDFT dict), return string.
str_m    -> Extract value from dict (PDMFT dict), return string.
str_i    -> Extract value from dict (PIMP dict), return string.
str_s    -> Extract value from dict (PSOLVER dict), return string.
```
=#

#
include("config.jl")
#
export inp_toml
export fil_dict
export rev_dict
export chk_dict
export _v
export get_c
export get_m
export get_a
export get_s

#
include("inout.jl")
#
export read_real_data
export read_complex_data

include("kernel.jl")
export build_kernel
export make_blur

include("model.jl")
include("maxent.jl")
include("sac.jl")
include("base.jl")




export read_param
export read_data

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
    #isinteractive() && _precompile()
    _precompile()
end

end # END OF MODULE
