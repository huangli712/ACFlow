#
# Project : Gardenia
# Source  : ACFlow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/05/05
#

"""
    ACFlow

ACFlow is a modern software package for solving the many-body analytical
continuation problem. It is usually used to convert the single-particle
or two-particle correlators from imaginary axis to real axis. Now this
package is under heavy development. **PLEASE USE IT AT YOUR OWN RISK**.

Nowadays the ACFlow package supports the following algorithms:

* Maximum Entropy Method (`MaxEnt` solver)
* Stochastic Analytical Continuation (`StochAC` solver)
* Stochastic Optimization Method (`StochOM` solver)

More algorithms will be implemented soon.

Note that ACFlow package has been designed to be integrated into the
`Zen` package. Actually, it is also compatible with the `iQIST` package.
In the other words, the end user can use it to analytically continue the
imaginary time (or Matsubara frequency) data generated by the various
quantum impurity solvers in the `iQIST` package. Of course, it is quite
easy to implement some kinds of interfaces for the other quantum impurity
solvers.

For more details about how to obtain, install and use the ACFlow package,
please visit the following website:

* `https://huangli712.github.io/projects/acflow/index.html`

Any suggestions, comments, and feedbacks are welcome. Enjoy it!
"""
module ACFlow

#=
### *Using Standard Libraries*
=#

using Distributed
using LinearAlgebra
using Statistics
using Random
using Dates
using Printf
using DelimitedFiles
using TOML

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
PBASE           -> Dict for general setup.
PMaxEnt         -> Dict for MaxEnt solver.
PStochAC        -> Dict for StochAC solver.
PStochOM        -> Dict for StochOM solver.
#
AbstractSolver  -> Abstract AC solver.
MaxEntSolver    -> It represents the MaxEnt solver.
StochACSolver   -> It represents the StochAC solver.
StochOMSolver   -> It represents the StochOM solver.
#
AbstractData    -> Abstract input data in imaginary axis.
RawData         -> Raw input data.
GreenData       -> Preprocessed input data.
#
AbstractGrid    -> Abstract mesh for input data.
FermionicImaginaryTimeGrid -> Grid in fermionic imaginary time axis.
FermionicMatsubaraGrid -> Grid in fermionic Matsubara frequency axis.
BosonicImaginaryTimeGrid -> Grid in bosonic imaginary time axis.
BosonicMatsubaraGrid -> Grid in bosonic Matsubara frequency axis.
#
AbstractMesh    -> Abstract grid for calculated spectral function.
LinearMesh      -> Linear mesh.
TangentMesh     -> Tangent mesh.
LorentzMesh     -> Lorentzian mesh.
HalfLorentzMesh -> Lorentzian mesh at half-positive axis.
#
AbstractMC      -> Abstract monte carlo engine.
StochACMC       -> Monte carlo engine used in the StochAC solver.
StochOMMC       -> Monte Carlo engine used in the StochOM solver.
```
=#

#
include("types.jl")
#
export DType
export ADT
export PBASE
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
export LorentzMesh
export HalfLorentzMesh
#
export AbstractMC
export StochACMC
export StochOMMC

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
#
require       -> Check julia envirnoment.
setup_args    -> Setup ARGS manually.
query_args    -> Query program's arguments.
welcome       -> Print welcome message.
overview      -> Print runtime information of ACFlow.
goodbye       -> Say goodbye.
sorry         -> Say sorry.
prompt        -> Print some messages or logs to the output devices.
line_to_array -> Convert a line to a string array.
```
=#

#
include("util.jl")
#
export @cswitch
export @time_call
export @pcs
#
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

#=
### *Includes And Exports* : *math.jl*
=#

#=
*Summary* :

To provide some numerical algorithms, such as numerical integrations,
interpolations, Einstein summation notation, and optimization method.

*Members* :

```text
secant    -> Root finding secant algorithm.
newton    -> Root finding newton algorithm.
#
trapz     -> Numerical integration (composite trapezoidal rule).
simpson   -> Numerical integration (simpson rule).
#
AbstractInterpolation -> Abstract struct for interpolation.
LinearInterpolation -> Linear interpolation.
QuadraticInterpolation -> Quadratic interpolation.
CubicSplineInterpolation -> Cubic spline interpolation.
#
@einsum   -> Macro for Einstein summation notation.
#
curve_fit -> Try to fit the given (x,y) data to a predefined model.
```
=#

#
include("math.jl")
#
export secant
export newton
#
export trapz
export simpson
#
export AbstractInterpolation
export LinearInterpolation
export QuadraticInterpolation
export CubicSplineInterpolation
#
export @einsum
#
export curve_fit

#=
### *Includes And Exports* : *grid.jl*
=#

#=
*Summary* :

To implement various grid for the input data.

*Members* :

```text
rebuild -> Rebuild the grid.
```
=#

#
include("grid.jl")
#
export rebuild

#=
### *Includes And Exports* : *mesh.jl*
=#

#=
*Summary* :

To implement various mesh for the calculated spectral functions.

*Members* :

```text
N/A
```
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
inp_toml -> Parse case.toml, return raw configuration information.
fil_dict -> Fill dicts for configuration parameters.
rev_dict -> Update dicts for configuration parameters.
chk_dict -> Check dicts for configuration parameters.
_v       -> Verify dict's values.
get_b    -> Extract value from dict (PBASE dict), return raw value.
get_m    -> Extract value from dict (PMaxEnt dict), return raw value.
get_a    -> Extract value from dict (PStochAC dict), return raw value.
get_s    -> Extract value from dict (PStochOM dict), return raw value.
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
export get_b
export get_m
export get_a
export get_s

#=
### *Includes And Exports* : *inout.jl*
=#

#=
*Summary* :

To read the input data or write the calculated results.

*Members* :

```text
read_real_data    -> Read data in imaginary time axis.
read_cmplx_data   -> Read data in Matsubara frequency axis.
#
write_spectrum    -> Write spectral functions.
write_model       -> Write default model function.
write_misfit      -> Write ????.
write_hamiltonian -> Write effective hamiltonian for StochAC solver.
write_probability -> Write Bayesian a-posteriori probability.
write_backward    -> Write reproduced input data in imaginary axis.
write_complete    -> Write full data in real axis.
write_statistics  -> Write statistical information for StochAC/StochOM solver.
```
=#

#
include("inout.jl")
#
export read_real_data
export read_cmplx_data
#
export write_spectrum
export write_model
export write_misfit
export write_hamiltonian
export write_probability
export write_backward
export write_complete
export write_statistics

#=
### *Includes And Exports* : *model.jl*
=#

#=
*Summary* :

To define some default model functions.

*Members* :

```text
build_flat_model         -> Construct a flat model.
build_gaussian_model     -> Construct a gaussian model.
build_1gaussian_model    -> Construct a shifted gaussian model.
build_2gaussians_model   -> Construct a two-gaussians model.
build_lorentzian_model   -> Construct a lorentzian model.
build_1lorentzian_model  -> Construct a shifted lorentzian model.
build_2lorentzians_model -> Construct a two-lorentzians model.
build_risedecay_model    -> Construct a rise-and-decay model.
build_file_model         -> Construct a model from file.
build_func_model         -> Construct a model by user-defined function.
```
=#

#
include("model.jl")
#
export build_flat_model
export build_gaussian_model
export build_1gaussian_model
export build_2gaussians_model
export build_lorentzian_model
export build_1lorentzian_model
export build_2lorentzians_model
export build_risedecay_model
export build_file_model
export build_func_model

#=
### *Includes And Exports* : *kernel.jl*
=#

#=
*Summary* :

To define various kernel functions.

*Members* :

```text
build_kernel        -> Build kernel function.
build_kernel_symm   -> Build kernel function for symmetric case.
#
make_blur           -> Add preblur effect to the spectral functions.
make_singular_space -> Perform singular value decomposition for kernel.
make_gauss_peaks    -> Generate a series gaussian peaks in a linear mesh.
```
=#

#
include("kernel.jl")
#
export build_kernel
export build_kernel_symm
#
export make_blur
export make_singular_space
export make_gauss_peaks

#=
### *Includes And Exports* : *maxent.jl*
=#

#=
*Summary* :

To implement the MaxEnt solver for analytically continuation problem.

*Members* :

```text
MaxEntContext -> Essential struct for the MaxEnt solver.
#
solve         -> Wrapper function for the MaxEnt solver.
init          -> Initialize maximum entropy simulation.
run           -> Perform maximum entropy simulation.
last          -> Postprocess the calculated results and write them.
#
historic      -> historic algorithm.
classic       -> Classic algorithm.
bryan         -> Bryan algorithm.
chi2kink      -> Chi2kink algorithm.
optimizer     -> Optimize the non-linear equation.
#
precompute    -> Precompute some key coefficients.
f_and_J       -> Define the function that need to be optimized.
f_and_J_offdiag -> Define the function that need to be optimized (offdiag version).
svd_to_real   -> From singular to real space.
svd_to_real_offdiag -> From singular to real space (offdiag version).
calc_entropy  -> Calculate entropy.
calc_entropy_offdiag -> Calculate entropy (offdiag version).
calc_bayes    -> Calculate Bayesian probability.
calc_bayes_offdiag -> Calculate Bayesian probability (offdiag version).
calc_chi2     -> Calculate ????.
```
=#

#
include("maxent.jl")
#
export MaxEntContext
#
export solve
export init
export run
export last
#
export historic
export classic
export bryan
export chi2kink
export optimizer
#
export precompute
export f_and_J
export f_and_J_offdiag
export svd_to_real
export svd_to_real_offdiag
export calc_entropy
export calc_entropy_offdiag
export calc_bayes
export calc_bayes_offdiag
export calc_chi2

#=
### *Includes And Exports* : *sac.jl*
=#

#=
*Summary* :

To implement the StochAC solver for analytically continuation problem.

*Members* :

```text
StochACElement -> A struct that contains monte carlo field configurations.
StochACContext -> Essential struct for the StochAC solver.
#
solve          -> Wrapper function for the StochAC solver.
init           -> Initialize stochastic analytical continuation simulation.
run (prun)     -> Perform stochastic analytical continuation simulation.
average        -> Evaluate the averaged results.
last           -> Postprocess the calculated results and write them.
#
warmup         -> Warmup monte carlo engine.
sample         -> Sample field configurations via metropolis algorithm.
measure        -> Measure spectral functions and internal energies.
#
init_mc        -> Create a StochACMC struct.
init_element   -> Create a StochACElement struct.
init_iodata    -> Preprocess the input data.
#
calc_fmesh     -> Build dense linear mesh in [wmin,wmax].
calc_xmesh     -> Build dense linear mesh in [0,1].
calc_phi       -> Calculate ?? function.
calc_delta     -> Precompute ?? functions.
calc_hamil     -> Calculate ??-resolved Hc.
calc_htau      -> Calculate ??-resolved h(??).
calc_alpha     -> Calculate ?? parameters.
#
try_mov1       -> Try to change the weights of ?? functions.
try_mov2       -> Try to shift the positions of ?? functions.
try_swap       -> Try to exchange configurations between two adjacent layers.
```
=#

#
include("sac.jl")
#
export StochACElement
export StochACContext
#
export solve
export init
export run
export prun
export average
export last
#
export warmup
export sample
export measure
#
export init_mc
export init_element
export init_iodata
#
export calc_fmesh
export calc_xmesh
export calc_phi
export calc_delta
export calc_hamil
export calc_htau
export calc_alpha
#
export try_mov1
export try_mov2
export try_swap

#=
### *Includes And Exports* : *som.jl*
=#

#=
*Summary* :

To implement the StochOM solver for analytically continuation problem.

*Members* :

```text
Box            -> A struct for describing the field configuration.
StochOMElement -> A struct that contains monte carlo field configurations.
StochOMContext -> Essential struct for the StochOM solver.
#
solve          -> Wrapper function for the StochOM solver.
init           -> Initialize stochastic optimization simulation.
run (prun)     -> Perform stochastic optimization simulation.
average        -> Evaluate the averaged results.
last           -> Postprocess the calculated results and write them.
#
update         -> Sample field configurations via metropolis algorithm.
#
init_mc        -> Create a StochOMMC struct.
init_element   -> Create a StochOMElement struct.
init_context   -> Prepare data for a StochOMContext struct.
init_iodata    -> Preprocess the input data.
#
calc_lambda    -> Build kernel function.
calc_error     -> Calculate ????.
calc_green     -> Reproduce green's function via the field configurations.
calc_norm      -> Calculate norm of the field configurations.
#
try_insert     -> Try to insert a new box in the configuration.
try_remove     -> Try to remove a box.
try_shift      -> Try to shift a box.
try_width      -> Try to change width of a box.
try_height     -> Try to change height of a box.
try_split      -> Try to split a box.
try_merge      -> Try to merge two boxes.
```
=#

#
include("som.jl")
#
export Box
export StochOMElement
export StochOMContext
#
export solve
export init
export run
export prun
export average
export last
#
export update
#
export init_mc
export init_element
export init_context
export init_iodata
#
export calc_lambda
export calc_error
export calc_green
export calc_norm
#
export try_insert
export try_remove
export try_shift
export try_width
export try_height
export try_split
export try_merge

#=
### *Includes And Exports* : *base.jl*
=#

#=
*Summary* :

To provide the basic workflow for end user of the ACFlow package.

*Members* :

```text
solve       -> Select solver to solve the analytical continuation problem.
#
reprod      -> Try to generate the input data via calculated spectrum.
kramers     -> Calculate real part of response function.
#
setup_param -> Setup parameters.
read_param  -> Read parameters from case.toml.
read_data   -> Read the input data.
#
make_data   -> Preprocess the input data.
make_grid   -> Generate grid for the input data.
make_mesh   -> Generate mesh for the calculated spectrum.
make_model  -> Generate default model function.
make_kernel -> Generate kernel function.
```
=#

#
include("base.jl")
#
export solve
#
export reprod
export kramers
#
export setup_param
export read_param
export read_data
#
export make_data
export make_grid
export make_mesh
export make_model
export make_kernel

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
