The ACFlow toolkit is designed to be flexible and easy-to-use. It provides three running modes to facilitate analytic continuation calculations, namely the interactive, script, and standard modes.

## Interactive Mode

With the ACFlow toolkit, the users can configure and carry out analytic continuation simulations interactively in Julia's REPL (Read-Eval-Print Loop) environment. For example,

```julia-repl
julia> using ACFlow
julia> setup_args("ac.toml")
julia> read_param()
julia> mesh, Aout, Gout = solve(read_data())
```

Here, `ac.toml` is a configuration file, which contains essential computational parameters. The return values of the `solve()` function (i.e., `mesh`, `Aout`, and `Gout`) are mesh at real axis ``\omega``, spectral density ``A(\omega)``, and reproduced Green's function ``\tilde{G}``, respectively. They can be further analyzed or visualized by the users.

## Script Mode

The core functionalities of the ACFlow toolkit are exposed to the users via a simple application programming interface. So, the users can write Julia scripts easily by themselves to perform analytic continuation simulations. A minimal Julia script (`acrun.jl`) is listed as follows:

```julia
#!/usr/bin/env julia

using ACFlow

setup_args("ac.toml")
read_param()
mesh, Aout, Gout = solve(read_data())
```

Of course, this script can be extended to finish complex tasks. Later, a realistic example will be provided to show how to complete an analytic continuation of Matsubara self-energy function via the script mode (See [`Matsubara Self-Energy Function`](@ref ex_sigma)).

## Standard Mode

In the standard mode, the users have to prepare the input data manually. In addition, a configuration file must be provided. Supposed that the configuration file is `ac.toml`, then the analytic continuation calculation is launched as follows:

```shell
$ /home/your_home/acflow/util/acrun.jl ac.toml
```

or

```shell
$ /home/your_home/acflow/util/acprun.jl ac.toml
```

Noted that the `acrun.jl` script runs sequentially, while the `acprun.jl` script supports parallel and distributed computing. The two scripts are in the `acflow/util` folder. As we can conclude from the filename extension of configuration file (`ac.toml`), it adopts the `TOML` specification. The users may edit it with any text-based editors. Next we will introduce syntax and format of the input data files and configuration files.

## Parallel Calculations

Besides the `MaxEnt` solver, the computational efficiencies of the `StochAC`, `StochSK`, `StochOM`, and `StochPX` solvers are rather low. So, these solvers are parallelized to accelerate the analytic continuation simulations. The ACFlow toolkit provides a script, namely `acprun.jl`, to drive parallel calculations. Now the users should specify the number of parallel workers in this script:

```julia
#!/usr/bin/env julia
...
using Distributed # Julia's package to support distributed computing
...
addprocs(8)       # Now the number of parallel workers is 8. A total of 9
                  # processes are launched (8 workers + 1 master process).
...
```

It is limited by the available computational resources. A minimal PBS script is shown as follows:

```shell
#!/bin/bash
#PBS -N ACFlow
#PBS -l nodes=1:ppn=9
#PBS -q score
...
/home/your_home/acflow/util/acprun.jl ac.toml > nohup.dat 2>&1 # Please fix acprun.jl's path.
```

It is used to submit parallel jobs to computer clusters. Be careful, in order to maintain load balancing, the number of allocated CPUs should be larger than the number of parallel workers.

!!! note

    Now the `MaxEnt`, `BarRat`, and `NevanAC` solvers don't support parallel calculations.
