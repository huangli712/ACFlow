## Powered by Julia

The ACFlow toolkit is developed with pure Julia language. Thanks to powerful type system and multiple dispatch paradigm of the Julia language, the seven different analytic continuation solvers are integrated into an united software architecture. Redundant codes are greatly reduced. It is quite easy to implement new analytic continuation solver or add new features to the existing solvers if necessary. Distributed computing is a built-in feature of Julia. So, it is straightforward to realize parallel calculations in the ACFlow toolkit. Now except for the `MaxEnt`, `BarRat`, and `NevanAC` solvers, all the other solvers are parallelized.

## Code Repository

The official code repository of the ACFlow toolkit is:

```text
https://github.com/huangli712/ACFlow
```

## Core Codes

The source codes of the ACFlow toolkit are placed in the `acflow/src` folder. Their functions are summarized in **Table 1**.

| Filename | Description |
| :------- | :---------- |
| `ACFlow.jl` | Entry of the ACFlow module. |
| `maxent.jl` | Maximum entropy method. |
| `rfa.jl`    | Barycentric rational function approximation. |
| `nac.jl`    | Nevanlinna analytical continuation. |
| `sac.jl`    | Stochastic analytic continuation (K. S. D. Beach's algorithm). |
| `san.jl`    | Stochastic analytic continuation (A. W. Sandvik's algorithm). |
| `som.jl`    | Stochastic optimization method. |
| `spx.jl`    | Stochastic pole expansion. |
| `global.jl` | Numerical and physical constants. |
| `types.jl`  | Basic data structures and computational parameters. |
| `base.jl`   | Driver for analytic continuation simulation. |
| `inout.jl`  | Read input data and write calculated results. |
| `config.jl` | Parse configuration file and extract computational parameters. |
| `math.jl`   | Root finding, numerical integration, interpolation, Einstein summation, and curve fitting. |
| `util.jl`   | Some utility functions. |
| `mesh.jl`   | Meshes for spectral density. |
| `grid.jl`   | Grids for input data. |
| `model.jl`  | Default model functions. |
| `kernel.jl` | Kernel functions. |

**Table 1 |** List of source codes of the ACFlow toolkit.

!!! note

    There are four more scripts in the `acflow/util` folder.

    * `acrun.jl` - It is used to launch the analytic continuation tasks sequentially.
    * `acprun.jl` - It is used to launch the analytic continuation tasks parallelly.
    * `gmesh.jl` -> It will generate a dynamical mesh for next analytic continuation simulation.
    * `ppole.jl` - It is used to postprocess the outputs by the StochPX solver.

## Documentation

 The documentation of the ACFlow toolkit is written by using the `Markdown` language and the `Documenter.jl` package. The source codes are placed in the `acflow/docs` folder. The users can build documentation by themselves. Please see [Installation](@ref install) for how to do that. Or they can read the latest documentation in the following website:

```text
https://huangli712.github.io/projects/acflow/index.html
```

## Tests and Examples

Forty-eight tests and four tutorials are also shipped with the ACFlow toolkit. The source codes for internal tests are placed in the `acflow/test` folder, while those for tutorials are saved in the `acflow/tutor` folder. See `acflow/test/test.md` and `acflow/tutor/tutor.md` for more details.

!!! note

    We also develop an individual package, namely `ACTest`, to benchmark various analytic continuation solvers or methods. Please check the following URL for more details:

    ```text
    https://github.com/huangli712/ACTest
    ```
