The input files for the ACFlow toolkit can be divided into two groups: data files and configuration files.

## Data Files

The input data should be stored in some text-based files, which adopt the space-separated data format. For imaginary time Green's function, the data file should contain three columns. They represent ``\tau``, ``\bar{G}(\tau)``, and standard deviation of ``\bar{G}(\tau)``. For fermionic Matsubara Green's function, the data file should contain five columns. They represent ``\omega_n``, Re``G(i\omega_n)``, Im``G(i\omega_n)``, standard deviation of Re``G(i\omega_n)``, and standard deviation of Im``G(i\omega_n)``. For bosonic correlation function ``\chi(i\omega_n)``, the data file should contain three columns. They represent ``\omega_n``, Re``\chi(i\omega_n)``, and standard deviation of Re``\chi(i\omega_n)``.

## Configuration Files

The configuration file adopts the TOML format. It is used to customize the computational parameters. It consists of one or more blocks. Possible blocks (or sections) of the configuration file include `[BASE]`, `[MaxEnt]`, `[BarRat]`, `[NevanAC]`, `[StochAC]`, `[StochSK]`, `[StochOM]`, and `[StochPX]`. The `[BASE]` block is mandatory, while the other blocks are optional. A schematic configuration file (`ac.toml`) is listed as follows:

```toml
[BASE]
finput = "giw.data"
solver = "StochOM"
...

[MaxEnt]
method = "chi2kink"
...

[BarRat]
atype  = "cont"
...

[NevanAC]
pick   = true
...

[StochAC]
nfine  = 10000
...

[StochSK]
method = "chi2min"
...

[StochOM]
ntry   = 100000
...

[StochPX]
method = "mean"
...
```

In the `[BASE]` block, the analytic continuation problem is defined. The solver used to solve the problem must be assigned. The types of mesh, grid, default model function, and kernel function are also determined. The `[MaxEnt]`, `[BarRat]`, `[NevanAC]`, `[StochAC]`, `[StochSK]`, `[StochOM]`, and `[StochPX]` blocks are used to customize the corresponding analytic continuation solvers further. In **Table 1**-**Table 8**, all the possible input parameters for these blocks are collected and summarized. As for detailed explanations of these parameters, please see [Parameters](@ref param).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`finput`  | string  | ''green.data'' | Filename for input data. |
|`solver`  | string  | ''MaxEnt''     | Solver for the analytic continuation problem. |
|`ktype`   | string  | ''fermi''      | Type of kernel function. |
|`mtype`   | string  | ''flat''       | Type of default model function. |
|`grid`    | string  | ''ffreq''      | Grid for input data (imaginary axis). |
|`mesh`    | string  | ''linear''     | Mesh for output data (real axis). |
|`ngrid`   | integer | 10             | Number of grid points. |
|`nmesh`   | integer | 501            | Number of mesh points. |
|`wmax`    | float   | 5.0            | Right boundary (maximum value) of mesh. |
|`wmin`    | float   | -5.0           | Left boundary (minimum value) of mesh. |
|`beta`    | float   | 10.0           | Inverse temperature. |
|`offdiag` | bool    | false          | Treat the off-diagonal part of matrix-valued function? |
|`fwrite`  | bool    | N/A            | Are the analytic continuation results written into files? |
|`pmodel`  | array   | N/A            | Additional parameters for customizing the default model. |
|`pmesh`   | array   | N/A            | Additional parameters for customizing the mesh. |
|`exclude` | array   | N/A            | Restriction of energy range of the spectrum. |

**Table 1 |** Possible parameters for the `[BASE]` block (See [[BASE] Block](@ref base_block)).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`method` | string  | ''chi2kink'' | How to determine the optimized ``\alpha`` parameter? |
|`stype`  | string  | ''sj''       | Type of the entropic factor. |
|`nalph`  | integer | 12           | Total number of the used ``\alpha`` parameters. |
|`alpha`  | float   | 1e9          | Starting value for the ``\alpha`` parameter. |
|`ratio`  | float   | 10.0         | Scaling factor for the ``\alpha`` parameter. |
|`blur`   | float   | -1.0         | Shall we preblur the kernel and spectrum? |

**Table 2 |** Possible input parameters for the `[MaxEnt]` block, which are used to configure the solver based on the maximum entropy method (See [[MaxEnt] Block](@ref maxent_block)).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`atype`    | string  | ''cont''  | Possible type of the spectrum. |
|`denoise`  | string  | ''prony'' | How to denoise the input data. |
|`epsilon`  | integer | 1e-10     | Threshold for the Prony approximation. |
|`pcut`     | float   | 1e-3      | Cutoff for unphysical poles. |
|`eta`      | float   | 1e-2      | Tiny distance from the real axis. |

**Table 3 |** Possible input parameters for the `[BarRat]` block, which are used to configure the solver based on the Barycentric rational function approximation (See [[BarRat] Block](@ref barrat_block)).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`pick`  | bool    | false | Check the Pick criterion or not. |
|`hardy` | bool    | false | Perform Hardy basis optimization or not. |
|`hmax`  | integer | 50    | Upper cut off of Hardy order. |
|`alpha` | float   | 1e-4  | Regulation parameter for smooth norm. |
|`eta`   | float   | 1e-2  | Tiny distance from the real axis. |

**Table 4 |** Possible input parameters for the `[NevanAC]` block, which are used to configure the solver based on the Nevanlinna analytical continuation (See [[NevanAC] Block](@ref nevanac_block)).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`nfine`  | integer | 10000   | Number of points of a very fine linear mesh. |
|`ngamm`  | integer | 512     | Number of ``\delta`` functions. |
|`nwarm`  | integer | 4000    | Number of Monte Carlo thermalization steps. |
|`nstep`  | integer | 4000000 | Number of Monte Carlo sweeping steps. |
|`ndump`  | integer | 40000   | Intervals for monitoring Monte Carlo sweeps. |
|`nalph`  | integer | 20      | Total number of the used ``\alpha`` parameters. |
|`alpha`  | float   | 1.0     | Starting value for the ``\alpha`` parameter. |
|`ratio`  | float   | 1.2     | Scaling factor for the ``\alpha`` parameter. |

**Table 5 |** Possible input parameters for the `[StochAC]` block, which are used to configure the solver based on the stochastic analytic continuation (Beach's algorithm. See [[StochAC] Block](@ref stochac_block)).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`method` | string  | ''chi2min'' | How to determine the optimized ``\Theta`` parameter? |
|`nfine`  | integer | 100000      | Number of points of a very fine linear mesh. |
|`ngamm`  | integer | 1000        | Number of ``\delta`` functions. |
|`nwarm`  | integer | 1000        | Number of Monte Carlo thermalization steps. |
|`nstep`  | integer | 20000       | Number of Monte Carlo sweeping steps. |
|`ndump`  | integer | 200         | Intervals for monitoring Monte Carlo sweeps. |
|`retry`  | integer | 10          | How often to recalculate the goodness-of-fit function. |
|`theta`  | float   | 1e6         | Starting value for the ``\Theta`` parameter. |
|`ratio`  | float   | 0.9         | Scaling factor for the ``\Theta`` parameter. |

**Table 6 |** Possible input parameters for the `[StochSK]` block, which are used to configure the solver based on the stochastic analytic continuation (Sandvik's algorithm. See [[StochSK] Block](@ref stochsk_block)).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`ntry`   | integer | 2000  | Number of attempts to figure out the solution. |
|`nstep`  | integer | 1000  | Number of Monte Carlo sweeping steps per try. |
|`nbox`   | integer | 100   | Number of rectangles to used construct the spectrum. |
|`sbox`   | float   | 0.005 | Minimum area of the randomly generated rectangles. |
|`wbox`   | float   | 0.02  | Minimum width of the randomly generated rectangles. |
|`norm`   | float   | -1.0  | Is the norm calculated? |

**Table 7 |** Possible input parameters for the `[StochOM]` block, which are used to configure the solver based on the stochastic optimization method (See [[StochOM] Block](@ref stochom_block)).

| Parameter | Type | Default | Description |
| :-------- | :--- | :------ | :---------- |
|`method` | string  | ''mean'' | How to evaluate the final spectral density? |
|`nfine`  | integer | 100000   | Number of points of a very fine linear mesh. |
|`npole`  | integer | 200      | Number of poles. |
|`ntry`   | integer | 1000     | Number of attempts to figure out the solution. |
|`nstep`  | integer | 1000000  | Number of Monte Carlo sweeping steps per attempt / try. |
|`theta`  | float   | 1e+6     | Artificial inverse temperature ``\Theta``. |
|`eta`    | float   | 1e-4     | Tiny distance from the real axis ``\eta``. |

**Table 8 |** Possible input parameters for the `[StochPX]` block, which are used to configure the solver based on the stochastic pole expansion (See [[StochPX] Block](@ref stochpx_block)).
