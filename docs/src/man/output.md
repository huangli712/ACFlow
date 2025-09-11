Once the analytic continuation simulation is finished, the final spectral function ``A(\omega)`` is outputted to `Aout.data`. Note that ``A(\omega)`` is equivalent to the imaginary part of real frequency Green's function Im``G(\omega)``. Then the ACFlow toolkit will automatically calculate the corresponding real part Re``G(\omega)`` via the [Kramers-Kronig transformation](@ref spectrum). The full Green's function at real axis ``G(\omega)`` is stored in `Gout.data`. The spectral function is also used to reconstruct the imaginary time or Matsubara Green's functions [``\tilde{G}(\tau)`` or ``\tilde{G}(i\omega_n)``], which is stored in `repr.data`. Besides the three output files, the ACFlow toolkit will generate quite a few output files, which can be used to analyze and diagnose the calculated results. All of the possible output files of the ACFlow toolkit are collected and explained in **Table 1**.

| Filename | Description |
| :------- | :---------- |
|`Aout.data` | Final spectral function ``A(\omega)``. |
|`Aout.data.alpha`_``i`` | ``\alpha``-resolved spectral function ``A_{\alpha}(\omega)`` for the `StochAC` solver. |
|`repr.data` | Reproduced Green's function ``\tilde{G}`` at imaginary time or frequency axis. |
|`Gout.data` | Full Green's function at real axis ``G(\omega)``. |
|`chi2.data` | ``\log_{10}(\chi^2)`` vs ``\log_{10}(\alpha)``. |
|`goodness.dat` | ``\log_{10}(\chi^2)`` vs ``\log_{10}(\Theta)`` for the `StochSK` solver. |
|`model.data` | Default model ``m(\omega)``. |
|`prony.data` | Prony approximation to the Matsubara data (for the `BarRat` solver). |
|`Gprony.data` | Preprocessed Matsubara data by Prony approximation (for the `BarRat` solver). |
|`barycentric.data` | Barycentric rational function approximation to the Matsubara data (for the `BarRat` solver) |
|`hamil.data` | ``U(\alpha)`` vs ``\alpha`` for the `StochAC` solver. |
|`passed.data`| Indices of selected solutions for the `StochOM` and the `StochPX` solvers. |
|`pole.data` | Amplitudes and positions of the poles for the `StochPX` solver. |
|`prob.data` | ``P[\alpha\|\bar{G}]`` vs ``\alpha`` for the `MaxEnt` solver (`bryan` algorithm). |
|`stat.data` | Monte Carlo statistical information for stochastic sampling methods. |
|`err.out`| Error or exception messages (stacktrace) collected during simulation. |

**Table 1 |** Possible output files of the ACFlow toolkit.

!!! warning

    For bosonic systems, the `MaxEnt`, `StochAC`, `StochSK`, and `StochOM` solvers will generate and output ``\tilde{A}(\omega)``, instead of traditional ``A(\omega)``. That is to say, in `Aout.data`, the data are actually ``\tilde{A}(\omega)``. If the users want to retrieve ``A(\omega)``, they have to do the transformation by themselves:

    ```math
    \tilde{A}(\omega) = \frac{A(\omega)}{\omega},
    ```

    or resort to `Gout.data`. On the other hand, the `BarRat` and `StochPX` solvers will always generate and output ``A(\omega)``, irrespective of bosonic and fermionic systems.

    The `NevanAC` solver doesn't support bosonic system directly.
