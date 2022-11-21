# Output Files

Once the analytical continuation simulation is finished, the final spectral function ``A(\omega)`` is outputted to `Aout.data`. As is shown in Eq.~(\ref{eq:ImG}), ``A(\omega)`` is equivalent to the imaginary part of real frequency Green's function Im``G(\omega)``. Then the ACFlow toolkit will automatically calculate the corresponding real part Re``G(\omega)`` via the Kramers-Kronig transformation [see Eq.~(\ref{eq:kk})]. The full Green's function at real axis ``G(\omega)`` is stored in `Gout.data`. The spectral function is also used to reconstruct the imaginary time or Matsubara Green's functions [``\tilde{G}(\tau)`` or ``\tilde{G}(i\omega_n)``], which is stored in `{repr.data}. Besides the three output files, the ACFlow toolkit will generate quite a few output files, which can be used to analyze and diagnose the calculated results. All of the possible output files of the ACFlow toolkit are collected and explained in Table~\ref{tab:output}. 

| Filename | Description |
| :------- | :---------- |
|`Aout.data` | Final spectral function ``A(\omega)``. |
|`Gout.data` | Full Green's function at real axis ``G(\omega)``. |
|`repr.data` | Reproduced Green's function ``\tilde{G}`` at imaginary time or frequency axis. |
|`model.data` | Default model function ``m(\omega)``. |
|`chi2.data` | ``\log_{10}(\chi^2)`` vs ``\log_{10}(\alpha)``. |
|`prob.data` | ``P[\alpha\|\bar{G}]`` vs ``\alpha`` for the `MaxEnt` solver (`bryan` algorithm). |
|`Aout.data.alpha`_``i`` | ``\alpha``-resolved spectral function ``A_{\alpha}(\omega)`` for the `StochAC` solver. |
|`hamil.data` | ``U(\alpha)`` vs ``\alpha`` for the `StochAC` solver. |
|`goodness.dat` | ``\log_{10}(\chi^2)`` vs ``\log_{10}(\Theta)`` for the `StochSK` solver. |
|`stat.data` | Monte Carlo statistical information for stochastic sampling method. |

**Table** Possible output files of the ACFlow toolkit.
