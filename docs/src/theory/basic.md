## Finite Temperature Green's Functions

Under the Wick's rotation ``t \to i\tau``, the time evolution operator in the Heisenberg picture ``e^{itH}`` will be replaced by ``e^{-\tau H}``. Such a transformation will increase efficiency of QMC random walking and suppress numerical oscillation (when ``t`` is large, the periodic oscillation of ``e^{itH}`` is quite obvious). This is an important reason why most of the finite temperature QMC algorithms are formulated in imaginary time axis. The outputs of finite temperature QMC simulations are usually single-particle or two-particle correlation functions. For example, the imaginary time Green's function ``G(\tau)`` is defined as follows:
```math
G(\tau) = \langle \mathcal{T}_{\tau} d(\tau) d^{\dagger}(0) \rangle,
```
where ``\tau`` denotes imaginary time, ``\mathcal{T}_{\tau}`` denotes time-ordered operator, and ``d`` and ``d^{\dagger}`` are annihilation and creation operators, respectively. The Matsubara Green's function ``G(i\omega_n)`` can be measured by QMC simulations or constructed from ``G(\tau)`` via direct Fourier transformation:
```math
G(i\omega_n) = \int^{\beta}_0 d\tau~e^{-i\omega_n \tau} G(\tau),
```
```math
G(\tau) = \frac{1}{\beta} \sum_n e^{i\omega_n \tau} G(i\omega_n).
```
Here, ``\beta`` means the inverse temperature (``\beta \equiv 1/T``) and ``\omega_n`` is the Matsubara frequency. Note that ``\omega_n`` is equal to ``(2n + 1) \pi / \beta`` for fermions and ``2n\pi/ \beta`` for bosons (``n`` is an integer).

---

## [Spectral Representation](@id spectrum)

Supposed that the spectral density of the single-particle Green's function is ``A(\omega)``, then we have:
```math
G(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 \pm e^{-\beta\omega}}
          A(\omega),
```
with the positive (negative) sign for fermionic (bosonic) operators. Similarly,
```math
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega
               \frac{1}{i\omega_n - \omega} A(\omega).
```
The two spectral equations denote the spectral representation of Green's function. We notice that the SPX method, as well as the other analytic continuation methods that are classified as ASM, are closely related to the spectral representation. Next we would like to make further discussions about this representation for the fermionic and bosonic correlators.

### Fermionic Correlators

The spectral density ``A(\omega)`` is defined on ``(-\infty,\infty)``. It is positive definite, i.e., ``A(\omega) \ge 0``. The spectral equations can be reformulated as:
```math
G(\tau) = \int^{+\infty}_{-\infty} d\omega~K(\tau,\omega) A(\omega),
```
and
```math
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega~K(\omega_n,\omega) A(\omega),
```
respectively. The kernel functions ``K(\tau,\omega)`` and ``K(\omega_n,\omega)`` are defined as follows:
```math
K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}},
```
and
```math
K(\omega_n,\omega) = \frac{1}{i\omega_n - \omega}.
```

### Bosonic Correlators

The spectral density ``A(\omega)`` obeys the following constraint: ``\text{sign}(\omega) A(\omega) \ge 0``. Thus, it is more convenient to define a new function ``\tilde{A}(\omega)`` where ``\tilde{A}(\omega) = A(\omega)/\omega``. Clearly, ``\tilde{A}(\omega)`` is always positive definite. As a result the two spectral equations can be rewritten as:
```math
G(\tau) = \int^{+\infty}_{-\infty} d\omega~
    K(\tau,\omega)\tilde{A}(\omega),
```
and
```math
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega~
    K(\omega_n,\omega) \tilde{A}(\omega),
```
respectively. Now the bosonic kernel ``K(\tau,\omega)`` becomes:
```math
K(\tau,\omega) = \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}}.
```
Especially, ``K(\tau,0) = 1/\beta``. As for ``K(\omega_n,\omega)``, its expression is:
```math
K(\omega_n,\omega) = \frac{\omega}{i\omega_n - \omega}.
```
Especially, ``K(0,0) = -1``. Besides the bosonic Green's function, typical correlator of this kind includes the transverse spin susceptibility ``\chi_{+-}(\tau) = \langle S_{+}(\tau) S_{-}(0) \rangle``, where ``S_{+} = S_x + iS_y`` and ``S_{-} = S_x - i S_y``.

### Bosonic Correlators of Hermitian Operators

There is a special case of the previous observable kind with ``d = d^{\dagger}``. Here, ``A(\omega)`` becomes an odd function, and equivalently, ``\tilde{A}(\omega)`` is an even function [i.e., ``\tilde{A}(\omega) = \tilde{A}(-\omega)``]. Therefore the limits of integrations in the two spectral equations are reduced from ``(-\infty,\infty)`` to ``(0,\infty)``. So the two equations can be transformed into:
```math
G(\tau) = \int^{+\infty}_{0} d\omega~K(\tau,\omega)\tilde{A}(\omega),
```
and
```math
G(i\omega_n) = \int^{+\infty}_{0} d\omega~K(\omega_n,\omega) \tilde{A}(\omega),
```
respectively. The corresponding ``K(\tau,\omega)`` reads:
```math
K(\tau,\omega) = \frac{\omega \left[e^{-\tau\omega} + e^{-(\beta - \tau)\omega}\right]}
                      {1 - e^{-\beta\omega}}.
```
Especially, ``K(\tau,0) = 2 / \beta``. And ``K(\omega_n,\omega)`` becomes:
```math
K(\omega_n, \omega) = \frac{-2\omega^2}{\omega_n^2 + \omega^2}.
```
Especially, ``K(0,0) = -2``. Perhaps the longitudinal spin susceptibility ``\chi_{zz}(\tau) = \langle S_z(\tau) S_z(0) \rangle`` and the charge susceptibility ``\chi_{ch}(\tau) = \langle N(\tau) N(0) \rangle`` are the most widely used observables of this kind.

As mentioned above, the kernel function is quite strange. The values of ``K(\tau,\omega)`` could change by tens of orders of magnitude. Especially, at large positive and negative frequencies, ``K(\tau,\omega)`` is exponentially small. It implies that at large ``|\omega|`` the features of ``A(\omega)`` depend upon the fine structures of ``G(\tau)``. However, the ``G(\tau)`` data provided by QMC simulations are always fluctuant and noisy. Tiny deviations in ``G(\tau)`` from its expected values can lead to enormous changes in ``A(\omega)``. Thus, analytic continuation is often characterized as an ill-posed problem.

In principle, for incomplete and noise ``G(\tau)`` or ``G(i\omega_n)``, the number of spectral functions ``A(\omega)`` that satisfy the spectral equations is infinite. So the question becomes which ``A(\omega)`` should be chosen. Now there are two different strategies to solve this problem. The first one is to choose the most likely ``A(\omega)``. The second one is to evaluate the average of all the candidate spectral functions.

---

## [Pole Representation](@id pole_repr)

It is well known that the finite temperature many-body Green's functions can be expressed within the Lehmann representation:
```math
G_{ab}(z) = \frac{1}{Z} \sum_{m,n}
\frac{\langle n | d_a | m \rangle \langle m | d_b^{\dagger} | n \rangle}{z + E_n - E_m}
\left(e^{-\beta E_n} \pm e^{-\beta E_m}\right),
```
where ``a`` and ``b`` are the band indices, ``d`` (``d^{\dagger}``) denote the annihilation (creation) operators, ``|n \rangle`` and ``|m \rangle`` are eigenstates of the Hamiltonian ``\hat{H}``, and ``E_n`` and ``E_m`` are the corresponding eigenvalues, ``Z`` is the partition function (``Z = \sum_n e^{-\beta E_n}``). The positive sign corresponds to fermions, while the negative sign corresponds to bosons. The domain of this function is on the complex plane, but the real axis is excluded (``z \in \{0\} \bigcup \mathbb{C}~\backslash~\mathbb{R}``). If ``z = i\omega_n \in i\mathbb{R}``, ``G_{ab}(i\omega_n)`` is the Matsubara Green's function. If ``z = \omega + i0^{+}``, ``G_{ab}(\omega + i0^{+}) = G_{ab}^{R}(\omega)`` is called the retarded Green's function.

At first we focus on the diagonal cases (``a = b``). For the sake of simplicity, the band indices are ignored in the following discussions. Let ``A_{mn} = \langle n | d | m \rangle \langle m | d^{\dagger} | n \rangle \left(e^{-\beta E_n} + e^{-\beta E_m}\right) / Z`` and ``P_{mn} = E_m - E_n``, then ``G(z) = \sum_{m,n} A_{mn} / (z - P_{mn})``. Clearly, only those nonzero elements of ``A_{mn}`` contribute to the Green's function. If the indices ``m`` and ``n`` are further compressed into ``\gamma`` (i.e, ``\gamma = \{m,n\}``), then the Lehmann representation is simplified to:
```math
G(z) = \sum^{N_p}_{\gamma = 1} \frac{A_{\gamma}}{z - P_{\gamma}}.
```
Here, ``A_{\gamma}`` and ``P_{\gamma}`` mean the amplitude and location of the ``\gamma``-th pole, respectively. ``N_p`` means the number of poles, which is equal to the total number of nonzero ``A_{mn}``. Such an analytic expression of Green's function is called the `pole expansion`. It is valid for both fermionic and bosonic correlators.

### [Fermionic Correlators](@id xi_fc)

For fermionic systems, the pole representation for Matsubara Green's function can be recast as:
```math
G(i\omega_n) = \sum^{N_p}_{\gamma = 1} \Xi(\omega_n, P_{\gamma}) A_{\gamma}.
```
Here, ``\Xi`` is the kernel matrix. It is evaluated by:
```math
\Xi(\omega_n, \omega) = \frac{1}{i\omega_n - \omega}.
```
Note that ``A_{\gamma}`` and ``P_{\gamma}`` should satisfy the following constraints:
```math
\forall \gamma, 0 \le A_{\gamma} \le 1, \sum_{\gamma} A_{\gamma} = 1, P_{\gamma} \in \mathbb{R}.
```

### Bosonic Correlators

For bosonic systems, the pole representation for Matsubara Green's function can be defined as follows:
```math
G(i\omega_n) = \sum^{N_p}_{\gamma=1} \Xi(\omega_n, P_{\gamma}) \tilde{A}_{\gamma}.
```
Here, ``\Xi`` is evaluated by:
```math
\Xi(\omega_n, \omega) = \frac{G_0 \omega}{i\omega_n - \omega},
```
where ``G_{0} = -G(i\omega_n = 0)``, which should be a positive real number. Be careful, ``\Xi(0,\omega) = -G_0``. ``\tilde{A}_{\gamma}`` is the renormalized amplitude of the ``\gamma``-th pole:
```math
\tilde{A}_{\gamma} = \frac{A_{\gamma}}{G_0 P_{\gamma}}.
```
Note that ``\tilde{A}_{\gamma}`` and ``P_{\gamma}`` should satisfy the following constraints:
```math
\forall \gamma, 0 \le \tilde{A}_{\gamma} \le 1, \sum_{\gamma} \tilde{A}_{\gamma} = 1, P_{\gamma} \in \mathbb{R}.
```

### Bosonic Correlators of Hermitian Operators

Its pole representation can be defined as follows (``\forall \gamma``, ``A_{\gamma} > 0`` and ``P_{\gamma} > 0``):
```math
\begin{align*}
G(i\omega_n) &= \sum^{N_p}_{\gamma = 1}
\left(
    \frac{A_{\gamma}}{i\omega_n - P_{\gamma}} -
    \frac{A_{\gamma}}{i\omega_n + P_{\gamma}}
\right) \\
&= \sum^{N_p}_{\gamma = 1}
\Xi(\omega_n, P_{\gamma}) \tilde{A}_{\gamma}.
\end{align*}
```
Thus, the kernel matrix ``\Xi`` reads:
```math
\Xi(\omega_n, \omega) = \frac{-G_0 \omega^2}{\omega^2_n + \omega^2}.
```
Especially, ``\Xi(0,0) = -G_0``. The renormalized weight ``\tilde{A}_{\gamma}`` reads:
```math
\tilde{A}_{\gamma} = \frac{2A_{\gamma}}{G_0 P_{\gamma}}.
```
The constraints for ``\tilde{A}_{\gamma}`` and ``P_{\gamma}`` are the same with those for general bosonic correlators.

### Matrix-Valued Green's Functions

As for the off-diagonal cases (``a \neq b``), it is lightly to prove that ``\sum_{\gamma} A_{\gamma} = 0``. It implies that there exist poles with negative weights. Hence we can split the poles into two groups according to the signs of their amplitudes. The Matsubara Green's function can be expressed as follows:
```math
\begin{align*}
G(i\omega_n) &= \sum^{N^{+}_p}_{\gamma = 1}
               \frac{A^{+}_{\gamma}}{i\omega_n - P^{+}_{\gamma}} -
               \sum^{N^{-}_p}_{\gamma = 1}
               \frac{A^{-}_{\gamma}}{i\omega_n - P^{-}_{\gamma}} \\
             &= \sum^{N^{+}_p}_{\gamma = 1}
               \Xi(\omega_n, P^{+}_{\gamma}) A^{+}_{\gamma} -
               \sum^{N^{-}_p}_{\gamma = 1}
               \Xi(\omega_n, P^{-}_{\gamma}) A^{-}_{\gamma}
.
\end{align*}
```
Here, ``\Xi(\omega_n, \omega)`` is already defined in the [Fermionic Correlators](@ref xi_fc) section. The ``A^{\pm}_{\gamma}`` and ``P^{\pm}_{\gamma}`` are restricted by the rules for fermionic correlators. In addition,
```math
N_p = N^{+}_p + N^{-}_p,
```
and
```math
\sum^{N^{+}_p}_{\gamma = 1} A^{+}_{\gamma} -
\sum^{N^{-}_p}_{\gamma = 1} A^{-}_{\gamma} = 0.
```

---

## Relevant Parameters

See [[BASE] Block](@ref base_block)
