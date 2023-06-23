# Basic Principles

## Finite Temperature Green's Functions

Under the Wick's rotation ``t \to i\tau``, the time evolution operator in the Heisenberg picture ``e^{itH}`` will be replaced by ``e^{-\tau H}``. Such a transformation will increase efficiency of QMC random walking and suppress numerical oscillation (when ``t`` is large, the periodic oscillation of ``e^{itH}`` is quite obvious). This is an important reason why most of the finite temperature QMC algorithms are formulated in imaginary time axis. The outputs of finite temperature QMC simulations are usually single-particle or two-particle correlation functions. For example, the imaginary time Green's function ``G(\tau)`` is defined as follows:
```math
\begin{equation}
G(\tau) = \langle \mathcal{T}_{\tau} d(\tau) d^{\dagger}(0) \rangle,
\end{equation}
```
where ``\tau`` denotes imaginary time, ``\mathcal{T}_{\tau}`` denotes time-ordered operator, and ``d`` and ``d^{\dagger}`` are annihilation and creation operators, respectively. The Matsubara Green's function ``G(i\omega_n)`` can be measured by QMC simulations or constructed from ``G(\tau)`` via direct Fourier transformation:
```math
\begin{equation}
G(i\omega_n) = \int^{\beta}_0 d\tau~e^{-i\omega_n \tau} G(\tau),
\end{equation}
```
```math
\begin{equation}
G(\tau) = \frac{1}{\beta} \sum_n e^{i\omega_n \tau} G(i\omega_n).
\end{equation}
```
Here, ``\beta`` means the inverse temperature (``\beta \equiv 1/T``) and ``\omega_n`` is the Matsubara frequency. Note that ``\omega_n`` is equal to ``(2n + 1) \pi / \beta`` for fermions and ``2n\pi/ \beta`` for bosons (``n`` is an integer).

## [Spectral representation](@id spectrum)

Supposed that the spectral density of the single-particle Green's function is ``A(\omega)``, then we have:
```math
\begin{equation}
G(\tau) = \int^{+\infty}_{-\infty} d\omega
          \frac{e^{-\tau\omega}}{1 \pm e^{-\beta\omega}}
          A(\omega),
\end{equation}
```
with the positive (negative) sign for fermionic (bosonic) operators. Similarly,
```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega
               \frac{1}{i\omega_n - \omega} A(\omega).
\end{equation}
```
The two equations denote the spectral representation of Green's function. We notice that the SPX method, as well as the other analytic continuation methods that are classified as ASM, are closely related to the spectral representation. Next we would like to make further discussions about this representation for the fermionic and bosonic correlators.  

### Fermionic correlators

The spectral density ``A(\omega)`` is defined on ``(-\infty,\infty)``. It is positive definite, i.e., ``A(\omega) \ge 0``. Eq.~(4) and Eq.~(5) can be reformulated as:
```math
\begin{equation}
G(\tau) = \int^{+\infty}_{-\infty} d\omega~K(\tau,\omega) A(\omega),
\end{equation}
```
and
```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega~K(\omega_n,\omega) A(\omega),
\end{equation}
```
respectively. The kernel functions ``K(\tau,\omega)`` and ``K(\omega_n,\omega)`` are defined as follows:
```math
\begin{equation}
K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 + e^{-\beta\omega}},
\end{equation}
```
and
```math
\begin{equation}
K(\omega_n,\omega) = \frac{1}{i\omega_n - \omega}.
\end{equation}
```

### Bosonic correlators

The spectral density ``A(\omega)`` obeys the following constraint: ``\text{sign}(\omega) A(\omega) \ge 0``. Thus, it is more convenient to define a new function ``\tilde{A}(\omega)`` where ``\tilde{A}(\omega) = A(\omega)/\omega``. Clearly, ``\tilde{A}(\omega)`` is always positive definite. As a result Eq.~(4) and Eq.~(5) can be rewritten as:
```math
\begin{equation}
G(\tau) = \int^{+\infty}_{-\infty} d\omega~
    K(\tau,\omega)\tilde{A}(\omega),
\end{equation}
```
and
```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega~
    K(\omega_n,\omega) \tilde{A}(\omega),
\end{equation}
```
respectively. Now the bosonic kernel ``K(\tau,\omega)`` becomes:
```math
\begin{equation}
K(\tau,\omega) = \frac{\omega e^{-\tau\omega}}{1 - e^{-\beta\omega}}.
\end{equation}
```
Especially, ``K(\tau,0) = 1/\beta``. As for ``K(\omega_n,\omega)``, its expression is:
```math
\begin{equation}
K(\omega_n,\omega) = \frac{\omega}{i\omega_n - \omega}.
\end{equation}
```
Especially, ``K(0,0) = -1``. Besides the bosonic Green's function, typical correlator of this kind includes the transverse spin susceptibility ``\chi_{+-}(\tau) = \langle S_{+}(\tau) S_{-}(0) \rangle``, where ``S_{+} = S_x + iS_y`` and ``S_{-} = S_x - i S_y``.

### Bosonic correlators of Hermitian operators

There is a special case of the previous observable kind with ``c = c^{\dagger}``. Here, ``A(\omega)`` becomes an odd function, and equivalently, ``\tilde{A}(\omega)`` is an even function [i.e., ``\tilde{A}(\omega) = \tilde{A}(-\omega)``]. Therefore the limits of integrations in Eq.~(4) and Eq.~(5) are reduced from ``(-\infty,\infty)`` to ``(0,\infty)``. So the two equations can be transformed into:
```math
\begin{equation}
G(\tau) = \int^{+\infty}_{0} d\omega~
    K(\tau,\omega)\tilde{A}(\omega),
\end{equation}
```
and
```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{0} d\omega~
    K(\omega_n,\omega) \tilde{A}(\omega),
\end{equation}
```
respectively. The corresponding ``K(\tau,\omega)`` reads:
```math
\begin{equation}
K(\tau,\omega) = \frac{\omega \left[e^{-\tau\omega} + e^{-(\beta - \tau)\omega}\right]}
                      {1 - e^{-\beta\omega}}.
\end{equation}
```
Especially, ``K(\tau,0) = 2 / \beta``. And ``K(\omega_n,\omega)`` becomes:
```math
\begin{equation}
K(\omega_n, \omega) = \frac{-2\omega^2}{\omega_n^2 + \omega^2}.
\end{equation}
```
Especially, ``K(0,0) = -2``. Perhaps the longitudinal spin susceptibility ``\chi_{zz}(\tau) = \langle S_z(\tau) S_z(0) \rangle`` and the charge susceptibility ``\chi_{ch}(\tau) = \langle N(\tau) N(0) \rangle`` are the most widely used observables of this kind.

As mentioned above, the kernel function is quite strange. The values of ``K(\tau,\omega)`` could change by tens of orders of magnitude. Especially, at large positive and negative frequencies, ``K(\tau,\omega)`` is exponentially small. It implies that at large ``|\omega|`` the features of ``A(\omega)`` depend upon the fine structures of ``G(\tau)``. However, the ``G(\tau)`` data provided by QMC simulations are always fluctuant and noisy. Tiny deviations in ``G(\tau)`` from its expected values can lead to enormous changes in ``A(\omega)``. Thus, analytical continuation is often characterized as an ill-posed problem.

In principle, for incomplete and noise ``G(\tau)`` or ``G(i\omega_n)``, the number of spectral functions ``A(\omega)`` that satisfy Eq.(11) and Eq.(12) is infinite. So the question becomes which ``A(\omega)`` should be chosen. Now there are two different strategies to solve this problem. The first one is to choose the most likely ``A(\omega)``. The second one is to evaluate the average of all the candidate spectral functions.
