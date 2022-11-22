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

## [Spectral Density](@id spectrum)

Clearly, neither ``G(\tau)`` nor ``G(i\omega_n)`` can be observed experimentally. We have to extract dynamical response function, i.e., the spectral density ``A(\omega)``, from them. ``A(\omega)`` is indeed an observable quantity. It is related to ``G(\tau)`` via the following Laplace transformation:
```math
\begin{equation}
G(\tau) = \int^{+\infty}_{-\infty} d\omega \frac{e^{-\tau\omega}}{1 \pm e^{-\beta\omega}} A(\omega),
\end{equation}
```
where +(-) in the denominator is for fermionic (bosonic) system. ``G(i\omega_n)`` and ``A(\omega)`` manifest similar relation:
```math
\begin{equation}
G(i\omega_n) = \int^{+\infty}_{-\infty} d\omega' \frac{A(\omega')}{i\omega_n - \omega'}.
\end{equation}
```
It is obvious that Eq.~(\ref{eq:spectral_density_1}) and Eq.~(\ref{eq:spectral_density_2}) are indeed special forms of the Fredholm integral equation of the first kind [see Eq.~(\ref{eq:fredholm})]. So, the central problem of analytical continuation is to search optimal ``A(\omega)`` for given ``G(\tau)`` or ``G(i\omega_n)``.

Sometimes the spectral density ``A(\omega)`` is called as spectral function in the references. It is tied to the imaginary part of real frequency Green's function ``G(\omega)``:
```math
\begin{equation}
A(\omega) = -\frac{1}{\pi} \rm{Im}G(\omega).
\end{equation}
```
From Im``G(\omega)``, Re``G(\omega)`` could be calculated via the Kramers-Kronig transformation:
```math
\begin{equation}
\mathrm{Re} G(\omega) = \frac{1}{\pi} \mathcal{P}
  \int_{-\infty}^{\infty} d\omega'~
  \frac{\mathrm{Im} G(\omega')}{\omega'-\omega},
\end{equation}
```
where ``\mathcal{P}`` means Cauchy principal value. Besides Eq.~(\ref{eq:spectral_density_1}) and Eq.~(\ref{eq:spectral_density_2}), ``A(\omega)`` has to obey some additional constraints or sum-rules. For fermionic systems, the spectral functions must be positive:
```math
\begin{equation}
A(\omega) \ge 0.
\end{equation}
```
While for bosonic systems, the above condition turns into:
```math
\begin{equation}
\text{sign}(\omega) A(\omega) \ge 0.
\end{equation}
```
In addition, the spectral function ``A(\omega)`` is always bounded,
```math
\begin{equation}
\int^{+\infty}_{-\infty} d\omega~A(\omega) < \infty.
\end{equation}
```
It can be utilized to normalize the final spectral function.

## Kernel functions

Eq.~(\ref{eq:spectral_density_1}) and Eq.~(\ref{eq:spectral_density_2}) can be reformulated as follows:
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
where ``K(\tau,\omega)`` and ``K(\omega_n, \omega')`` are the so-called kernel functions. Their definitions are as follows:
```math
\begin{equation}
K(\tau,\omega) = \frac{e^{-\tau\omega}}{1 \pm e^{-\beta\omega}},
\end{equation}
```
and
```math
\begin{equation}
K(\omega_n,\omega) = \frac{1}{i\omega_n - \omega},
\end{equation}
```
where +(-) in the denominator of Eq.~(\ref{eq:ktau}) stands for fermions (bosons).

As mentioned above, the kernel function is quite strange. The values of ``K(\tau,\omega)`` could change by tens of orders of magnitude. Especially, at large positive and negative frequencies, ``K(\tau,\omega)`` is exponentially small. It implies that at large ``|\omega|`` the features of ``A(\omega)`` depend upon the fine structures of ``G(\tau)``. However, the ``G(\tau)`` data provided by QMC simulations are always fluctuant and noisy. Tiny deviations in ``G(\tau)`` from its expected values can lead to huge changes in ``A(\omega)``. Thus, analytical continuation is often characterized as an ill-posed problem.

In principle, for incomplete and noise ``G(\tau)`` or ``G(i\omega_n)``, the number of spectral functions ``A(\omega)`` that satisfy Eq.~(\ref{eq:kernel_t}) and Eq.~(\ref{eq:kernel_w}) is infinite. So the question becomes which ``A(\omega)`` should be chosen. Now there are two different strategies to solve this problem. The first one is to choose the most likely ``A(\omega)``. The second one is to evaluate the average of all the candidate spectral functions. In next section, we will introduce three primary analytical continuation methods that follow the two strategies and have been implemented in the ACFlow toolkit. We will concentrate on analytical continuation of imaginary time Green's functions in main text.