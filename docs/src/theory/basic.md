# Basic Principles

## Finite temperature Green's functions

Under the Wick's rotation $t \to i\tau$, the time evolution operator in the Heisenberg picture $e^{itH}$ will be replaced by $e^{-\tau H}$. Such a transformation will increase efficiency of QMC random walking and suppress numerical oscillation (when $t$ is large, the periodic oscillation of $e^{itH}$ is quite obvious). This is an important reason why most of the finite temperature QMC algorithms are formulated in imaginary time axis. The outputs of finite temperature QMC simulations are usually single-particle or two-particle correlation functions. For example, the imaginary time Green's function $G(\tau)$ is defined as follows: 
```math
\begin{equation}
G(\tau) = \langle \mathcal{T}_{\tau} d(\tau) d^{\dagger}(0) \rangle,
\end{equation}
```
where $\tau$ denotes imaginary time, $\mathcal{T}_{\tau}$ denotes time-ordered operator, and $d$ and $d^{\dagger}$ are annihilation and creation operators, respectively. The Matsubara Green's function $G(i\omega_n)$ can be measured by QMC simulations or constructed from $G(\tau)$ via direct Fourier transformation:
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
Here, $\beta$ means the inverse temperature ($\beta \equiv 1/T$) and $\omega_n$ is the Matsubara frequency. Note that $\omega_n$ is equal to $(2n + 1) \pi / \beta$ for fermions and $2n\pi/ \beta$ for bosons ($n$ is an integer).
