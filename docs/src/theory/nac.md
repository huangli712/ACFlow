!!! warning

    The NevanAC solver is experimental. Please use it at your own risk.

!!! info

    This page is written by Dr. Shuang Liang.

## [Nevanlinna Interpolation](@id nac)

It is well known that the retarded Green's function, denoted as ``G^{R}(\omega + i0^{+})``, and the Matsubara Green's function, denoted as ``G(i\omega_n)``, can both be consistently represented as ``G(z)``, where ``z \in \mathbb{C} \backslash \mathbb{R}``. The Nevanlinna analytical continuation (NAC) method utilizes the fact that the negative fermionic Green's function, denoted as ``f(z) = -G(z)``, belongs to the class of Nevanlinna functions. By applying the invertible M``\ddot{\text{o}}``bius transform ``h(z) = (z-i)/(z+i)`` to the function value of ``f(z)``, the Nevanlinna function is mapped in a one-to-one fashion to a contractive function ``\theta(z) = h[f(z)]``. This contractive function ``\theta(z)`` can be expressed in the form of a continued fraction expansion, and an iterative algorithm can be constructed accordingly. The recursion relation between two steps $\theta_j(z)$ and $\theta_{j+1}(z)$ is given by:
```math
\theta_j(z) = \frac{ \theta_{j+1}(z) + \gamma_j}{\gamma_j^*h_j(z) \theta_{j+1}(z) +1}.
```
In this equation, ``h_j(z) = (z-Y_j)/(z+Y_j)``, ``Y_j = i\omega_j`` represents the ``j``-th Matsubara frequency used, and ``\gamma_j = \theta_j(Y_j)`` represents the function value of the ``j``-th contractive function at the point ``Y_j``. The final expression of the recursive function $\theta(z)$ can be written as:
```math
\theta(z)[z;\theta_{N_s+1}(z)] = \frac{a(z)\theta_{N_s+1}(z) + b(z)}{c(z)\theta_{N_s+1}(z) + d(z)},
```
where
```math
  \left(
    \begin{matrix}
      a(z) & b(z) \\
      c(z) & d(z)
    \end{matrix}
  \right) = \prod_{j=1}^{N_s}
  \left(
    \begin{matrix}
      h_j(z)           & \gamma_j \\
      \gamma_j^*h_j(z) & 1
    \end{matrix}
  \right),
```
with ``j`` increasing from left to right. Here ``N_s`` is the overall iteration step, which is equivalent to the number of data points. After obtaining ``\theta(z)``, one can immediately get the Green's function by an inverse ``\textrm{M\"{o}bius}`` transform as ``G(z)= -h^{-1}[\theta(z)]``. Note that the Pick criterion should be fulfilled for the existence of the Nevanlinna interpolation.

## Hardy Basis Functions

Additionally, it is worth noting that there is flexibility in choosing ``\theta_{N_s+1}(z)``, which can be used to select the most desirable spectral function. In original reference for the NAC method, ``\theta_{N_s+1}(z)`` is expanded in the Hardy basis and chosen in such a way that it achieves the smoothest possible spectral function. The loss function employed in this selection process is given by:
```math
\mathcal{L} = \left\vert 1-\int \frac{{\textrm{d}} \omega}{2\pi} \rho_{\theta_{N_s+1}}(\omega) \right\vert^2 + \lambda \left\Vert \frac{{\textrm{d}}^2 \rho_{\theta_{N_s+1}}(\omega)}{{\textrm{d}} \omega^2} \right\Vert^2_F.
```
This loss function consists of two terms. The first term enforces the proper sum rule, while the second term incorporates the smoothness condition. ``\lambda`` is an adjustable parameter. By preserving the ``Nevanlinna'' analytic structure of Green's functions, the NAC method automatically generates positive and normalized spectral functions. However, it is important to emphasize that the method is sensitive to noise, and either a large number of data points $N$ or a high Hardy order $H$ can potentially lead to numerical instabilities.

## Bosonic Cases

Although the NAC method has been extended to support the analytic continuation of matrix-valued Green's functions, it cannot be directly applied to bosonic systems in its original formalism. Quite recently, Nogaki *et al.* suggest an ingenious trick to work around this limitation. Their basic idea is to introduce an auxiliary fermionic function. Let us start with a bosonic Green's function ``G(\tau)`` that satisfies the periodic condition ``G(\tau+\beta) = G(\tau)``. One can construct an artificial anti-periodic fermionic Green's function ``\tilde{G}(\tau)`` as follows:
```math
\tilde{G}(\tau) = \begin{cases}
G(\tau) &\quad (0<\tau<\beta) \\
-G(\tau + \beta) &\quad (-\beta<\tau<0)
\end{cases}
```
Clearly, this auxiliary fermionic Green's function exhibits the same value as the bosonic Green's function in the range ``0<\tau<\beta``. It is easy to prove that the relation between the bosonic spectral function ``\rho(\omega)`` and the auxiliary fermionic spectral function ``\tilde{\rho}(\omega)`` is as follows:
```math
\rho(\omega) = \tilde{\rho}(\omega) \tanh(\beta\omega/2).
```
Furthermore, the sum rule for ``\tilde{\rho}(\omega)`` is given by:
```math
\int_{-\infty}^{\infty} {\textrm{d}} \omega~\tilde{\rho}(\omega) = G(\tau=0^+) + G(\tau=\beta-0^-).
```
Given ``\tilde{G}(\tau)``, it is easy to construct ``\tilde{G}(i\nu_n)`` via direct Fourier transformation, where ``\nu_n = (2n+1)\pi/\beta`` are the fermionic Matsubara frequencies. Since
```math
\tilde{G}(i\nu_n) = \int_{-\infty}^{\infty}
    {\textrm{d}} \omega \frac{\tilde{\rho}(\omega)}{i\nu_n - \omega},
```
one can perform analytic continuation for ``\tilde{G}(i\nu_n)`` via the standard NAC method to get ``\tilde{\rho}(\omega)``. And then the bosonic spectral function ``\rho(\omega)`` can be derived according to the above equations.

## Relevant Parameters

See [[NevanAC] Block](@ref nevanac_block)

## References

**[1]** Jiani Fei, Chia-Nan Yeh, and Emanuel Gull, Nevanlinna Analytical Continuation, *Phys. Rev. Lett.* **126**, 056402 (2021).
