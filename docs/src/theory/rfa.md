!!! warning

    The BarRat solver is experimental. Please use it at your own risk.

## [Barycentric Rational Function](@id rfa)

The barycentric formula takes the form of a quotient of two partial fractions,

```math
r(z) = \frac{n(z)}{d(z)}
     = \sum^m_{j=1} \frac{w_j f_j}{z - z_j}
     {\huge/} \sum^m_{j=1} \frac{w_j}{z - z_j},
```

where ``m \ge 1`` is an integer, ``z_1, \cdots, z_m`` are a set of real or complex distinct support points (`nodes`), ``f_1, \cdots, f_m`` are a set of real or complex data `values`, and ``w_1, \cdots, w_m`` are a set of real or complex `weights`. As indicated in this equation, we just let ``n(z)`` and ``d(z)`` stand for the partial fractions in the numerator and the denominator.

## Prony's Interpolation

Our input data consists of an odd number ``2N + 1`` of Matsubara points ``G(i\omega_n)`` that are uniformly spaced. Prony's interpolation method interpolates ``G_k`` as a sum of exponentials

```math
G_k = \sum^{N-1}_{i=0} w_i \gamma^k_i,
```

where ``0 \le k \le 2N``, ``w_i`` denote complex weights and ``\gamma_i`` corresponding nodes.

## Prony's Approximation

Prony's interpolation method is unstable. We therefore employs a Prony approximation, rather than an interpolation of ``G``. For the physical Matsubara functions, which decay in magnitude to zero for ``i\omega_n \to i\infty``, only ``K \propto \log{1/\varepsilon}`` out of
all ``N`` nodes in the Prony approximation have weights ``|w_i| > \varepsilon``. Thus, we have

```math
\left|G_k - \sum^{K-1}_{i=0} w_i \gamma^k_i\right| \le \varepsilon,
```

for all ``0 \le k \le 2N``.

## Relevant Parameters

See [[BarRat] Block](@ref barrat_block)

## References

**[1]** Li Huang, Changming Yue, Barycentric rational function approximation made simple: A fast analytic continuation method for Matsubara Green's functions, *Phys. Rev. B* **111**, 125139 (2025).
