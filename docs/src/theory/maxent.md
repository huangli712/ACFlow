Perhaps the maximum entropy method is the most frequently used approach for analytic continuation problems because of its high computational efficiency. Next, we will discuss the basic principles and several variants of it.

## [Bayesian Inference](@id mem)

Bayes's theorem is the cornerstone of the maximum entropy method. Given two events ``a`` and ``b``, Bayes's theorem says:
```math
P[a|b]P[b] = P[b|a]P[a],
```
where ``P[a]`` is the probability of event ``a``, ``P[a|b]`` is the conditional probability of event ``a`` with given event ``b``. In the scenario of analytic continuation problem, ``\bar{G}(\tau)`` and ``A(\omega)`` are treated as two events, where ``\bar{G}(\tau)`` denotes the measured value of ``G(\tau)``. So the best solution for ``A(\omega)`` is of course the one that maximizes ``P[A|\bar{G}]``, which is called the posterior probability. According to the Bayes's theorem, we get
```math
P[A|\bar{G}] = \frac{P[\bar{G}|A]P[A]}{P[\bar{G}]},
```
where ``P[\bar{G}|A]`` is the likelihood function, ``P[A]`` is the prior probability, and ``P[\bar{G}]`` is the evidence. Since the evidence is a normalization constant depending on the prior probability and the likelihood function only, it is ignored in the following discussions. Thus,
```math
P[A|\bar{G}] \propto P[\bar{G}|A]P[A].
```

## Posterior Probability

In the maximum entropy method, the likelihood function ``P[\bar{G}|A]`` is assumed to be in direct proportion to ``e^{-\chi^2/2}``. ``\chi^2`` is named as goodness-of-fit function, which measures distance between ``\bar{G}(\tau)`` and reconstructed imaginary time Green's function ``\tilde{G}(\tau)``:
```math
\chi^2 = \sum^{L}_{i = 1} \left[\frac{\bar{G}_i(\tau) - \tilde{G}_i(\tau)}{\sigma_i}\right]^2,
```
```math
\tilde{G}_i = \sum_j K_{ij} A_j.
```
Here, ``L`` is number of imaginary time points, ``\sigma`` denotes the error bar (standard deviation) of ``\bar{G}(\tau)``. ``K_{ij}`` and ``A_j`` are discrete kernel and spectral functions, respectively. On the other hand, the prior probability ``P[A]`` is supposed to be in direct proportion to ``e^{\alpha S}``, where ``\alpha`` is a regulation parameter and ``S`` means entropy. Sometimes ``S`` is also known as the Kullback-Leibler distance (or the Shannon-Jaynes entropy). Its formula is as follows:
```math
S= \int d\omega \left(A(\omega) - m(\omega) - A(\omega)\log\left[\frac{A(\omega)}{m(\omega)}\right]\right),
```
where ``m(\omega)`` is the default model function. The ACFlow toolkit also supports another kind of entropy, i.e., the Bayesian Reconstruction entropy. It reads:
```math
S = \int d\omega \left(1 - \frac{A(\omega)}{m(\omega)} + \log \left[\frac{A(\omega)}{m(\omega)}\right]\right).
```

According to the Bayes's theorem, the posterior probability ``P[A|\bar{G}] \propto e^{Q}`` and
```math
Q = \alpha S - \frac{\chi^2}{2}.
```

## Various Algorithms

Now the original analytic continuation problem becomes how to figure out the optimal ``A(\omega)`` that maximizes ``Q``. In other words, we have to solve the following equation:
```math
\frac{\partial Q}{\partial A} \bigg|_{A = \hat{A}} = 0,
```
where ``\hat{A}(\omega)`` is the optimal ``A(\omega)``. This equation can be easily solved by using standard Newton method. However, the obtained ``\hat{A}(\omega)`` is ``\alpha``-dependent. That is to say, for a given ``\alpha``, there is always a ``\hat{A}(\omega)`` that satisfies this equation. So, new problem arises because we have to find out a way to generate the final spectral function from these ``\alpha``-resolved ``\hat{A}(\omega)``. Now there exist four algorithms, namely `historic`, `classic`, `bryan`, and `chi2kink`. Next we will introduce them one by one.

### Historic Algorithm

The historic algorithm is quite simple. The ``\alpha`` parameter will be adjusted iteratively to meet the following criterion:
```math
\chi^2 = N,
```
where ``N`` is the number of mesh points for spectral density ``A(\omega)``.

### Classic Algorithm

The basic equation for the classic algorithm reads
```math
-2 \alpha S(A_{\alpha}) = \text{Tr}
\left[
\frac{\Lambda(A_{\alpha})}{\alpha I + \Lambda(A_{\alpha})}
\right],
```
where ``I`` is an identity matrix. The elements of ``\Lambda`` matrix are calculated as follows:
```math
\Lambda_{ij} = \sqrt{A_i} \left(\sum_{kl} K_{ki} [C^{-1}]_{kl} K_{lj}\right) \sqrt{A_j},
```
where ``C`` is the covariance matrix. The basic equation will be iteratively solved until the optimal ``\alpha`` and ``\hat{A}(\omega)`` are figured out.

### Bryan Algorithm

In both historic and classic algorithms, the spectral function ``\hat{A}(\omega)`` is always related to an optimal ``\alpha`` parameter. However, the spirit of the bryan algorithm is completely different. It tries to yield a series of ``\alpha`` parameters and calculate the corresponding ``A_{\alpha}(\omega)``. Then the final spectral function ``A(\omega)`` is obtained by evaluating the following integration:
```math
\overline{A(\omega)} = \int d\alpha~A_{\alpha}(\omega) P[\alpha | \bar{G}].
```

### Chi2kink Algorithm

This algorithm was proposed by Bergeron and Tremblay recently. The first step is to generate a series of ``\alpha`` parameters and evaluate the corresponding spectral functions ``A_{\alpha}(\omega)`` and goodness-of-fit functions ``\chi^{2}[A_{\alpha}]``. Then we plot ``\log_{10}(\chi^{2})`` as a function of ``\log_{10}(\alpha)``. Usually this plot is split into three different regions: (1) *Default model region*. In the limit of ``\alpha \to \infty``, ``\chi^{2}`` goes to a constant high value. It means that the likelihood function ``e^{-\chi^2/2}`` has negligible weight, such that the prior probability ``e^{\alpha S}`` becomes dominant and minimizes ``Q[A]``. At that time, the calculated ``A(\omega)`` resembles the default model function ``m(\omega)``. (2) *Noise-fitting region*. In the limit of ``\alpha \to 0``, ``\chi^2`` is relatively flat and approaches its global minimum. In this region, the minimization algorithm tends to fit the noise in ``G(\tau)``. (3) *Information-fitting region*. ``\alpha S`` is comparable with ``\chi^2/2``, so that ``\chi^{2}`` is strongly dependent on ``\alpha``. Bergeron *et al.* suggested that the optimal ``\alpha`` parameter situates in the crossover between noise-fitting region and information-fitting region. So the second derivative of ``\chi^{2}`` with respect to ``\alpha`` is calculated, and the maximum value in the resulting curve indicates the optimal value of ``\alpha``. Quite recently, Kaufmann and Held proposed a more numerically stable and flexible approach to compute the optimal ``\alpha``. They use the following function to fit ``\log_{10}[\chi^{2}(\alpha)]-\log_{10}(\alpha)``:
```math
\phi(x;a,b,c,d) = a + \frac{b}{1 + e^{-d(x-c)}},
```
where ``a``, ``b``, ``c``, and ``d`` are fitting parameters. Then the optimal ``\alpha`` is approximated by ``10^{c-f/d}``, where ``f`` is an empirical constant (The favorite value of ``f`` lies in ``[2,2.5]``).

## Relevant Parameters

See [[MaxEnt] Block](@ref maxent_block)

## References

**[1]** M. Jarrell and J. E. Gubernatis, Bayesian inference and the analytic continuation of imaginary-time quantum Monte Carlo data, *Phys. Rep.* **269**, 133 (1996).

**[2]** D. Bergeron and A.-M. S. Tremblay, Algorithms for optimized maximum entropy and diagnostic tools for analytic continuation, *Phys. Rev. E* **94**, 023303 (2016).
