## Analytic Continuation Problem

It is well-known that quantum Monte Carlo (QMC) method is a powerful and exact numerical approach, and has been widely used in many research fields, such as nuclear physics, condense matter physics, and many-body physics. Here, we just focus on the finite temperature QMC algorithms, which are used to solve interacting lattice models or quantum impurity models. Generally speaking, the simulated results of QMC methods are some sorts of single-particle or two-particle correlation functions, which are usually defined on imaginary time axis (``\tau \equiv -it``) or Matsubara frequency axis (``i\omega_n``). Therefore, they can't be compared directly with the correspondingly experimental results, including but not limited to electronic density of states ``A(\omega)``, optical conductivity ``\sigma(\omega)``, dynamical structure factor ``S(\mathbf{q},\omega)``, and so on. It is necessary to convert the QMC simulated results from imaginary time axis or Matsubara frequency axis to real axis (i.e. ``\tau \to \omega`` or ``i\omega_n \to \omega``), which is the origin of the analytic continuation problem.

## Fredholm Integral Equation

Let's concentrate on the following Fredholm integral equation of the first kind:
```math
g(y) = \int K(y,x) f(x)~\text{d}x.
```
Here, $K(y,x)$ is the known kernel function, ``f(x)`` is the model function, and ``g(y)`` denotes raw data. Given ``f(x)``, it is quite easy to get ``g(y)`` via numerical integration. However, given ``g(y)``, solving the Fredholm integral equation reversely to get ``f(x)`` is not as easy as expected. There is no universal solution. In some cases, even the existence of solution can not be guaranteed.

## Analytic Continuation Methods

The so-called analytic continuation problem can be reformulated in terms of the Fredholm integral equation. Thus, its objective is to seek a reasonable ``f(x)`` to satisfy the above equation. The QMC simulated data ``g(y)`` are noisy and the kernel function ``K(y,x)`` is ill conditioned, which make analytic continuation of QMC simulated data a huge challenge. In order to solve this problem, in the past decades peoples have developed numerous methods, including the least square fitting method, singular value decomposition, Pad``\text{\'{e}}`` approximation, Tikhonov-Philips regularization method, maximum entropy method, stochastic analytic continuation, stochastic optimization method, sparse modelling method, and machine learning method, etc. However, each method has its pros and cons. None of these methods can override the others. The analytic continuation problem is still far away from being completely solved.
