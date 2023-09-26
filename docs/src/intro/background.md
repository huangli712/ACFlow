# Background

## Analytic Continuation Problem

It is well-known that quantum Monte Carlo (QMC) method is a powerful and exact numerical approach, and has been widely used in many research fields, such as nuclear physics, condense matter physics, and many-body physics. Here, we just focus on the finite temperature QMC algorithms, which are used to solve interacting lattice models or quantum impurity models. Generally speaking, the simulated results of QMC methods are some sorts of single-particle or two-particle correlation functions, which are usually defined on imaginary time axis (``\tau \equiv -it``) or Matsubara frequency axis (``i\omega_n``). Therefore, they can't be compared directly with the correspondingly experimental results, including but not limited to electronic density of states ``A(\omega)``, optical conductivity ``\sigma(\omega)``, dynamical structure factor ``S(\mathbf{q},\omega)``, and so on. It is necessary to convert the QMC simulated results from imaginary time axis or Matsubara frequency axis to real axis (i.e. ``\tau \to \omega`` or ``i\omega_n \to \omega``), which is the origin of the analytic continuation problem.

## Fredholm Integral Equation

Let's concentrate on the following Fredholm integral equation of the first kind:
```math
g(y) = \int K(y,x) f(x)~dx.
```
Here, $K(y,x)$ is the known kernel function, ``f(x)`` is the model function, and ``g(y)`` denotes raw data. Given ``f(x)``, it is quite easy to get ``g(y)`` via numerical integration. However, given ``g(y)``, solving the Fredholm integral equation reversely to get ``f(x)`` is not as easy as expected. There is no universal solution. In some cases, even the existence of solution can not be guaranteed.

## Available Analytic Continuation Methods

The so-called analytic continuation problem can be reformulated in terms of the Fredholm integral equation. Thus, its objective is to seek a reasonable ``f(x)`` to satisfy the above equation. The QMC simulated data ``g(y)`` are noisy and the kernel function ``K(y,x)`` is ill conditioned, which make analytic continuation of QMC simulated data a huge challenge. In order to solve this problem, in the past decades peoples have developed numerous methods, including the least square fitting method, singular value decomposition, Pad``\text{\'{e}}`` approximation, Tikhonov-Philips regularization method, maximum entropy method, stochastic analytic continuation, stochastic optimization method, sparse modelling method, and machine learning method, etc. However, each method has its pros and cons. None of these methods can override the others. The analytic continuation problem is still far away from being completely solved.

## Our Motivations

In recent years, quite a few analytic continuation codes have been released, including maxent (by Mark Jarrell), ``\Omega``Maxent, ana\_cont, ALPSCore/maxent, TRIQS/som, ALF, just to name a few. We note that the maximum entropy method has dominated this field for quite a long time. Thus most of these codes only support the maximum entropy method. It is rather difficult to crosscheck the simulated results obtained by various analytic continuation methods. In addition, the features of the available codes are quite limited and hard to be extended. In order to fill in this gap, we would like to present a new open source toolkit, called ACFlow, for analytic continuation. This toolkit implements three primary analytic continuation methods, namely the maximum entropy method, stochastic analytic continuation, and stochastic optimization method, within an united framework. It provides an easy-to-used library and application interface. Some diagnostic and analytic tools are also available. With ACFlow, the users can easily setup and execute analytic continuation calculations, and validate the calculated results. We believe that this toolkit will play a vital role in solving analytic continuation problems.

!!! info

    Quite recently, a new analytic continuation method, namely the stochastic pole expansion, has been implemented in the ACFlow toolkit. So, now it supports four analytic continuation methods :-).
