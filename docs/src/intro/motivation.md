## Available Codes

In recent years, quite a few analytic continuation codes have been released, including

* maxent (by Mark Jarrell)
* ``\Omega``Maxent
* ana\_cont
* ALPSCore/maxent
* TRIQS/som
* ALF
* Nevanlinna.jl
* PronyAC
* Minipole

just to name a few. We note that the maximum entropy method has dominated this field for quite a long time. Thus most of these codes only support the maximum entropy method. It is rather difficult to crosscheck the simulated results obtained by various analytic continuation methods. In addition, the features of the available codes are quite limited and hard to be extended.

## New Toolkit: ACFlow

In order to fill in this gap, we would like to present a new open source toolkit, called **ACFlow**, for analytic continuation. This toolkit implements three primary analytic continuation methods, namely the

* Maximum entropy method
* Stochastic analytic continuation
* Stochastic optimization method

within an united framework. It provides an easy-to-used library and application interface. Some diagnostic and analytic tools are also available. With ACFlow, the users can easily setup and execute analytic continuation calculations, and validate the calculated results. We believe that this toolkit will play a vital role in solving analytic continuation problems.

!!! info

    Quite recently, three new analytic continuation methods, namely

    * Nevanlinna analytical continuation
    * Stochastic pole expansion
    * Barycentric rational function approximation

    have been implemented in the ACFlow toolkit. So, now it supports six analytic continuation methods :-).
