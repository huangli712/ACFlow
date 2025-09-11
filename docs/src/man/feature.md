Now the ACFlow toolkit supports six analytic continuation methods as introduced before. It includes seven different analytic continuation solvers, namely

* `MaxEnt`
* `BarRat`
* `NevanAC`
* `StochAC`
* `StochSK`
* `StochOM`
* `StochPX`

Just as their names suggested, the `MaxEnt` solver implements the maximum entropy method. The `BarRat` solver implements the barycentric rational function approximation. The `NevanAC` solver implements the Nevanlinna analytical continuation. The `StochAC` and `StochSK` solvers implement the K. S. D. Beach's algorithm and A. W. Sandvik's algorithm of the stochastic analytic continuation, respectively. The `StochOM` solver implements the stochastic optimization method. The `StochPX` solver implements the stochastic pole expansion method. The ACFlow toolkit also provides a convenient library, which can be used to prepare and carry out analytic continuation calculations flexibly. The major features of the present ACFlow toolkit (*v2.0.0* and above) are summarized in **Table 1**.

| Features | MaxEnt | BarRat | NevanAC | StochAC | StochSK | StochOM | StochPX |
| :------- | :----: | :----: | :-----: | :-----: | :-----: | :-----: | :-----: |
|Matrix-valued Green's function | Y | Y | N | N | N | N | Y |
|Fragment input grid            | Y | N | Y | Y | Y | Y | Y |
|Imaginary time grid            | Y | N | N | Y | Y | Y | N |
|Matsubara frequency grid       | Y | Y | Y | Y | Y | Y | Y |
|Linear mesh                    | Y | Y | Y | Y | Y | Y | Y |
|Nonlinear mesh                 | Y | Y | Y | Y | Y | Y | Y |
|Fermionic kernel               | Y | Y | Y | Y | Y | Y | Y |
|Bosonic kernel                 | Y | Y | N | Y | Y | Y | Y |
|Self-defined model function    | Y | N | N | N | N | N | N |
|Constrained analytic continuation | N | N | N | Y | Y | Y | Y |
|Self-adaptive parameterization | N | N | N | Y | N | N | Y |
|Regeneration of input data     | Y | Y | Y | Y | Y | Y | Y |
|Kramers-Kronig transformation  | Y | Y | Y | Y | Y | Y | Y |
|Parallel computing             | N | N | N | Y | Y | Y | Y |
|Parallel tempering             | N | N | N | Y | N | N | N |
|Interactive mode               | Y | Y | Y | Y | Y | Y | Y |
|Script mode                    | Y | Y | Y | Y | Y | Y | Y |
|Standard mode                  | Y | Y | Y | Y | Y | Y | Y |

**Table 1 |** Major features of the ACFlow toolkit. `MaxEnt`, `BarRat`, `NevanAC`, `StochAC`, `StochSK`, `StochOM`, and `StochPX` are the seven analytic continuation solvers as implemented in this toolkit.

In **Table 1**, `Y` means yes while `N` means no. `Interactive mode`, `Script mode`, and `Standard model` are three running modes supported by the ACFlow toolkit. We will introduce them later. The `MaxEnt` solver supports the `historic`, `classic`, `bryan`, and `chi2kink` algorithms to determine the ``\alpha`` parameter. The `StochAC` solver is only compatible with a flat model function, while the `BarRat`, `NevanAC`, `StochSK`, `StochOM`, and `StochPX` solvers don't rely on any default model functions. For the `StochAC` and `StochPX` solvers, they can refine the nonlinear mesh for parameterization of ``\delta``-like peaks or poles, according to the spectral function generated in any previous runs or by some `a prioir` information. This is the so-called `self-adaptive parameteration`.

!!! info

    Note that analytic continuation problem is a hotspot in computational physics and many-body physics all the time. Many efforts have been devoted to solve it in recent years. Noticeable achievements include maximum quantum entropy method, Nevanlinna analytical continuation, blocked-mode sampling and grid point sampling in stochastic analytic continuation, constrained stochastic analytic continuation, machine learning assisted analytic continuation, and so on. We would like to incorporate these new progresses into the ACFlow toolkit in the near future. BTW, contributions from the other users are always welcomed.
