# Main Features

Now the ACFlow toolkit supports four analytical continuation methods as introduced before. It includes five different analytical continuation solvers, namely 

* `MaxEnt`
* `StochAC`
* `StochSK`
* `StochOM`
* `StochPX`

Just as their names suggested, the `MaxEnt` solver implements the maximum entropy method. The `StochAC` and `StochSK` solvers implement the K. S. D. Beach's algorithm and A. W. Sandvik's algorithm of the stochastic analytical continuation, respectively. The `StochOM` solver implements the stochastic optimization method. The `StochPX` solver implements the stochastic pole expansion method. The ACFlow toolkit also provides a convenient library, which can be used to prepare and carry out analytical continuation calculations flexibly. The major features of the present ACFlow toolkit (*v1.2.0*) are summarized in **Table 1**.

| Features | MaxEnt | StochAC | StochSK | StochOM | StochPX |
| :------- | :----: | :-----: | :-----: | :-----: | :-----: |
|Matrix-valued Green's function | Y | N | N | N | N |
|Imaginary time grid            | Y | Y | Y | Y | N |
|Matsubara frequency grid       | Y | Y | Y | Y | Y |
|Linear mesh                    | Y | Y | Y | Y | Y |
|Nonlinear mesh                 | Y | Y | Y | Y | Y |
|Fermionic kernel               | Y | Y | Y | Y | Y |
|Bosonic kernel                 | Y | Y | Y | Y | N |
|Self-defined model function    | Y | N | N | N | N |
|Constrained analytical continuation | N | Y | Y | Y | Y |
|Regeneration of input data     | Y | Y | Y | Y | Y |
|Kramers-Kronig transformation  | Y | Y | Y | Y | Y |
|Parallel computing             | N | Y | Y | Y | Y |
|Parallel tempering             | N | Y | N | N | N |
|Interactive mode               | Y | Y | Y | Y | Y |
|Script mode                    | Y | Y | Y | Y | Y |
|Standard mode                  | Y | Y | Y | Y | Y |

**Table 1 |** Major features of the ACFlow toolkit. `MaxEnt`, `StochAC`, `StochSK`, `StochOM`, and `StochPX` are the five analytical continuation solvers as implemented in this toolkit.

In **Table 1**, `Y` means yes while `N` means no. `Interactive mode`, `Script mode`, and `Standard model` are three running modes supported by the ACFlow toolkit. We will introduce them later. The `MaxEnt` solver supports the `historic`, `classic`, `bryan`, and `chi2kink` algorithms to determine the ``\alpha`` parameter. The `StochAC` solver is only compatible with a flat model function, while the `StochSK`, `StochOM`, and `StochPX` solvers don't rely on any default model functions. The `StochOM` solver does not support analytical continuation of fermionic imaginary time Green's function for the moment. 

!!! info

    Note that analytical continuation problem is a hotspot in computational physics and many-body physics all the time. Many efforts have been devoted to solve it in recent years. Noticeable achievements include maximum quantum entropy method, Nevanlinna analytical continuation, blocked-mode sampling and grid point sampling in stochastic analytical continuation, constrained stochastic analytical continuation, machine learning assisted analytical continuation, and so on. We would like to incorporate these new progresses into the ACFlow toolkit in the near future. BTW, contriibutions from the users are always welcomed.