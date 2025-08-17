#
# Project : Gardenia
# Source  : maxent.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/08/18
#

#=
### *Customized Structs* : *MaxEnt Solver*
=#

"""
    MaxEntContext

Mutable struct. It is used within the MaxEnt solver only.

### Members
* Gᵥ     -> Input data for correlator.
* σ²     -> Actually 1.0 / σ².
* grid   -> Grid for input data.
* mesh   -> Mesh for output spectrum.
* model  -> Default model function.
* kernel -> Default kernel function.
* Vₛ     -> Matrix from singular value decomposition.
* W₂     -> Precomputed array.
* W₃     -> Precomputed array.
* Bₘ     -> Precomputed array.
"""
mutable struct MaxEntContext
    Gᵥ     :: Vector{F64}
    σ²     :: Vector{F64}
    grid   :: AbstractGrid
    mesh   :: AbstractMesh
    model  :: Vector{F64}
    kernel :: Array{F64,2}
    hess   :: Array{F64,2}
    Vₛ     :: Array{F64,2}
    W₂     :: Array{F64,2}
    W₃     :: Array{F64,3}
    Bₘ     :: Vector{F64}
end

#=
### *Global Drivers*
=#

"""
    solve(S::MaxEntSolver, rd::RawData)

Solve the analytic continuation problem by the maximum entropy method. It
is the driver for the MaxEnt solver.

If the input correlators are bosonic, this solver will return A(ω) / ω
via `Aout`, instead of A(ω). At this time, `Aout` is not compatible with
`Gout`. If the input correlators are fermionic, this solver will return
A(ω) in `Aout`. Now it is compatible with `Gout`. These behaviors are just
similar to the StochAC, StochSK, and StochOM solvers.

It seems that the MaxEnt solver is hard to create δ-like spectra.

### Arguments
* S -> A MaxEntSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* mesh -> Real frequency mesh, ω.
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
"""
function solve(S::MaxEntSolver, rd::RawData)
    println("[ MaxEnt ]")
    #
    mec = init(S, rd)
    darr, sol = run(mec)
    gout = last(mec, darr, sol)
    #
    return mec.mesh.mesh, sol[:A], gout
end

"""
    init(S::MaxEntSolver, rd::RawData)

Initialize the MaxEnt solver and return a MaxEntContext struct.

### Arguments
* S -> A MaxEntSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* mec -> A MaxEntContext struct.
"""
function init(S::MaxEntSolver, rd::RawData)
    # Prepera input data
    G = make_data(rd)
    Gᵥ = G.value
    σ² = 1.0 ./ G.covar
    println("Postprocess input data: ", length(σ²), " points")

    # Prepare grid for input data
    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    # Prepare mesh for output spectrum
    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    # Prepare default model function
    model = make_model(mesh)
    println("Build default model: ", get_b("mtype"))

    # Prepare kernel function
    kernel = make_kernel(mesh, grid)
    println("Build default kernel: ", get_b("ktype"))

    # Prepare some essential intermediate variables
    Vₛ, W₂, W₃, Bₘ, hess = precompute(Gᵥ, σ², mesh, model, kernel)
    println("Precompute key coefficients")

    return MaxEntContext(Gᵥ, σ², grid, mesh, model,
                         kernel, hess, Vₛ, W₂, W₃, Bₘ)
end

"""
    run(mec::MaxEntContext)

Perform maximum entropy simulation with different algorithms. Now it
supports the `historic`, `classic`, `bryan`, and `chi2kink` algorithms.

### Arguments
* mec -> A MaxEntContext struct.

### Returns
* svec -> A vector of dictionaries. It contains the intermediate solutions.
* sol -> Dictionary. It contains the final solution.
"""
function run(mec::MaxEntContext)
    stype = get_m("stype")
    method = get_m("method")

    # Note that the Bayesian Reconstruction entropy is compatible with
    # all the four algorithms so far.
    if stype == "br"
        println("Bayesian Reconstruction entropy is used!")
    else
        println("Shannon–Jaynes entropy is used!")
    end

    @cswitch method begin
        @case "historic"
            return historic(mec)
            break

        @case "classic"
            return classic(mec)
            break

        @case "bryan"
            return bryan(mec)
            break

        @case "chi2kink"
            return chi2kink(mec)
            break
    end
end

"""
    last(mec::MaxEntContext, svec::Vector, sol::Dict)

Postprocess the results generated during the maximum entropy simulations.
Here `sol` is the final solution for the analytic continuation problem,
while `svec` contains all the intermediate results (it is a vector of
dictionary actually).

### Arguments
* mec -> A MaxEntContext struct.
* svec -> See above explanations.
* sol -> See above explanations.

### Returns
* G -> Retarded Green's function, G(ω).
"""
function last(mec::MaxEntContext, svec::Vector, sol::Dict)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Write the spectral function
    fwrite && write_spectrum(mec.mesh, sol[:A])

    # Write the model function
    fwrite && write_model(mec.mesh, mec.model)

    # Write α-χ² data
    α_vec = map(x -> x[:α], svec)
    χ_vec = map(x -> x[:χ²], svec)
    fwrite && write_misfit(α_vec, χ_vec)

    # Write P[α|A] for bryan algorithm
    if haskey(svec[end], :prob)
        p_vec = map(x -> x[:prob], svec)
        fwrite && write_probability(α_vec, p_vec)
    end

    # Regenerate the input data and write them
    #
    # Perhaps the spectral function (sol[:A]) is blurred.
    # But we need the pure spectral function (sol[:Araw]).
    Aout = haskey(sol, :Araw) ? sol[:Araw] : sol[:A]
    G = reprod(mec.mesh, mec.kernel, Aout)
    fwrite && write_backward(mec.grid, G)

    # Calculate full response function on real axis and write them
    if get_b("ktype") == "fermi"
        _G = kramers(mec.mesh, Aout)
    else
        _G = kramers(mec.mesh, Aout .* mec.mesh)
    end
    fwrite && write_complete(mec.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

"""
    historic(mec::MaxEntContext)

Apply the historic algorithm to solve the analytic continuation problem.
It choose α in a way that χ² ≈ N.

For the historic algorithm, `alpha` is usually 10⁶, and `ratio` is 10.0.
It is compatible with the Bayesian Reconstruction entropy.

### Arguments
* mec -> A MaxEntContext struct.

### Returns
* svec -> A vector of dictionaries. It contains the intermediate solutions.
* sol -> Dictionary. It contains the final solution.

See also: [`MaxEntContext`](@ref).
"""
function historic(mec::MaxEntContext)
    function root_fun(_α, _u)
        res = optimizer(mec, _α, _u, use_bayes)
        @. _u = res[:u]
        return length(mec.σ²) / res[:χ²] - 1.0
    end

    println("Apply historic algorithm to determine optimized α")

    use_bayes = false
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    n_svd = length(mec.Bₘ)

    u_vec = zeros(F64, n_svd)
    s_vec = []

    conv = 0.0
    while conv < 1.0
        sol = optimizer(mec, alpha, u_vec, use_bayes)
        push!(s_vec, sol)
        alpha = alpha / ratio
        conv = length(mec.σ²) / sol[:χ²]
    end

    u_vec = s_vec[end-1][:u]
    alpha = s_vec[end][:α]
    α_opt = secant(root_fun, alpha, u_vec)

    sol = optimizer(mec, α_opt, u_vec, use_bayes)
    println("Optimized α : $α_opt log10(α) : $(log10(α_opt))")

    return s_vec, sol
end

"""
    classic(mec::MaxEntContext)

Apply the classic algorithm to solve the analytic continuation problem.

Classic algorithm uses Bayes statistics to approximately determine the
most probable value of α. We always start at a large value of α, where
the optimization yields basically the default model, therefore `u_vec`
is only a few steps away from 0 (= default model). And then we gradually
decrease α, step by step moving away from the default model towards data
fitting. Using `u_vec` as start for the next (smaller) α brings a great
speedup into this procedure.

For the classic algorithm, `alpha` is usually 10⁶, and `ratio` is 10.0.
It is incompatible with the Bayesian Reconstruction entropy.

### Arguments
* mec -> A MaxEntContext struct.

### Returns
* svec -> A vector of dictionaries. It contains the intermediate solutions.
* sol -> Dictionary. It contains the final solution.

See also: [`MaxEntContext`](@ref).
"""
function classic(mec::MaxEntContext)
    function root_fun(_α, _u)
        res = optimizer(mec, _α, _u, use_bayes)
        @. _u = res[:u]
        return res[:conv] - 1.0
    end

    println("Apply classic algorithm to determine optimized α")

    use_bayes = true
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    n_svd = length(mec.Bₘ)

    u_vec = zeros(F64, n_svd)
    s_vec = []

    conv = 0.0
    while conv < 1.0
        sol = optimizer(mec, alpha, u_vec, use_bayes)
        push!(s_vec, sol)
        alpha = alpha / ratio
        @. u_vec = sol[:u]
        conv = sol[:conv]
    end

    c_vec = [x[:conv] for x in s_vec]
    α_vec = [x[:α] for x in s_vec]
    exp_opt = log10(α_vec[end] / α_vec[end-1])
    exp_opt = exp_opt / log10(c_vec[end] / c_vec[end-1])
    exp_opt = log10(α_vec[end-1]) - log10(c_vec[end-1]) * exp_opt

    # Starting from the predicted value of α, and starting optimization
    # at the solution for the next-lowest α, we find the optimal α by
    # secant root finding method.
    u_vec = s_vec[end-1][:u]
    alpha = 10.0 ^ exp_opt
    α_opt = secant(root_fun, alpha, u_vec)

    sol = optimizer(mec, α_opt, u_vec, use_bayes)
    println("Optimized α : $α_opt log10(α) : $(log10(α_opt))")

    return s_vec, sol
end

"""
    bryan(mec::MaxEntContext)

Apply the bryan algorithm to solve the analytic continuation problem.

Bryan's maxent calculates an average of spectral functions, weighted by
their Bayesian probability.

For the bryan algorithm, `alpha` is usually 500, and `ratio` is 1.1.
It is incompatible with the Bayesian Reconstruction entropy.

### Arguments
* mec -> A MaxEntContext struct.

### Returns
* svec -> A vector of dictionaries. It contains the intermediate solutions.
* sol -> Dictionary. It contains the final solution.

See also: [`MaxEntContext`](@ref).
"""
function bryan(mec::MaxEntContext)
    println("Apply bryan algorithm to determine optimized α")

    use_bayes = true
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    n_svd = length(mec.Bₘ)
    nmesh = length(mec.mesh)

    u_vec = zeros(F64, n_svd)
    s_vec = []

    maxprob = 0.0
    while true
        sol = optimizer(mec, alpha, u_vec, use_bayes)
        push!(s_vec, sol)
        alpha = alpha / ratio
        @. u_vec = sol[:u]
        prob = sol[:prob]
        if prob > maxprob
            maxprob = prob
        elseif prob < 0.01 * maxprob
            break
        end
    end

    α_vec = map(x->x[:α], s_vec)
    p_vec = map(x->x[:prob], s_vec)
    p_vec = -p_vec ./ trapz(α_vec, p_vec)
    A_vec = map(x->x[:A], s_vec)

    nprob = length(p_vec)
    A_opt = zeros(F64, nmesh)
    spectra = zeros(F64, nmesh, nprob)
    for i = 1:nprob
        spectra[:,i] = A_vec[i] * p_vec[i]
    end
    for j = 1:nmesh
        A_opt[j] = -trapz(α_vec, spectra[j,:])
    end

    sol = Dict(:A => A_opt)

    return s_vec, sol
end

"""
    chi2kink(mec::MaxEntContext)

Apply the chi2kink algorithm to solve the analytic continuation problem.

We start with an optimization at a large value of α, where we should get
only the default model. And then, α is decreased step-by-step, until the
minimal value of α is reached. Then, we fit a function

`ϕ(x; a, b, c, d) = a + b / [1 + exp(-d*(x-c))]`,

from which the optimal α is determined by

`x_opt = c - fit_position / d`,

and

`alpha_opt = 10^x_opt`.

For the chi2kink algorithm, `alpha` is usually 10⁹, `ratio` is 10.0, the
number of alpha parameters is 12. It is compatible with the Bayesian
Reconstruction entropy.

### Arguments
* mec -> A MaxEntContext struct.

### Returns
* svec -> A vector of dictionaries. It contains the intermediate solutions.
* sol -> Dictionary. It contains the final solution.

See also: [`MaxEntContext`](@ref).
"""
function chi2kink(mec::MaxEntContext)
    function fitfun(x, p)
        return @. p[1] + p[2] / (1.0 + exp(-p[4] * (x - p[3])))
    end

    println("Apply chi2kink algorithm to determine optimized α")

    use_bayes = false
    alpha = get_m("alpha")
    ratio = get_m("ratio")
    nalph = get_m("nalph")
    α_end = alpha / (ratio^nalph)
    n_svd = length(mec.Bₘ)

    u_vec = zeros(F64, n_svd)
    s_vec = []
    χ_vec = []
    α_vec = []

    while true
        sol = optimizer(mec, alpha, u_vec, use_bayes)
        push!(s_vec, sol)
        push!(α_vec, alpha)
        push!(χ_vec, sol[:χ²])
        @. u_vec = sol[:u]
        alpha = alpha / ratio
        if alpha < α_end
            break
        end
    end

    good = isfinite.(χ_vec)
    guess = [0.0, 5.0, 2.0, 0.0]
    fit = curve_fit(fitfun, log10.(α_vec[good]), log10.(χ_vec[good]), guess)
    _, _, c, d = fit.param

    # `fit_pos` is a control parameter for under/overfitting.
    # Good values are usually between 2 and 2.5. Smaller values usually
    # lead to underfitting, which is sometimes desirable. Larger values
    # lead to overfitting, which should be avoided.
    fit_pos = 2.5
    α_opt = c - fit_pos / d
    close = argmin( abs.( log10.(α_vec) .- α_opt ) )
    u_vec = s_vec[close][:u]
    α_opt = 10.0 ^ α_opt

    sol = optimizer(mec, α_opt, u_vec, use_bayes)
    println("Optimized α : $α_opt log10(α) : $(log10(α_opt))")

    return s_vec, sol
end

"""
    optimizer(
        mec::MaxEntContext,
        α::F64,
        us::Vector{F64},
        use_bayes::Bool
    )

Optimization of maxent functional for a given value of `α`. Since a priori
the best value of `α` is unknown, this function has to be called several
times in order to find a good value.

`α` means a weight factor of the entropy. `us` is a vector in singular
space. It is used as a starting value for the optimization. For the very
first optimization, done at large α, we use zeros, which corresponds to
the default model. Then we use the result of the previous optimization
as a starting value. `use_bayes` determines whether to use the Bayesian
inference parameters for `α`.

This function will return a dictionary object that holds the results of
the optimization, e.g. spectral function, χ² deviation.

### Arguments
* mec -> A MaxEntContext struct.
* α -> See above explanations.
* us -> See above explanations.
* use_bayes -> See above explanations.

### Returns
* dict -> A dictionary, the solution to analytic continuation problem.
"""
function optimizer(
    mec::MaxEntContext,
    α::F64,
    us::Vector{F64},
    use_bayes::Bool
    )
    blur = get_m("blur")
    offdiag = get_b("offdiag")

    if offdiag
        solution, call = newton(f_and_J_od, us, mec, α)
        u = copy(solution)
        A = svd_to_real_od(mec, solution)
        S = calc_entropy_od(mec, A)
    else
        solution, call = newton(f_and_J, us, mec, α)
        u = copy(solution)
        A = svd_to_real(mec, solution)
        S = calc_entropy(mec, A, u)
    end

    χ² = calc_chi2(mec, A)
    norm = trapz(mec.mesh, A)

    dict = Dict{Symbol,Any}(
        :u => u,
        :α => α,
        :S => S,
        :χ² => χ²,
        :norm => norm,
        :Q => α * S - 0.5 * χ²,
        :Araw => deepcopy(A),
    )

    if use_bayes
        if offdiag
            ng, tr, conv, prob = calc_bayes_od(mec, A, S, χ², α)
        else
            ng, tr, conv, prob = calc_bayes(mec, A, S, χ², α)
        end
        dict[:ngood] = ng
        dict[:trace] = tr
        dict[:conv] = conv
        dict[:prob] = prob
    end

    if blur > 0.0
        make_blur(mec.mesh, A, blur)
    end
    dict[:A] = A

    @printf("log10(α) = %8.4f ", log10(α))
    @printf("χ² = %8.4e ", χ²)
    @printf("S = %8.4e ", S)
    @printf("call = %4i ", call)
    @printf("norm = %8.4f\n", norm)

    return dict
end

#=
### *Service Functions*
=#

#=
*Remarks* :

Try to calculate some key variables by using the Einstein summation trick.

```math
\begin{equation}
B_m = \sum^{N}_{n = 1} \frac{1}{\sigma^2_n} \xi_m U_{nm} G_n,
\end{equation}
```

```math
\begin{equation}
W_{ml} = \sum_{pn} \frac{1}{\sigma^2_n}
    U_{nm}\xi_m U_{np} \xi_p V_{lp} \Delta_l D_l,
\end{equation}
```

```math
\begin{equation}
W_{mli} = W_{ml} V_{li}.
\end{equation}
```

Note that these variables do not depend on the spectral function
``A(\omega)``, so they could be computed at advance to improve the
computational efficiency.

---

The `hessian matrix` is also calculated here.

```math
\begin{equation}
L = \frac{1}{2} \chi^2,
\end{equation}
```

```math
\begin{equation}
\frac{\partial^2 L}{\partial A_i \partial A_j} =
\sum_n \frac{K_{ni} \Delta_i K_{nj} \Delta_j}{\sigma^2_n}.
\end{equation}
```
=#

"""
    precompute(
        Gᵥ::Vector{F64},
        σ²::Vector{F64},
        am::AbstractMesh,
        D::Vector{F64},
        K::Matrix{F64}
    )

Precompute some key coefficients. Here `Gᵥ` and `σ²` are input data, `am`
is the mesh for spectrum, `D` is the default model, and `K` is the kernel
function.

### Arguments
* Gᵥ -> Input correlator.
* σ² -> Error bar for input correlator.
* am -> See above explanations.
* D -> See above explanations.
* K -> See above explanations.

### Returns
* V -> An orthogonal matrix from singular value decomposition of kernel.
* W₂ -> The Wₘₗ matrix.
* W₃ -> The Wₘₗᵢ tensor.
* Bₘ -> The Bₘ vector.
* hess -> The Hessian matrix.
"""
function precompute(
    Gᵥ::Vector{F64},
    σ²::Vector{F64},
    am::AbstractMesh,
    D::Vector{F64},
    K::Matrix{F64}
    )
    # Create singular value space
    U, V, S = make_singular_space(K)

    # Evaluate sizes of the arrays
    nmesh = length(am)
    n_svd = length(S)

    # Allocate memories
    W₂ = zeros(F64, n_svd, nmesh)
    W₃ = zeros(F64, n_svd, n_svd, nmesh)
    Bₘ = zeros(F64, n_svd)
    hess = zeros(F64, nmesh, nmesh)

    # Get weight of the mesh, Δωₗ.
    Δ = am.weight

    # Compute Wₘₗ
    @einsum W₂[m,l] = σ²[k] * U[k,m] * S[m] * U[k,n] * S[n] * V[l,n] * Δ[l] * D[l]

    # Compute Wₘₗᵢ
    @einsum W₃[m,k,l] = W₂[m,l] * V[l,k]

    # Compute Bₘ
    @einsum Bₘ[m] = S[m] * U[k,m] * σ²[k] * Gᵥ[k]

    # Compute the Hessian matrix
    @einsum hess[i,j] = Δ[i] * Δ[j] * K[k,i] * K[k,j] * σ²[k]
    @. hess = (hess + hess') / 2.0

    return V, W₂, W₃, Bₘ, hess
end

#=
*Remarks* :

For Shannon-Jaynes entropy,

```math
\begin{equation}
w_l = \exp \left(\sum_m V_{lm} u_m\right).
\end{equation}
```

```math
\begin{equation}
f_m = \alpha u_m + \sum_l W_{ml} w_l - B_m.
\end{equation}
```

```math
\begin{equation}
J_{mi} = \alpha \delta_{mi} + \sum_l W_{mli} w_l.
\end{equation}
```

---

For Bayesian Reconstruction entropy,

```math
\begin{equation}
w_l = \frac{1}{1 - D_l \sum_m V_{lm} u_m}.
\end{equation}
```

```math
\begin{equation}
f_m = \alpha u_m + \sum_l W_{ml} w_l - B_m.
\end{equation}
```

```math
\begin{equation}
J_{mi} = \alpha \delta_{mi} + \sum_l W_{mli} D_l w_l w_l.
\end{equation}
```
=#

"""
    f_and_J(u::Vector{F64}, mec::MaxEntContext, α::F64)

This function evaluates the function whose root we want to find. Here
`u` is a singular space vector that parametrizes the spectral function,
and `α` is a (positive) weight factor of the entropy.

It returns `f`, value of the function whose zero we want to find, and
`J`, jacobian at the current position.

### Arguments
See above explanations.

### Returns
See above explanations.

See also: [`f_and_J_od`](@ref).
"""
function f_and_J(u::Vector{F64}, mec::MaxEntContext, α::F64)
    stype = get_m("stype")

    n_svd = length(mec.Bₘ)
    J = diagm([α for i = 1:n_svd])

    # For Shannon–Jaynes entropy
    if stype == "sj"
        w = exp.(mec.Vₛ * u)
        #
        for j = 1:n_svd
            for i = 1:n_svd
                J[i,j] = J[i,j] + dot(mec.W₃[i,j,:], w)
            end
        end
        #
        f = α * u + mec.W₂ * w - mec.Bₘ
    # For Bayesian Reconstruction entropy
    else
        w = mec.Vₛ * u
        w₁ = 1.0 ./ (1.0 .- mec.model .* w)
        w₂ = w₁ .* w₁ .* mec.model
        #
        for j = 1:n_svd
            for i = 1:n_svd
                J[i,j] = J[i,j] + dot(mec.W₃[i,j,:], w₂)
            end
        end
        #
        f = α * u + mec.W₂ * w₁ - mec.Bₘ
    end

    return f, J
end

#=
*Remarks* :

For Shannon-Jaynes entropy,

```math
\begin{equation}
w_l = \exp \left(\sum_m V_{lm} u_m\right).
\end{equation}
```

```math
\begin{equation}
f_m = \alpha u_m +
      \sum_l W_{ml}\left(w_l - \frac{1}{w_l}\right) - B_m.
\end{equation}
```

```math
\begin{equation}
J_{mi} = \alpha \delta_{mi} +
         \sum_{l} W_{mli} \left(w_l + \frac{1}{w_l}\right).
\end{equation}
```

---

For Bayesian Reconstruction entropy,

```math
\begin{equation}
w^+_l = \frac{1}{ 1 - D_l \sum_m V_{lm} u_m}.
\end{equation}
```

```math
\begin{equation}
w^-_l = \frac{1}{ 1 + D_l \sum_m V_{lm} u_m}.
\end{equation}
```

```math
\begin{equation}
f_m = \alpha u_m + \sum_l W_{ml} (w^+_l - w^-_l) - B_m.
\end{equation}
```

```math
\begin{equation}
J_{mi} = \alpha \delta_{mi} + \sum_l W_{mli} D_l (w^+_l w^+_l + w^-_l w^-_l).
\end{equation}
```
=#

"""
    f_and_J_od(u::Vector{F64}, mec::MaxEntContext, α::F64)

This function evaluates the function whose root we want to find. Here
`u` is a singular space vector that parametrizes the spectral function,
and `α` is a (positive) weight factor of the entropy.

It returns `f`, value of the function whose zero we want to find, and
`J`, jacobian at the current position.

This function is similar to `f_and_J`, but for offdiagonal elements.

### Arguments
See above explanations.

### Returns
See above explanations.

See also: [`f_and_J`](@ref).
"""
function f_and_J_od(u::Vector{F64}, mec::MaxEntContext, α::F64)
    stype = get_m("stype")

    n_svd = length(mec.Bₘ)
    J = diagm([α for i = 1:n_svd])

    # For Shannon–Jaynes entropy
    if stype == "sj"
        w = exp.(mec.Vₛ * u)
        #
        a⁺ = 1.0 .* w
        a⁻ = 1.0 ./ w
        a₁ = a⁺ - a⁻
        a₂ = a⁺ + a⁻
        #
        for j = 1:n_svd
            for i = 1:n_svd
                J[i,j] = J[i,j] + dot(mec.W₃[i,j,:], a₂)
            end
        end
        #
        f = α * u + mec.W₂ * a₁ - mec.Bₘ
    # For Bayesian Reconstruction entropy
    else
        w = mec.Vₛ * u
        #
        a⁺ = 1.0 ./ (1.0 .- mec.model .* w)
        a⁻ = 1.0 ./ (1.0 .+ mec.model .* w)
        a₁ = a⁺ - a⁻
        a₂ = (a⁺ .* a⁺ + a⁻ .* a⁻) .* mec.model
        #
        for j = 1:n_svd
            for i = 1:n_svd
                J[i,j] = J[i,j] + dot(mec.W₃[i,j,:], a₂)
            end
        end
        #
        f = α * u + mec.W₂ * a₁ - mec.Bₘ
    end

    return f, J
end

#=
*Remarks* :

For Shannon-Jaynes entropy,

```math
\begin{equation}
A_l = D_l \exp \left(\sum_m V_{lm} u_m\right).
\end{equation}
```

---

For Bayesian Reconstruction entropy,

```math
\begin{equation}
A_l = \frac{D_l}{ 1 - D_l \sum_m V_{lm} u_m}.
\end{equation}
```
=#

"""
    svd_to_real(mec::MaxEntContext, u::Vector{F64})

Go from singular value space to real space. It will transform the singular
space vector `u` into real-frequency space (to get the spectral function)
by `A(ω) = D(ω) eⱽᵘ`, where `D(ω)` is the default model, `V` is the matrix
from the singular value decomposition. The argument `u` means a singular
space vector that parametrizes the spectral function.

### Arguments
See above explanations.

### Returns
See above explanations.

See also: [`svd_to_real_od`](@ref).
"""
function svd_to_real(mec::MaxEntContext, u::Vector{F64})
    stype = get_m("stype")
    #
    # For Shannon–Jaynes entropy
    if stype == "sj"
        w = exp.(mec.Vₛ * u)
        return mec.model .* w
    # For Bayesian Reconstruction entropy
    else
        w = mec.Vₛ * u
        return mec.model ./ (1.0 .- mec.model .* w)
    end
end

#=
*Remarks* :

For Shannon-Jaynes entropy,

```math
\begin{equation}
A_l = D_l \exp \left(\sum_m V_{lm} u_m\right) -
      D_l \exp \left(-\sum_m V_{lm} u_m\right).
\end{equation}
```

---

For Bayesian Reconstruction entropy,

```math
\begin{equation}
w^+_l = \frac{1}{ 1 - D_l \sum_m V_{lm} u_m}.
\end{equation}
```

```math
\begin{equation}
w^-_l = \frac{1}{ 1 + D_l \sum_m V_{lm} u_m}.
\end{equation}
```

```math
\begin{equation}
A_l = D_l (w^+_l - w^-_l).
\end{equation}
```
=#

"""
    svd_to_real_od(mec::MaxEntContext, u::Vector{F64})

Go from singular value space to real space. It will transform the singular
space vector `u` into real-frequency space in the case of an offdiagonal
element. It will return the spectral function.

### Arguments
* mec -> A MaxEntContext struct.
* u -> A singular space vector that parametrizes the spectral function.

### Returns
See above explanations.

See also: [`svd_to_real`](@ref).
"""
function svd_to_real_od(mec::MaxEntContext, u::Vector{F64})
    stype = get_m("stype")
    #
    # For Shannon–Jaynes entropy
    if stype == "sj"
        w = exp.(mec.Vₛ * u)
        w⁺ = w
        w⁻ = 1.0 ./ w
        return mec.model .* (w⁺ .- w⁻)
    # For Bayesian Reconstruction entropy
    else
        w = mec.Vₛ * u
        w⁺ = 1.0 ./ (1.0 .- mec.model .* w)
        w⁻ = 1.0 ./ (1.0 .+ mec.model .* w)
        return mec.model .* (w⁺ .- w⁻)
    end
end

#=
*Remarks* :

Shannon–Jaynes entropy

```math
\begin{equation}
S[A] = \int^{\infty}_0 d\omega
\left[
    A - m -
    A \log{\left(\frac{A}{m}\right)}
\right],
\end{equation}
```

```math
\begin{equation}
S[A^{+},A^{-}] = \int^{+\infty}_0 d\omega
\left[
    \sqrt{A^2 + 4m^2} - 2m -
    A\log{\left(\frac{\sqrt{A^2 + 4m^2} + A}{2m}\right)}
\right].
\end{equation}
```

---

Bayesian Reconstruction entropy

```math
\begin{equation}
S[A] = \int^{\infty}_0 d\omega
\left[
    1 - \frac{A}{m} + \log{\left(\frac{A}{m}\right)}
\right],
\end{equation}
```

```math
\begin{equation}
S[A^{+},A^{-}] = \int^{+\infty}_0 d\omega
\left[
    2 - \frac{\sqrt{A^2 + m^2} + m}{m} +
    \log{\left(\frac{\sqrt{A^2 + m^2} + m}{2m}\right)}
\right].
\end{equation}
```
=#

"""
    calc_entropy(mec::MaxEntContext, A::Vector{F64}, u::Vector{F64})

It computes entropy for positive definite spectral function. Here the
arguments `A` means spectral function and `u` means a singular space
vector that parametrizes the spectral function.

### Arguments
See above explanations.

### Returns
* S -> Entropy.

See also: [`calc_entropy_od`](@ref).
"""
function calc_entropy(mec::MaxEntContext, A::Vector{F64}, u::Vector{F64})
    stype = get_m("stype")
    #
    # For Shannon–Jaynes entropy
    if stype == "sj"
        f = A - mec.model - A .* (mec.Vₛ * u)
    # For Bayesian Reconstruction entropy
    else
        𝑅 = A ./ mec.model
        #
        if any(x -> x < 0.0, 𝑅)
            @info "Negative spectrum occurs!"
            @info "The results might be questionable."
            @info "Perhaps you should switch to the Shannon–Jaynes entropy."
            f = 1.0 .- 𝑅 + log.(abs.(𝑅))
        else
            f = 1.0 .- 𝑅 + log.(𝑅)
        end
    end
    #
    return trapz(mec.mesh, f)
end

"""
    calc_entropy_od(mec::MaxEntContext, A::Vector{F64})

It compute *positive-negative entropy* for spectral function with norm 0.
Here the argument `A` means spectral function.

### Arguments
See above explanations.

### Returns
* S -> Entropy.

See also: [`calc_entropy`](@ref).
"""
function calc_entropy_od(mec::MaxEntContext, A::Vector{F64})
    stype = get_m("stype")
    #
    # For Shannon–Jaynes entropy
    if stype == "sj"
        root = sqrt.(A .^ 2.0 + 4.0 .* mec.model .* mec.model)
        f = root - 2.0 .* mec.model
        f = f - A .* log.((root + A) ./ (2.0 .* mec.model))
    # For Bayesian Reconstruction entropy
    else
        root = sqrt.(A .^ 2.0 + mec.model .^ 2.0) + mec.model
        f = 2.0 .- (root ./ mec.model) + log.(root ./ (2.0 .* mec.model))
    end
    #
    return trapz(mec.mesh, f)
end

#=
*Remarks* :

**Posterior distribution of ``\alpha``**

Because
```math
\begin{equation}
\text{Pr}[\alpha | \bar{G}] =
\text{Pr}[\alpha] \frac{e^Q}{Z_L Z_S}
\frac{(2\pi)^{N/2}}{\sqrt{\det{[\alpha I + \Lambda]}}},
\end{equation}
```

so
```math
\begin{equation}
\log \text{Pr}[\alpha | \bar{G}] =
\text{constant} + \log \text{Pr} [\alpha] +
\frac{1}{2} \text{Tr} \log \left[\frac{\alpha I}{\alpha I + A}\right] +
\alpha S - \frac{1}{2}\chi^2.
\end{equation}
```

The defining equation for the `classic MaxEnt` equation reads:

```math
-2\alpha S = \text{Tr} \left(\frac{\Lambda}{\alpha I + \Lambda}\right).
```

The summation on the right hand side of the above equation is defined to
be ``N_g``, the number of good measurements:

```math
\begin{equation}
N_g = \text{Tr} \left(\frac{\Lambda}{\alpha I + \Lambda}\right).
\end{equation}
```

If ``\lambda_i`` are the eigenvalues of ``\Lambda``, then

```math
\begin{equation}
N_g = \sum_i \frac{\lambda_i}{\alpha + \lambda_i}.
\end{equation}
```

**``\Lambda`` matrix**

For Shannon-Jaynes entropy,

```math
\begin{equation}
\frac{\partial^2 S[A]}{\partial A_i \partial A_j} =
-\delta_{ij}\frac{\Delta_i}{A_i} =
-\delta_{ij}\frac{\sqrt{\Delta_i \Delta_j}}{\sqrt{A_i A_j}},
\end{equation}
```

```math
\begin{equation}
\Lambda_{ij} = \sqrt{\frac{A_i}{\Delta_i}}
               \frac{\partial^2 L}{\partial A_i \partial A_j}
               \sqrt{\frac{A_j}{\Delta_j}}.
\end{equation}
```

```math
\begin{equation}
\frac{\partial^2 S[A^+,A^-]}{\partial A_i \partial A_j} =
-\delta_{ij}\frac{\Delta_i}{\sqrt{A_i^2 + 4m_i^2}} =
-\delta_{ij}\frac{\sqrt[4]{\Delta^2_i}\sqrt[4]{\Delta^2_j}}
                 {\sqrt[4]{A_i^2 + 4m_i^2} \sqrt[4]{A_j^2 + 4m_j^2}},
\end{equation}
```

```math
\begin{equation}
\Lambda_{ij} = \sqrt[4]{\frac{A^2_i + 4m_i^2}{\Delta^2_i}}
               \frac{\partial^2 L}{\partial A_i \partial A_j}
               \sqrt[4]{\frac{A^2_j + 4m_j^2}{\Delta^2_j}}.
\end{equation}
```

---

For Bayesian Reconstruction entropy,

```math
\begin{equation}
\frac{\partial^2 S[A]}{\partial A_i \partial A_j} =
-\delta_{ij} \frac{\Delta_i}{A^2_i} =
-\delta_{ij} \frac{\sqrt{\Delta_i \Delta_j}}{A_i A_j},
\end{equation}
```

```math
\begin{equation}
\Lambda_{ij} = \frac{A_i}{\sqrt{\Delta_i}}
               \frac{\partial^2 L}{\partial A_i \partial A_j}
               \frac{A_j}{\sqrt{\Delta_j}}.
\end{equation}
```

```math
\begin{equation}
\frac{\partial^2 S[A^+,A^-]}{\partial A_i \partial A_j} =
-\delta_{ij} X_{ij} Y_{ij},
\end{equation}
```

```math
\begin{equation}
X_{ij} = \frac{\sqrt{2\Delta_i} \sqrt{2\Delta_j}}{
        \left(\sqrt{A^2_i + m^2_i} + m_i + A_i\right)
        \left(\sqrt{A^2_j + m^2_j} + m_j + A_j\right)
    },
\end{equation}
```

```math
\begin{equation}
Y_{ij} =
    \frac{
        \sqrt{A_i + \sqrt{A^2_i + m^2_i}}
        \sqrt{A_j + \sqrt{A^2_j + m^2_j}}
    }{
        \sqrt[4]{A^2_i + m^2_i}
        \sqrt[4]{A^2_j + m^2_j}
    },
\end{equation}
```

```math
\begin{equation}
\Lambda_{ij} = Z_i \frac{\partial^2 L}{\partial A_i \partial A_j} Z_j
\end{equation}
```

```math
\begin{equation}
Z_i = \frac{\left(\sqrt{A^2_i + m^2_i} + m_i + A_i\right)}{\sqrt{2\Delta_i}}
      \frac{\sqrt[4]{A^2_i + m^2_i}}{\sqrt{A_i + \sqrt{A^2_i + m^2_i}}}
\end{equation}
```

```math
\begin{equation}
Z_j = \frac{\left(\sqrt{A^2_j + m^2_j} + m_j + A_j\right)}{\sqrt{2\Delta_j}}
      \frac{\sqrt[4]{A^2_j + m^2_j}}{\sqrt{A_j + \sqrt{A^2_j + m^2_j}}}
\end{equation}
```

**Reference:**

[1] G. J. Kraberger, *et al.*, Phys. Rev. B **96**, 155128 (2017).

[2] M. Jarrell, *et al.*, Phys. Rep. **269**, 133 (1996).
=#

"""
    calc_bayes(
        mec::MaxEntContext,
        A::Vector{F64},
        S::F64,
        χ²::F64,
        α::F64
    )

It calculates Bayesian convergence criterion (`ng`, `tr`, and `conv`) for
classic maxent (maximum of probablility distribution) and then Bayesian
a-posteriori probability (`log_prob`) for `α` after optimization of `A`.

Here, `A` is the spectral function, `S` the entropy, `χ²` the deviation,
and `α` weight factor of the entropy.

### Arguments
See above explanations.

### Returns
* ng -> -2.0αS.
* tr -> Tr(Λ / (αI + Λ)).
* conv -> Ratio between `ng` and `tr`.
* prob -> Pr[α | \\bar{G}].

See also: [`calc_bayes_od`](@ref).
"""
function calc_bayes(
    mec::MaxEntContext,
    A::Vector{F64},
    S::F64,
    χ²::F64,
    α::F64
    )
    stype = get_m("stype")
    mesh = mec.mesh

    if stype == "sj"
        T = sqrt.(A ./ mesh.weight)
    else
        T = A ./ sqrt.(mesh.weight)
    end
    Λ = (T * T') .* mec.hess

    nsvd = size(mec.Vₛ, 2)
    λ = eigvals(Hermitian(Λ))[end - nsvd + 1 : end]
    filter!(x -> x > 0.0, λ)
    ng = -2.0 * α * S
    tr = sum(λ ./ (α .+ λ))
    conv = tr / ng

    eig_sum = sum(log.(α ./ (α .+ λ)))
    log_prob = α * S - 0.5 * χ² + log(α) + 0.5 * eig_sum

    return ng, tr, conv, exp(log_prob)
end

"""
    calc_bayes_od(
        mec::MaxEntContext,
        A::Vector{F64},
        S::F64,
        χ²::F64,
        α::F64
    )

It calculates Bayesian convergence criterion (`ng`, `tr`, and `conv`) for
classic maxent (maximum of probablility distribution) and then Bayesian
a-posteriori probability (`log_prob`) for `α` after optimization of `A`.

Here, `A` is the spectral function, `S` the entropy, `χ²` the deviation,
and `α` weight factor of the entropy.

It is just a offdiagonal version of `calc_bayes()`.

### Arguments
See above explanations.

### Returns
* ng -> -2.0αS.
* tr -> Tr(Λ / (αI + Λ)).
* conv -> Ratio between `ng` and `tr`.
* prob -> Pr[α | \\bar{G}].

See also: [`calc_bayes`](@ref).
"""
function calc_bayes_od(
    mec::MaxEntContext,
    A::Vector{F64},
    S::F64,
    χ²::F64,
    α::F64
    )
    stype = get_m("stype")
    mesh = mec.mesh

    if stype == "sj"
        R = (A .^ 2.0 + 4.0 * mec.model .^ 2.0) ./ (mesh.weight .^ 2.0)
        T = R .^ 0.25
    else
        R = sqrt.(A .^ 2.0 + mec.model .^ 2.0)
        X = (R .+ mec.model .+ A) ./ sqrt.(2.0 * mesh.weight)
        Y = sqrt.(R) ./ sqrt.(A .+ R)
        T = X .* Y
    end
    Λ = (T * T') .* mec.hess

    nsvd = size(mec.Vₛ, 2)
    λ = eigvals(Hermitian(Λ))[end - nsvd + 1 : end]
    filter!(x -> x > 0.0, λ)
    ng = -2.0 * α * S
    tr = sum(λ ./ (α .+ λ))
    conv = tr / ng

    eig_sum = sum(log.(α ./ (α .+ λ)))
    log_prob = α * S - 0.5 * χ² + log(α) + 0.5 * eig_sum

    return ng, tr, conv, exp(log_prob)
end

"""
    calc_chi2(mec::MaxEntContext, A::Vector{F64})

It computes the χ²-deviation of the spectral function `A`.

### Arguments
* mec -> A MaxEntContext struct.
* A -> Spectral function.

### Returns
* χ² -> Goodness-of-fit functional.
"""
function calc_chi2(mec::MaxEntContext, A::Vector{F64})
    Gₙ = reprod(mec.mesh, mec.kernel, A)
    χ² = sum(mec.σ² .* ((mec.Gᵥ - Gₙ) .^ 2.0))
    return χ²
end
