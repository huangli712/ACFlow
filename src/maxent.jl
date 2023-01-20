#
# Project : Gardenia
# Source  : maxent.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/01/21
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

Solve the analytical continuation problem by the maximum entropy method.
"""
function solve(S::MaxEntSolver, rd::RawData)
    println("[ MaxEnt ]")
    mec = init(S, rd)
    darr, sol = run(mec)
    gout = last(mec, darr, sol)
    return mec.mesh.mesh, sol[:A], gout
end

"""
    init(S::MaxEntSolver, rd::RawData)

Initialize the MaxEnt solver and return a MaxEntContext struct.
"""
function init(S::MaxEntSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value
    σ² = 1.0 ./ G.covar
    println("Postprocess input data: ", length(σ²), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    model = make_model(mesh)
    println("Build default model: ", get_b("mtype"))

    kernel = make_kernel(mesh, grid)
    println("Build default kernel: ", get_b("ktype"))

    Vₛ, W₂, W₃, Bₘ, hess = precompute(Gᵥ, σ², mesh, model, kernel)
    println("Precompute key coefficients")

    return MaxEntContext(Gᵥ, σ², grid, mesh, model, kernel, hess, Vₛ, W₂, W₃, Bₘ)
end

"""
    run(mec::MaxEntContext)

Perform maximum entropy simulation with different algorithms. Now it
supports the `historic`, `classic`, `bryan`, and `chi2kink` algorithms.
"""
function run(mec::MaxEntContext)
    method = get_m("method")

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
Here `sol` is the final solution for the analytical continuation problem,
while `svec` contains all the intermediate results (it is a vector of
dictionary actually).
"""
function last(mec::MaxEntContext, svec::Vector, sol::Dict)
    # Write the spectral function
    write_spectrum(mec.mesh, sol[:A])

    # Write the model function
    write_model(mec.mesh, mec.model)

    # Write α-χ² data
    α_vec = map(x -> x[:α], svec)
    χ_vec = map(x -> x[:χ²], svec)
    write_misfit(α_vec, χ_vec)

    # Write P[α|A] for bryan algorithm
    if haskey(svec[end], :prob)
        p_vec = map(x -> x[:prob], svec)
        write_probability(α_vec, p_vec)
    end

    # Regenerate the input data and write them
    Aout = haskey(sol, :Araw) ? sol[:Araw] : sol[:A]
    G = reprod(mec.mesh, mec.kernel, Aout)
    write_backward(mec.grid, G)

    # Calculate full response function on real axis and write them
    if get_b("ktype") == "fermi"
        _G = kramers(mec.mesh, Aout)
    else
        _G = kramers(mec.mesh, Aout .* mec.mesh)
    end
    write_complete(mec.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

"""
    historic(mec::MaxEntContext)

Apply the historic algorithm to solve the analytical continuation problem.
It choose α in a way that χ² ≈ N.

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

Apply the classic algorithm to solve the analytical continuation problem.

Classic algorithm uses Bayes statistics to approximately determine the
most probable value of α. We always start at a large value of α, where
the optimization yields basically the default model, therefore `u_vec`
is only a few steps away from 0 (= default model). And then we gradually
decrease α, step by step moving away from the default model towards data
fitting. Using `u_vec` as start for the next (smaller) α brings a great
speedup into this procedure.

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

Apply the bryan algorithm to solve the analytical continuation problem.

Bryan's maxent calculates an average of spectral functions, weighted by
their Bayesian probability.

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

Apply the chi2kink algorithm to solve the analytical continuation problem.

We start with an optimization at a large value of α, where we should get
only the default model. And then, α is decreased step-by-step, until the
minimal value of α is reached. Then, we fit a function

`ϕ(x; a, b, c, d) = a + b / [1 + exp(-d*(x-c))]`,

from which the optimal α is determined by

`x_opt = c - fit_position / d`,

and

`alpha_opt = 10^x_opt`.

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
    optimizer(mec::MaxEntContext, α::F64, us::Vector{F64}, use_bayes::Bool)

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
"""
function optimizer(mec::MaxEntContext, α::F64, us::Vector{F64}, use_bayes::Bool)
    blur = get_m("blur")
    offdiag = get_b("offdiag")

    if offdiag
        solution, call = newton(f_and_J_offdiag, us, mec, α)
        u = copy(solution)
        A = svd_to_real_offdiag(mec, solution)
        S = calc_entropy_offdiag(mec, A)
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
            ng, tr, conv, prob = calc_bayes_offdiag(mec, A, S, χ², α)
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

"""
    precompute(Gᵥ::Vector{F64}, σ²::Vector{F64},
               am::AbstractMesh,
               D::Vector{F64},
               K::Matrix{F64})

Precompute some key coefficients. Here `Gᵥ` and `σ²` are input data, `am`
is the mesh for spectrum, `D` is the default model, and `K` is the kernel
function.
"""
function precompute(Gᵥ::Vector{F64}, σ²::Vector{F64},
                    am::AbstractMesh,
                    D::Vector{F64},
                    K::Matrix{F64})
    offdiag = get_b("offdiag")

    U, V, S = make_singular_space(K)

    nmesh = length(am)
    weight = am.weight
    n_svd = length(S)

    W₂ = zeros(F64, n_svd, nmesh)
    W₃ = zeros(F64, n_svd, n_svd, nmesh)
    Bₘ = zeros(F64, n_svd)
    hess = zeros(F64, nmesh, nmesh)

    if offdiag
        @einsum W₂[m,l] = σ²[k] * U[k,m] * S[m] * U[k,n] * S[n] * V[l,n] * weight[l]
    else
        @einsum W₂[m,l] = σ²[k] * U[k,m] * S[m] * U[k,n] * S[n] * V[l,n] * weight[l] * D[l]
    end

    @einsum W₃[j,k,i] = W₂[j,i] * V[i,k]

    @einsum Bₘ[m] = S[m] * U[k,m] * σ²[k] * Gᵥ[k]

    @einsum hess[i,j] = weight[i] * weight[j] * K[k,i] * K[k,j] * σ²[k]

    return V, W₂, W₃, Bₘ, hess
end

"""
    f_and_J(u::Vector{F64}, mec::MaxEntContext, α::F64)

This function evaluates the function whose root we want to find. Here
`u` is a singular-space vector that parametrizes the spectral function,
and `α` is a (positive) weight factor of the entropy.

It returns `f`, value of the function whose zero we want to find, and
`J`, jacobian at the current position.

See also: [`f_and_J_offdiag`](@ref).
"""
function f_and_J(u::Vector{F64}, mec::MaxEntContext, α::F64)
    w = exp.(mec.Vₛ * u)

    n_svd = length(mec.Bₘ)

    J = diagm([α for i = 1:n_svd])
    for j = 1:n_svd
        for i = 1:n_svd
            J[i,j] = J[i,j] + dot(mec.W₃[i,j,:], w)
        end
    end

    f = α * u + mec.W₂ * w - mec.Bₘ

    return f, J
end

"""
    f_and_J_offdiag(u::Vector{F64}, mec::MaxEntContext, α::F64)

This function evaluates the function whose root we want to find. Here
`u` is a singular-space vector that parametrizes the spectral function,
and `α` is a (positive) weight factor of the entropy.

It returns `f`, value of the function whose zero we want to find, and
`J`, jacobian at the current position.

This function is similar to `f_and_J`, but for offdiagonal elements.

See also: [`f_and_J`](@ref).
"""
function f_and_J_offdiag(u::Vector{F64}, mec::MaxEntContext, α::F64)
    w = exp.(mec.Vₛ * u)

    a_plus = mec.model .* w
    a_minus = mec.model ./ w
    a1 = a_plus - a_minus
    a2 = a_plus + a_minus

    n_svd = length(mec.Bₘ)

    J = diagm([α for i = 1:n_svd])
    for j = 1:n_svd
        for i = 1:n_svd
            J[i,j] = J[i,j] + dot(mec.W₃[i,j,:], a2)
        end
    end

    f = α * u + mec.W₂ * a1 - mec.Bₘ

    return f, J
end

"""
    svd_to_real(mec::MaxEntContext, u::Vector{F64})

Go from singular value space to real space. It will transform the singular
space vector `u` into real-frequency space (to get the spectral function)
by `A(ω) = D(ω) eⱽᵘ`, where `D(ω)` is the default model `V` is the matrix
from the singular value decomposition. The argument `u` means a singular
space vector that parametrizes the spectral function.

See also: [`svd_to_real_offdiag`](@ref).
"""
function svd_to_real(mec::MaxEntContext, u::Vector{F64})
    return mec.model .* exp.(mec.Vₛ * u)
end

"""
    svd_to_real_offdiag(mec::MaxEntContext, u::Vector{F64})

Go from singular value space to real space. It will transform the singular
space vector `u` into real-frequency space in the case of an offdiagonal
element. It will return the spectral function.

See also: [`svd_to_real`](@ref).
"""
function svd_to_real_offdiag(mec::MaxEntContext, u::Vector{F64})
    w = exp.(mec.Vₛ * u)
    return (mec.model .* w) - (mec.model ./ w)
end

"""
    calc_entropy(mec::MaxEntContext, A::Vector{F64}, u::Vector{F64})

It computes entropy for positive definite spectral function. Here the
arguments `A` means spectral function and `u` means a singular space
vector that parametrizes the spectral function.

See also: [`calc_entropy_offdiag`](@ref).
"""
function calc_entropy(mec::MaxEntContext, A::Vector{F64}, u::Vector{F64})
    f = A - mec.model - A .* (mec.Vₛ * u)
    return trapz(mec.mesh, f)
end

"""
    calc_entropy_offdiag(mec::MaxEntContext, A::Vector{F64})

It compute *positive-negative entropy* for spectral function with norm 0.
Here the argument `A` means spectral function.

See also: [`calc_entropy`](@ref).
"""
function calc_entropy_offdiag(mec::MaxEntContext, A::Vector{F64})
    root = sqrt.(A .^ 2.0 + 4.0 .* mec.model .* mec.model)
    f = root - 2.0 .* mec.model - A .* log.((root + A) ./ (2.0 .* mec.model))
    return trapz(mec.mesh, f)
end

"""
    calc_bayes(mec::MaxEntContext,
               A::Vector{F64},
               S::F64, χ²::F64, α::F64)

It calculates Bayesian convergence criterion (`ng`, `tr`, and `conv`) for
classic maxent (maximum of probablility distribution) and then Bayesian
a-posteriori probability (`log_prob`) for `α` after optimization of `A`.

Here, `A` is the spectral function, `S` the entropy, `χ²` the deviation,
and `α` weight factor of the entropy.

See also: [`calc_bayes_offdiag`](@ref).
"""
function calc_bayes(mec::MaxEntContext,
                    A::Vector{F64},
                    S::F64, χ²::F64, α::F64)
    mesh = mec.mesh

    T = sqrt.(A ./ mesh.weight)
    Λ = (T * T') .* mec.hess

    λ = eigvals(Hermitian(Λ))
    ng = -2.0 * α * S
    tr = sum(λ ./ (α .+ λ))
    conv = tr / ng

    eig_sum = sum(log.(α ./ (α .+ λ)))
    log_prob = α * S - 0.5 * χ² + log(α) + 0.5 * eig_sum

    return ng, tr, conv, exp(log_prob)
end

"""
    calc_bayes_offdiag(mec::MaxEntContext,
                       A::Vector{F64},
                       S::F64, χ²::F64, α::F64)

It calculates Bayesian convergence criterion (`ng`, `tr`, and `conv`) for
classic maxent (maximum of probablility distribution) and then Bayesian
a-posteriori probability (`log_prob`) for `α` after optimization of `A`.

Here, `A` is the spectral function, `S` the entropy, `χ²` the deviation,
and `α` weight factor of the entropy.

It is just a offdiagonal version of `calc_bayes()`.

See also: [`calc_bayes`](@ref).
"""
function calc_bayes_offdiag(mec::MaxEntContext,
                            A::Vector{F64},
                            S::F64, χ²::F64, α::F64)
    mesh = mec.mesh

    T = (( A .^ 2.0 + 4.0 * mec.model .* mec.model ) / (mesh.weight .^ 2.0)) .^ 0.25
    Λ = (T * T') .* mec.hess

    λ = eigvals(Hermitian(Λ))
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
"""
function calc_chi2(mec::MaxEntContext, A::Vector{F64})
    Gₙ = reprod(mec.mesh, mec.kernel, A)
    χ² = sum(mec.σ² .* ((mec.Gᵥ - Gₙ) .^ 2.0))
    return χ²
end
