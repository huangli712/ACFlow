#
# Project : Gardenia
# Source  : rfa.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/08/01
#

#
# Note:
#
# The following codes for the Barycentric rational function approximation
# are mostly adapted from
#
#     https://github.com/complexvariables/RationalFunctionApproximation.jl
#
# See
#
#     The AAA algorithm for rational approximation
#     Yuji Nakatsukasa, Olivier Sète, Lloyd N. Trefethen
#     SIAM Journal on Scientific Computing 40, A1494 (2018)
#
# and
#
#     AAA Rational Approximation on a Continuum
#     Tobin A. Driscoll, Yuji Nakatsukasa, Lloyd N. Trefethen
#     SIAM Journal on Scientific Computing 46, A929 (2024)
#
# for more details.
#

#=
### *Customized Structs* : *BarycentricFunction*
=#

#=
*Remarks* :

**Rational Barycentric Representation**

The barycentric formula takes the form of a quotient of two partial
fractions,

```math
r(z) = \frac{n(z)}{d(z)}
     = \sum^m_{j=1} \frac{w_j f_j}{z - z_j}
     {\huge/} \sum^m_{j=1} \frac{w_j}{z - z_j},
```

where ``m \ge 1`` is an integer, ``z_1, \cdots, z_m`` are a set of real
or complex distinct support points (`nodes`), ``f_1, \cdots, f_m`` are a
set of real or complex data `values`, and ``w_1, \cdots, w_m`` are a set
of real or complex `weights`. As indicated in this equation, we just let
``n(z)`` and ``d(z)`` stand for the partial fractions in the numerator
and the denominator.
=#

"""
    BarycentricFunction

Mutable struct. Barycentric representation of a rational function.

### Members

* nodes     -> Nodes of the rational function, ``z_i``.
* values    -> Values of the rational function, ``r(z_i)``.
* weights   -> Weights of the rational function, ``w_i``.
* w_times_f -> Weighted values of the rational function, ``w_i f_i``.
"""
mutable struct BarycentricFunction <: Function
    nodes     :: Vector{C64}
    values    :: Vector{C64}
    weights   :: Vector{C64}
    w_times_f :: Vector{C64}
end

"""
    BarycentricFunction(
        nodes   :: Vector{C64},
        values  :: Vector{C64},
        weights :: Vector{C64}
        )

Construct a `BarycentricFunction` type rational function.

### Arguments
* nodes   -> Interpolation nodes, ``z_i``.
* values  -> Values at the interpolation nodes, ``r(z_i)``.
* weights -> Barycentric weights, ``w_i``.
"""
function BarycentricFunction(
    nodes   :: Vector{C64},
    values  :: Vector{C64},
    weights :: Vector{C64}
    )
    @assert length(nodes) == length(values) == length(weights)
    w_times_f = values .* weights
    return BarycentricFunction(nodes, values, weights, w_times_f)
end

"""
    bc_nodes(r::BarycentricFunction)

Returns the nodes of the rational interpolant `r` as a vector. Actually,
they are ``z_i``.
"""
bc_nodes(r::BarycentricFunction) = r.nodes

"""
    bc_values(r::BarycentricFunction)

Returns the nodal values of the rational interpolant `r` as a vector. They
are ``r(z_i) \equiv f_i``.
"""
bc_values(r::BarycentricFunction) = r.values

"""
    bc_weights(r::BarycentricFunction)

Returns the weights of the rational interpolant `r` as a vector. Actually,
they are ``w_i``.
"""
bc_weights(r::BarycentricFunction) = r.weights

"""
    bc_degree(r::BarycentricFunction)

Returns the degree of the numerator and denominator of the rational `r`.
"""
bc_degree(r::BarycentricFunction) = length(r.nodes) - 1

"""
    bc_poles(r::BarycentricFunction)

Return the poles of the rational function `r`.
"""
function bc_poles(r::BarycentricFunction)
    w = bc_weights(r)
    z = bc_nodes(r)
    nonzero = @. !iszero(w)
    z, w = z[nonzero], w[nonzero]
    #
    m = length(w)
    B = diagm( [zero(F64); ones(F64, m)] )
    E = [zero(F64) transpose(w); ones(F64, m) diagm(z) ];
    #
    pole = [] # Put it into scope
    try
        pole = filter( isfinite, eigvals(E, B) )
    catch
        # Generalized eigen not available in extended precision, so:
        λ = filter( z->abs(z)>1e-13, eigvals(E\B) )
        pole = 1 ./ λ
    end

    return pole
end

"""
    (r::BarycentricFunction)(z::Number)

Evaluate the Barycentric rational function at `z`.
"""
function (r::BarycentricFunction)(z::Number)
    if isinf(z)
        return sum(r.w_times_f) / sum(r.weights)
    end
    #
    # Try to determine whether z is a valid node
    k = findfirst(z .== r.nodes)
    #
    if isnothing(k) # Not at a node
        C = @. 1 / (z - r.nodes)
        return sum(C .* r.w_times_f) / sum(C .* r.weights)
    else            # Interpolation at node
        return r.values[k]
    end
end

#=
### *Adaptive Antoulas-Anderson Algorithm*
=#

"""
    aaa(z, y)

Adaptively compute a Barycentric rational interpolant.

### Arguments
* z::AbstractVector{<:Number} -> Interpolation nodes.
* y::AbstractVector{<:Number} -> Values at nodes.
* max_degree::Integer=150 -> Maximum numerator/denominator degree to use.
* float_type::Type=F64 -> Floating point type to use for the computation.
* tol::Real=1000*eps(float_type) -> Tolerance for stopping.
* lookahead::Integer=10 -> Number of iterations to determines stagnation.
* stats::Bool=false -> Return convergence statistics.

### Returns
* `r::BarycentricFunction` -> The rational interpolant.
* `stats::NamedTuple` -> Convergence statistics, if keyword `stats = true`.

### Examples
```julia-repl
julia> z = 1im * range(-10, 10, 500);

julia> y = @. exp(z);

julia> r = aaa(z, y);

julia> bc_degree(r) # both numerator and denominator
12

julia> first(nodes(r), 4)
4-element Vector{ComplexF64}:
 0.0 - 6.272545090180361im
 0.0 + 9.43887775551102im
 0.0 - 1.1022044088176353im
 0.0 + 4.909819639278557im

julia> r(1im * π / 2)
-2.637151617496356e-15 + 1.0000000000000002im
```
"""
function aaa(
    z::AbstractVector{<:Number},
    y::AbstractVector{<:Number};
    max_degree = 150,
    float_type = F64,
    tol = 1000 * eps(float_type),
    lookahead = 10,
    stats = false
    )

    @assert float_type <: AbstractFloat
    T = float_type
    fmax = norm(y, Inf) # For scaling
    m = length(z)
    iteration = NamedTuple[]
    err = T[]
    besterr, bestidx, best = Inf, NaN, nothing

    # Allocate space for Cauchy matrix, Loewner matrix, and residual
    C = similar(z, (m, m))
    L = similar(z, (m, m))
    R = complex(zeros(size(z)))

    ȳ = sum(y) / m
    s, idx = findmax(abs(y - ȳ) for y in y)
    push!(err, s)

    # The ordering of nodes matters,
    # while the order of test points does not.
    node_index = Int[]
    push!(node_index, idx)
    test_index = Set(1:m)
    delete!(test_index, idx)

    n = 0 # Number of poles
    while true
        n += 1
        σ = view(z, node_index)
        fσ = view(y, node_index)
        # Fill in matrices for the latest node
        @inbounds @fastmath for i in test_index
            δ = z[i] - σ[n]
            # δ can be zero if there are repeats in z
            C[i, n] = iszero(δ) ? 1 / eps() : 1 / δ
            L[i, n] = (y[i] - fσ[n]) * C[i, n]
        end
        istest = collect(test_index)
        _, _, V = svd( view(L, istest, 1:n) )
        w = V[:, end] # Barycentric weights

        CC = view(C, istest, 1:n)
        num = CC * (w.*fσ)
        den = CC * w
        @. R[istest] = y[istest] - num / den
        push!(err, norm(R, Inf))
        push!(iteration, (; weights=w, active=copy(node_index)))

        if (Base.last(err) < besterr)
            besterr = Base.last(err)
            bestidx = length(iteration)
            best = Base.last(iteration)
        end

        # Are we done?
        if (besterr <= tol*fmax) ||
            (n == max_degree + 1) ||
            ((length(iteration) - bestidx >= lookahead) && (besterr < 1e-2*fmax))
            break
        end

        _, j = findmax(abs, R)
        push!(node_index, j)
        delete!(test_index, j)
        R[j] = 0
    end

    idx, w = best.active, best.weights
    r = BarycentricFunction(z[idx], y[idx], w)
    if stats
        return r, (;err, iteration)
    else
        return r
    end
end

#
# Note:
#
# The following codes for the Prony approximation are mostly adapted from
#
#     https://github.com/Green-Phys/PronyAC
#
# See
#
#     Minimal Pole Representation and Controlled Analytic Continuation
#     of Matsubara Response Functions
#     Lei Zhang and Emanuel Gull
#     arXiv:2312.10576 (2024)
#
# for more details.
#

#=
### *Customized Structs* : *PronyApproximation*
=#

#=
*Remarks* :

**Prony interpolation**

Our input data consists of an odd number ``2N + 1`` of Matsubara points
``G(i\omega_n)`` that are uniformly spaced. Prony's interpolation method
interpolates ``G_k`` as a sum of exponentials

```math
G_k = \sum^{N-1}_{i=0} w_i \gamma^k_i,
```

where ``0 \le k \le 2N``, ``w_i`` denote complex weights and ``\gamma_i``
corresponding nodes.

---

**Prony approximation**

Prony's interpolation method is unstable. We therefore employs a Prony
approximation, rather than an interpolation of ``G``. For the physical
Matsubara functions, which decay in magnitude to zero for
``i\omega_n \to i\infty``, only ``K \propto \log{1/\varepsilon}`` out of
all ``N`` nodes in the Prony approximation have weights
``|w_i| > \varepsilon``. Thus, we have

```math
\left|G_k - \sum^{K-1}_{i=0} w_i \gamma^k_i\right| \le \varepsilon,
```

for all ``0 \le k \le 2N``.
=#

"""
    PronyApproximation

Mutable struct. Prony approximation to a complex-valued Matsubara function.

### Members

* 𝑁ₚ -> Number of nodes for Prony approximation.
* ωₚ -> Non-negative Matsubara frequency.
* 𝐺ₚ -> Complex values at ωₚ.
* Γₚ -> Nodes for Prony approximation, ``γ_i``.
* Ωₚ -> Weights for Prony approximation, ``w_i``.
"""
mutable struct PronyApproximation <: Function
    𝑁ₚ :: I64
    ωₚ :: Vector{F64}
    𝐺ₚ :: Vector{C64}
    Γₚ :: Vector{C64}
    Ωₚ :: Vector{C64}
end

"""
    PronyApproximation(
        𝑁ₚ :: I64,
        ωₚ :: Vector{F64},
        𝐺ₚ :: Vector{C64},
        v  :: Vector{C64}
        )

Construct a `PronyApproximation` type interpolant function. Once it is
available, then it can be used to produce a smooth G at given ω.

This function should not be called directly by the users.

### Arguments
* 𝑁ₚ -> Number of nodes for Prony approximation.
* ωₚ -> Non-negative Matsubara frequency (postprocessed).
* 𝐺ₚ -> Complex values at ωₚ (postprocessed).
* v  -> Selected vector from the orthogonal matrix `V`.
"""
function PronyApproximation(
    𝑁ₚ :: I64,
    ωₚ :: Vector{F64},
    𝐺ₚ :: Vector{C64},
    v  :: Vector{C64}
    )
    # Evaluate cutoff for Γₚ
    Λ = 1.0 + 0.5 / 𝑁ₚ

    # Evaluate Γₚ and Ωₚ
    Γₚ = prony_gamma(v, Λ)
    Ωₚ = prony_omega(𝐺ₚ, Γₚ)

    # Sort Γₚ and Ωₚ
    idx_sort = sortperm(abs.(Ωₚ))
    reverse!(idx_sort)
    Ωₚ = Ωₚ[idx_sort]
    Γₚ = Γₚ[idx_sort]

    # Return a PronyApproximation object
    return PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, Γₚ, Ωₚ)
end

"""
    PronyApproximation(ω₁::Vector{F64}, 𝐺₁::Vector{C64}, ε::F64)

Construct a `PronyApproximation` type interpolant function. Once it is
available, then it can be used to produce a smooth G at ω.

If the noise level of the input data is known, this function is a good
choice. The parameter `ε` can be set to the noise level.

### Arguments
* ω₁ -> Non-negative Matsubara frequency (raw).
* 𝐺₁ -> Complex values at ωₚ (raw).
* ε  -> Threshold for the Prony approximation.
"""
function PronyApproximation(ω₁::Vector{F64}, 𝐺₁::Vector{C64}, ε::F64)
    # Preprocess the input data to get the number of nodes, frequency
    # points ωₚ, and Matsubara data 𝐺ₚ.
    𝑁ₚ, ωₚ, 𝐺ₚ = prony_data(ω₁, 𝐺₁)

    # Perform singular value decomposition and select reasonable `v`.
    S, V = prony_svd(𝑁ₚ, 𝐺ₚ)
    v = prony_v(S, V, ε)

    return PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, v)
end

"""
    PronyApproximation(ω₁::Vector{F64}, 𝐺₁::Vector{C64})

Construct a `PronyApproximation` type interpolant function. Once it is
available, then it can be used to produce a smooth G at ω. Note that this
function employs a smart and iterative algorithm to determine the optimal
Prony approximation.

This function is time-consuming. But if the noise level of the input data
is unknown, this function is useful.

### Arguments
* ω₁ -> Non-negative Matsubara frequency (raw).
* 𝐺₁ -> Complex values at ωₚ (raw).
"""
function PronyApproximation(ω₁::Vector{F64}, 𝐺₁::Vector{C64})
    # Preprocess the input data to get the number of nodes, frequency
    # points ωₚ, and Matsubara data 𝐺ₚ.
    𝑁ₚ, ωₚ, 𝐺ₚ = prony_data(ω₁, 𝐺₁)

    # Perform singular value decomposition
    S, V = prony_svd(𝑁ₚ, 𝐺ₚ)

    # Next we should determine the optimal `v`
    #
    # (1) Find maximum index for the exponentially decaying region.
    exp_idx = prony_idx(S)
    #
    # (2) Find minimum index
    ε = 1000 * S[exp_idx]
    new_idx = findfirst(x -> x < ε, S)
    #
    # (3) Create lists for chosen indices and the corresponding errors.
    idxrange = range( new_idx, min(exp_idx + 10, length(S)) )
    idx_list = collect(idxrange)
    err_list = zeros(F64, length(idx_list))
    #
    # (4) Create a list of pseudo-PronyApproximations, and then evaluate
    # their reliabilities and accuracies.
    for i in eachindex(idx_list)
        idx = idx_list[i]
        #
        # Extract `v`
        v = prony_v(V, idx)
        #
        # Reproduce G using pseudo PronyApproximation
        𝐺ₙ = PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, v)(ωₚ)
        #
        # Evaluate the difference and record it
        err_ave = mean(abs.(𝐺ₙ - 𝐺ₚ))
        err_list[i] = err_ave
        #
        @show i, idx, err_ave
    end
    #
    # (5) Find the optimal `v`, which should minimize |𝐺ₙ - 𝐺ₚ|
    idx = idx_list[argmin(err_list)]
    v = prony(V, idx)

    return PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, v)
end

"""
    prony_data(ω₁, 𝐺₁)

Prepare data for later Prony approximation. It will return the number
of nodes, frequency mesh ωₚ, and Green's function data 𝐺ₚ at this mesh.
"""
function prony_data(ω₁, 𝐺₁)
    # We have to make sure the number of data points is odd.
    osize = length(ω₁)
    nsize = iseven(osize) ? osize - 1 : osize
    #
    𝑁ₚ = div(nsize, 2) # Number of nodes for Prony approximation
    ωₚ = ω₁[1:nsize]   # Matsubara frequency, ωₙ
    𝐺ₚ = 𝐺₁[1:nsize]   # Matsubara Green's function, G(iωₙ)
    #
    return 𝑁ₚ, ωₚ, 𝐺ₚ
end

"""
    prony_svd(𝑁ₚ, 𝐺ₚ)

Perform singular value decomposition for the matrix ℋ that is constructed
from 𝐺ₚ. It will return the singular values `S` and orthogonal matrix `V`.
"""
function prony_svd(𝑁ₚ, 𝐺ₚ)
    ℋ = zeros(C64, 𝑁ₚ + 1, 𝑁ₚ + 1)
    #
    for i = 1 : 𝑁ₚ + 1
        ℋ[i,:] = 𝐺ₚ[i:i+𝑁ₚ]
    end
    #
    _, S, V = svd(ℋ)

    return S, V
end

function prony_idx(S::Vector{F64}, ε::F64)
    # Write singular values
    println("List of singular values:")
    for i in eachindex(S)
        @printf("%4i %16.12e\n", i, S[i])
    end

    # Determine idx, such that S[idx] < ε.
    idx = findfirst(x -> x < ε, S)

    # Check idx
    if isnothing(idx)
        error("Please increase ε and try again!")
    end

    return idx
end

function prony_idx(S::Vector{F64})
    n_max = min(3 * floor(I64, log(1.0e12)), floor(I64, 0.8 * length(S)))
    #@show n_max
    idx_fit = collect(range(ceil(I64, 0.8*n_max), n_max))
    val_fit = S[idx_fit]
    A = hcat(idx_fit, ones(I64, length(idx_fit)))
    a, b = pinv(A) * log.(val_fit)
    S_approx = exp.(a .* collect(range(1,n_max)) .+ b)
    idx = count(S[1:n_max] .> 5.0 * S_approx) + 1
    return idx
end

"""
    prony_v(V, idx::I64)

Extract suitable vector `v` from orthogonal matrix `V` according to the
threshold `ε`. The diagonal matrix (singular values) `S` is used to test
whether the threshold `ε` is reasonable and figure out the index for
extracting `v` from `V`.
"""
function prony_v(V, idx::I64)
    # Extract v from V
    println("Selected vector from orthogonal matrix V: ", idx)
    v = V[:,idx]

    return reverse!(v)
end

"""
    prony_gamma(v, Λ)

Try to calculate Γₚ. Actually, Γₚ are eigenvalues of a matrix constructed
by `v`. `Λ` is a cutoff for Γₚ. Only those Γₚ that are smaller than `Λ`
are kept.
"""
function prony_gamma(v, Λ)
    # The following codes actually calculate the roots of a polynominal
    # with coefficients given in v. The roots are Γₚ.
    non_zero = findall(!iszero, v)
    trailing_zeros = length(v) - non_zero[end]
    #
    vnew = v[non_zero[1]:non_zero[end]]
    N = length(vnew)
    #
    if N > 1
        A = diagm(-1=>ones(C64, N - 2))
        @. A[1,:] = -vnew[2:end] / vnew[1]
        roots = eigvals(A)
    else
        roots = []
    end
    #
    Γₚ = vcat(roots, zeros(C64, trailing_zeros))

    filter!(x -> abs(x) < Λ, Γₚ)

    return Γₚ
end

"""
    prony_omega(𝐺ₚ, Γₚ)

Try to calculate Ωₚ.
"""
function prony_omega(𝐺ₚ, Γₚ)
    A = zeros(C64, length(𝐺ₚ), length(Γₚ))
    #
    for i in eachindex(𝐺ₚ)
        A[i,:] = Γₚ .^ (i - 1)
    end
    #
    return pinv(A) * 𝐺ₚ
end

"""
    (p::PronyApproximation)(w::Vector{F64})

Evaluate the Prony approximation at `w`.
"""
function (p::PronyApproximation)(w::Vector{F64})
    x₀ = @. (w - w[1]) / (w[end] - w[1])
    A = zeros(C64, length(x₀), length(p.Ωₚ))
    #
    for i in eachindex(x₀)
        @. A[i,:] = p.Γₚ ^ (2.0 * p.𝑁ₚ * x₀[i])
    end
    #
    return A * p.Ωₚ
end

#=
### *Customized Structs* : *BarRat Solver*
=#

"""
    BarRatContext

Mutable struct. It is used within the BarRat solver only.

### Members

* Gᵥ   -> Input data for correlator.
* grid -> Grid for input data.
* mesh -> Mesh for output spectrum.
* 𝒫    -> Prony approximation for the input data.
* ℬ    -> Barycentric rational function approximation for the input data.
"""
mutable struct BarRatContext
    Gᵥ   :: Vector{C64}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    𝒫    :: Union{Missing,PronyApproximation}
    ℬ    :: Union{Missing,BarycentricFunction}
end

#=
### *Global Drivers*
=#

"""
    solve(S::BarRatSolver, rd::RawData)

Solve the analytic continuation problem by the Barycentric rational
function method.
"""
function solve(S::BarRatSolver, rd::RawData)
    println("[ BarRat ]")
    #
    brc = init(S, rd)
    run(brc)
    Aout, Gout = last(brc)
    #
    return brc.mesh.mesh, Aout, Gout
end

"""
    init(S::BarRatSolver, rd::RawData)

Initialize the BarRat solver and return a BarRatContext struct.
"""
function init(S::BarRatSolver, rd::RawData)
    # Prepera input data
    Gᵥ = rd.value
    println("Postprocess input data: ", length(Gᵥ), " points")

    # Prepare grid for input data
    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    # Prepare mesh for output spectrum
    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    return BarRatContext(Gᵥ, grid, mesh, missing, missing)
end

"""
    run(brc::BarRatContext)

At first, it will try to construct a Prony approximation for the input
Matsubara data. Then the Prony approximation is used to build smooth data
set (data denoising). Finally, the Barycentric rational function for this
data set is constructed. The member `ℬ` of the BarRatContext object
(`brc`) should be updated in this function.
"""
function run(brc::BarRatContext)
    denoise = get_r("denoise")
    ε = get_r("epsilon")

    ω = brc.grid.ω
    iω = ω * im
    G = brc.Gᵥ

    if denoise == "prony"
        println("Activate Prony approximation to denoise the input data")
        pa = PronyApproximation(ω, G, ε)
        #pa = PronyApproximation(ω, G)
        brc.𝒫 = pa
        #
        println("Construct Barycentric rational function approximation")
        brc.ℬ = aaa(iω, pa(ω))
    else
        brc.ℬ = aaa(iω, G)
    end
end

"""
    last(brc::BarRatContext)

It will process and write the calculated results by the BarRat solver,
including correlator at real axis, final spectral function, reproduced
correlator. The information about Prony approximation and Barycentric
rational function approximation will be written as well.
"""
function last(brc::BarRatContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Write information about Prony approximation
    fwrite && (get_r("denoise") == "prony") && begin
        write_prony(brc.𝒫.𝑁ₚ, brc.𝒫.Γₚ, brc.𝒫.Ωₚ)
    end

    # Write information about Barycentric rational function
    fwrite && write_barycentric(brc.ℬ.nodes, brc.ℬ.values, brc.ℬ.weights)

    # Calculate full response function on real axis and write them
    _G = brc.ℬ.(brc.mesh.mesh)
    fwrite && write_complete(brc.mesh, _G)

    # Calculate and write the spectral function
    Aout = -imag.(_G) ./ π
    fwrite && write_spectrum(brc.mesh, Aout)

    # Regenerate the input data and write them
    kernel = make_kernel(brc.mesh, brc.grid)
    G = reprod(brc.mesh, kernel, Aout)
    fwrite && write_backward(brc.grid, G)

    return Aout, _G
end
