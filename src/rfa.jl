#
# Project : Gardenia
# Source  : rfa.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/07/29
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
set of real or complex data values, and ``w_1, \cdots, w_m`` are a set
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
    BarycentricFunction(nodes, values, weights)

Construct a `BarycentricFunction` type rational function.

### Arguments
* `nodes::Vector`   -> Interpolation nodes, ``z_i``.
* `values::Vector`  -> Values at the interpolation nodes, ``r(z_i)``.
* `weights::Vector` -> Barycentric weights, ``w_i``.
"""
function Barycentric(
    nodes   :: Vector{C64},
    values  :: Vector{C64},
    weights :: Vector{C64}
    )
    @assert length(nodes) == length(values) == length(weights)
    w_times_f = values .* weights
    return BarycentricFunction(nodes, values, weights, w_times_f)
end

"""
    bc_nodes(r::Barycentric)

Returns the nodes of the rational interpolant `r` as a vector.
"""
bc_nodes(r::Barycentric) = r.nodes

"""
    bc_values(r::Barycentric)

Returns the nodal values of the rational interpolant `r` as a vector.
"""
bc_values(r::Barycentric) = r.values

"""
    bc_weights(r::Barycentric)

Returns the weights of the rational interpolant `r` as a vector.
"""
bc_weights(r::Barycentric) = r.weights

"""
    bc_degree(r::Barycentric)

Returns the degree of the numerator and denominator of the rational `r`.
"""
bc_degree(r::Barycentric) = length(r.nodes) - 1

"""
    bc_poles(r::Barycentric)

Return the poles of the rational function `r`.
"""
function bc_poles(r::Barycentric)
    T = F64
    w = bc_weights(r)
    nonzero = @. !iszero(w)
    z, w = bc_nodes(r)[nonzero], w[nonzero]
    m = length(w)
    B = diagm( [zero(T); ones(T, m)] )
    E = [zero(T) transpose(w); ones(T, m) diagm(z) ];

    pol = [] # Put it into scope
    try
        pol = filter( isfinite, eigvals(E, B) )
    catch
        # Generalized eigen not available in extended precision, so:
        λ = filter( z->abs(z)>1e-13, eigvals(E\B) )
        pol = 1 ./ λ
    end

    return pol
end

"""
    (r::Barycentric)(z::Number)

Evaluate the rational function at `z`.
"""
function (r::Barycentric)(z::Number)
    if isinf(z)
        return sum(r.w_times_f) / sum(r.weights)
    end
    #
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

Adaptively compute a rational interpolant.

### Arguments
* `z::AbstractVector{<:Number}` -> Interpolation nodes.
* `y::AbstractVector{<:Number}` -> Values at nodes.
* `max_degree::Integer=150` -> Maximum numerator/denominator degree to use.
* `float_type::Type=Float64` -> Floating point type to use for the computation.
* `tol::Real=1000*eps(float_type)` -> Tolerance for stopping.
* `lookahead::Integer=10` -> Number of iterations to determines stagnation.
* `stats::Bool=false` -> Return convergence statistics.

### Returns
* `r::Barycentric` -> The rational interpolant.
* `stats::NamedTuple` -> Convergence statistics, if keyword `stats=true`.

### Examples
```julia-repl
julia> z = 1im * range(-10, 10, 500);

julia> y = @. exp(z);

julia> r = aaa(z, y);

julia> degree(r)   # both numerator and denominator
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
    float_type = Float64,
    tol = 1000 * eps(float_type),
    lookahead = 10,
    stats = false
    )

    @assert float_type <: AbstractFloat
    T = float_type
    fmax = norm(y, Inf) # for scaling
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
    r = Barycentric(z[idx], y[idx], w)
    if stats
        return r, (;err, iteration)
    else
        return r
    end
end

mutable struct PronyApproximation
    𝑁ₚ :: Int64
    ωₚ :: Vector{Float64}
    𝐺ₚ :: Vector{ComplexF64}
    Γₚ :: Vector{ComplexF64}
    Ωₚ :: Vector{ComplexF64}
end

function PronyApproximation(winp, Ginp, err)
    N, w, G = prony_data(winp, Ginp)

    S, V = prony_svd(N, G)
    
    v = prony_v(S, V, err)

    cutoff = 1.0 + 0.5 / N
    gamma = prony_gamma(v, cutoff)
    omega = prony_omega(G, gamma)

    idx_sort = sortperm(abs.(omega))
    reverse!(idx_sort)
    omega = omega[idx_sort]
    gamma = gamma[idx_sort]

    return PronyApproximation(N, w, G, gamma, omega)
end

function prony_data(w, G)
    #data = readdlm("giw.data")
    #w = data[:,1]
    #gre = data[:,2]
    #gim = data[:,3]
    #G = gre + gim * im

    osize = length(w)
    nsize = iseven(osize) ? osize - 1 : osize
    N_ = div(nsize, 2)
    w_ = w[1:nsize]
    G_ = G[1:nsize]

    return N_, w_, G_
end

function prony_svd(N, G)
    H = zeros(ComplexF64, N + 1, N + 1)

    for i = 1 : N + 1
        H[i,:] = G[i:i+N]
    end

    _, S, V = svd(H)

    return S, V
end

function prony_v(S, V, err)
    idx = 1
    for i in eachindex(S)
        if S[i] < err
            idx = i
            break
        end
    end

    if S[idx] >= err
        @error "err is set to be too small!"
    end

    v = V[:, idx]
    return reverse!(v)
end

function prony_gamma(u, cutoff)
    non_zero = findall(!iszero, u)
    trailing_zeros = length(u) - non_zero[end]
    unew = u[non_zero[1]:non_zero[end]]
    N = length(unew)
    if N > 1
        A = diagm(-1=>ones(ComplexF64, N - 2))
        @. A[1,:] = -unew[2:end] / unew[1]
        roots = eigvals(A)
    else
    end
    gamma = vcat(roots, zeros(ComplexF64, trailing_zeros))
    filter!(x -> abs(x) < cutoff, gamma)
    return gamma
end

function prony_omega(G, gamma)
    A = zeros(ComplexF64, length(G), length(gamma))
    for i = 1:length(G)
        A[i,:] = gamma .^ (i - 1)
    end
    return pinv(A) * G
end

function (pa::PronyApproximation)(w::Vector{Float64})
    x0 = @. (w - w[1]) / (w[end] - w[1])
    A = zeros(ComplexF64, length(x0), length(pa.Ωₚ))
    for i in eachindex(x0)
        @. A[i,:] = pa.Γₚ ^ (2.0 * pa.𝑁ₚ * x0[i])
    end
    return A * pa.Ωₚ
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
"""
mutable struct  BarRatContext
    Gᵥ   :: Vector{C64}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    ℬ    :: Union{Missing,Barycentric}
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

    return BarRatContext(Gᵥ, grid, mesh, missing)
end

function run(brc::BarRatContext)
    err = 1.0e-12
    pa = PronyApproximation(brc.grid.ω, brc.Gᵥ, err)
    @show pa(brc.grid.ω) .- brc.Gᵥ
    #r = aaa(brc.grid.ω * im, brc.Gᵥ)
    r = aaa(brc.grid.ω * im, pa(brc.grid.ω))
    @show bc_poles(r)
    brc.ℬ = r
end

function last(brc::BarRatContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

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

