using LinearAlgebra
using DelimitedFiles
using Printf
using Statistics

const I64 = Int64
const F64 = Float64
const C64 = ComplexF64

"""
    PronyApproximation

Mutable struct. Prony approximation to a complex-valued Matsubara function.

### Members

* 𝑁ₚ -> Number of nodes.
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
    PronyApproximation(ω₁, 𝐺₁, ε)

Construct a `PronyApproximation` type interpolant function. Once it is
available, then it can be used to produce a smooth G at ω.

### Arguments
* ω₁::Vector{F64} -> Non-negative Matsubara frequency (raw).
* 𝐺₁::Vector{C64} -> Complex values at ωₚ (raw).
* ε::F64 -> Threshold for the Prony approximation.
"""
function PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, v)
    # Evaluate Γₚ and Ωₚ
    Λ = 1.0 + 0.5 / 𝑁ₚ
    Γₚ = prony_gamma(v, Λ)
    Ωₚ = prony_omega(𝐺ₚ, Γₚ)

    # Sort Γₚ and Ωₚ
    idx_sort = sortperm(abs.(Ωₚ))
    reverse!(idx_sort)
    Ωₚ = Ωₚ[idx_sort]
    Γₚ = Γₚ[idx_sort]

    return PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, Γₚ, Ωₚ)
end

function PronyApproximation(ω₁, 𝐺₁, ε)
    # Get number of nodes, frequency points ωₚ, and Matsubara data 𝐺ₚ.
    𝑁ₚ, ωₚ, 𝐺ₚ = prony_data(ω₁, 𝐺₁)

    # Singular value decomposition
    S, V = prony_svd(𝑁ₚ, 𝐺ₚ)
    v = prony_v(S, V, ε)

    return PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, v)
end

function PronyApproximation(ω₁, 𝐺₁)
    𝑁ₚ, ωₚ, 𝐺ₚ = prony_data(ω₁, 𝐺₁)
    S, V = prony_svd(𝑁ₚ, 𝐺ₚ)

    exp_idx = find_idx_with_exp_decay(S)
    ε = 1000 * S[exp_idx]
    new_idx = findfirst(x -> x < ε, S)

    idx_list = collect(range(new_idx, min(exp_idx + 10, length(S))))
    err_list = zeros(F64, length(idx_list))

    for i in eachindex(idx_list)
        idx = idx_list[i]
        v = V[:,idx]
        reverse!(v)

        Gn = PronyApproximation(𝑁ₚ, ωₚ, 𝐺ₚ, v)(ωₚ)
        err_ave = mean(abs.(Gn - 𝐺ₚ))
        @show i, idx, err_ave
        err_list[i] = err_ave
    end

    idx = idx_list[argmin(err_list)]
    #@show idx, S[idx]
    v = V[:,idx]
    reverse!(v)

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

"""
    prony_v(S, V, ε)

Extract suitable vector `v` from orthogonal matrix `V` according to the
threshold `ε`. The diagonal matrix (singular values) `S` is used to test
whether the threshold `ε` is reasonable and figure out the index for
extracting `v` from `V`.
"""
function prony_v(S, V, ε)
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

function find_idx_with_exp_decay(S::Vector{F64})
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

data = readdlm("w.data")
ω₁ = data[:]

data = readdlm("giw.data")
gre = data[:,1]
gim = data[:,2]
𝐺₁ = gre + gim * im

PronyApproximation(ω₁, 𝐺₁)