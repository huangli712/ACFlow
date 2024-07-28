using DelimitedFiles
using LinearAlgebra
using Statistics

function prony_data()
    data = readdlm("giw.data")
    w = data[:,1]
    gre = data[:,2]
    gim = data[:,3]
    G = gre + gim * im

    osize = length(w)
    nsize = iseven(osize) ? osize - 1 : osize
    N = div(nsize, 2)
    w = w[1:nsize]
    G = G[1:nsize]

    return N, w, G
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

mutable struct PronyApproximation
    𝑁ₚ :: Int64
    ωₚ :: Vector{Float64}
    𝐺ₚ :: Vector{ComplexF64}
    Γₚ :: Vector{ComplexF64}
    Ωₚ :: Vector{ComplexF64}
end

function PronyApproximation(err)
    N, w, G = prony_data()

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

function get_value(omega, gamma, w, N)
    x0 = @. (w - w[1]) / (w[end] - w[1])
    A = zeros(ComplexF64, length(x0), length(omega))
    for i in eachindex(x0)
        @. A[i,:] = gamma ^ (2.0 * N * x0[i])
    end
    return A * omega
end

function (pa::PronyApproximation)(w::Vector{Float64})
    x0 = @. (w - w[1]) / (w[end] - w[1])
    A = zeros(ComplexF64, length(x0), length(pa.Ωₚ))
    for i in eachindex(x0)
        @. A[i,:] = pa.Γₚ ^ (2.0 * pa.𝑁ₚ * x0[i])
    end
    return A * pa.Ωₚ
end

err = 1.0e-3
pa = PronyApproximation(err)
value = pa(pa.ωₚ) #get_value(pa.Ωₚ, pa.Γₚ, pa.ωₚ, pa.𝑁ₚ)
@show maximum(abs.(pa.𝐺ₚ - value))
@show mean(abs.(pa.𝐺ₚ - value))
