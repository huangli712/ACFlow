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

function prony_approx(err)
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

    return N, w, G, omega, gamma
end

function get_value(omega, gamma, w, N)
    x0 = @. (w - w[1]) / (w[end] - w[1])
    A = zeros(ComplexF64, length(x0), length(omega))
    for i in eachindex(x0)
        @. A[i,:] = gamma ^ (2.0 * N * x0[i])
    end
    return A * omega
end

err = 1.0e-3
N, w, G, omega, gamma = prony_approx(err)
value = get_value(omega, gamma, w, N)
@show maximum(abs.(G - value))
@show mean(abs.(G - value))
