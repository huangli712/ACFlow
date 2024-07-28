using DelimitedFiles
using LinearAlgebra
using Statistics

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

function get_svd(N, w, G)
    a = w[1]
    b = w[end]
    x_k = range(a, b, 2 * N + 1)
    H = zeros(ComplexF64, N + 1, N + 1)

    for i = 1 : N + 1
        H[i,:] = G[i:i+N]
    end

    _, S, V = svd(H)

    return a, b, x_k, S, V
end

function find_idx_with_err(S, err)
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

    return idx
end

function find_v_with_idx(S, V, idx)
    if idx >= length(S)
        @error "index is invalid!"
    end

    sigma = S[idx]
    v = V[:, idx]
    return sigma, v
end

function roots(u)
    non_zero = findall(!iszero, u)
    trailing_zeros = length(u) - non_zero[end]
    #@show trailing_zeros
    unew = u[non_zero[1]:non_zero[end]]
    #@show length(unew)
    N = length(unew)
    if N > 1
        A = diagm(-1=>ones(ComplexF64, N - 2))
        #println(A[end,:])
        @. A[1,:] = -unew[2:end] / unew[1]
        #for i = 1:49
        #    println(i, " ", A[1,i])
        #end
        #println(A[1,:])
        #println(A[end,:])
        #println(size(A))
        roots = eigvals(A)
        #@show size(roots)
    else
    end
    return vcat(roots, zeros(ComplexF64, trailing_zeros))
end

function find_omega(G, gamma)
    #@show size(G), size(gamma)
    A = zeros(ComplexF64, length(G), length(gamma))
    #@show size(A)
    for i = 1:length(G)
        A[i,:] = gamma .^ (i - 1)
    end
    return pinv(A) * G
end

function get_value(omega, gamma, x, a, b, N)
    x0 = @. (x - a) / (b - a)
    #@show x0 
    A = zeros(ComplexF64, length(x0), length(omega))
    for i = 1:length(x0)
        @. A[i,:] = gamma ^ (2.0 * N * x0[i])
    end
    #@show A[1,:]
    #@show A[end,:]
    value = A * omega
end

err = 1.0e-3
a, b, x_k, S, V = get_svd(N, w, G)
idx = find_idx_with_err(S, err)
sigma, v = find_v_with_idx(S, V, idx)
#@show sigma
#@show v

cutoff = 1.0 + 0.5 / N
#@show cutoff

#using Roots
vinv = reverse(v)
#@show v
#@show vinv

gamma = roots(vinv)
#@show gamma
filter!(x -> abs(x) < cutoff, gamma)
#@show length(gamma)
#@show gamma

omega = find_omega(G, gamma)

#@show gamma
#@show omega

idx_sort = sortperm(abs.(omega))
#@show idx_sort
reverse!(idx_sort)
omega = omega[idx_sort]
gamma = gamma[idx_sort]
#@show omega
#@show gamma
value = get_value(omega, gamma, x_k, a, b, N)
@show length(value)
@show length(G)
@show maximum(abs.(G - value))
@show mean(abs.(G - value))
