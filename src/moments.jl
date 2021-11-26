#
# Project : Pansy
# Source  : ZenCore.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/19
#

function calc_moments(ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    niw = length(𝐺.value)

    N_FIT_MAX = 200
    
    n_c = 3
    n_min = 2
    n_max = niw - 2 * n_c  - 4

    V𝑀₀ = F64[]
    V𝑀₁ = F64[]
    V𝑀₂ = F64[]
    V𝑀₃ = F64[]
    for n = n_min:n_max
        n_fit = N_FIT_MAX
        if niw - n + 1 < N_FIT_MAX
            n_fit = niw - n + 1
        end

        𝐶 = diagm(𝐺.covar)[2 * n - 1 : 2 * (n + n_fit - 1), 2 * n - 1 : 2 * (n + n_fit - 1)]
        𝑋 = zeros(F64, 2 * n_fit, 2 * n_c)
        for j = 1:n_c
            for i = n:(n + n_fit - 1)
                𝑋[2 * (i - n) + 1, 2 * j - 0] = (-1)^j / (ω.grid[i])^(2*j)
                𝑋[2 * (i - n) + 2, 2 * j - 1] = (-1)^j / (ω.grid[i])^(2*j-1)
            end 
        end
        𝐴 = 𝑋' * inv(𝐶) * 𝑋

        𝐻 = F64[]
        for i = n:(n + n_fit - 1)
            v = 𝐺.value[i]
            push!(𝐻, real(v))
            push!(𝐻, imag(v))
        end
        𝐵 = 𝑋' * inv(𝐶) * 𝐻

        𝑉 = 𝐴 \ 𝐵

        push!(V𝑀₀, 𝑉[1])
        push!(V𝑀₁, 𝑉[2])
        push!(V𝑀₂, 𝑉[3])
        push!(V𝑀₃, 𝑉[4])
    end

    var𝑀₀ = F64[]
    var𝑀₁ = F64[]
    var𝑀₂ = F64[]
    var𝑀₃ = F64[]

    n_v = Int(niw / 16)
    if n_v < 2
        n_v = 2
    end

    for j = n_v : n_max - n_min - n_v
        push!(var𝑀₀, var(V𝑀₀[j - n_v + 1: j + n_v + 1]))
        push!(var𝑀₁, var(V𝑀₁[j - n_v + 1: j + n_v + 1]))
        push!(var𝑀₂, var(V𝑀₂[j - n_v + 1: j + n_v + 1]))
        push!(var𝑀₃, var(V𝑀₃[j - n_v + 1: j + n_v + 1]))
    end

    _, j₀ = findmin(var𝑀₀)
    _, j₁ = findmin(var𝑀₁)
    _, j₂ = findmin(var𝑀₂)
    _, j₃ = findmin(var𝑀₃)

    j₀ = j₀ + n_v
    j₁ = j₁ + n_v
    j₂ = j₂ + n_v
    j₃ = j₃ + n_v

    𝑀₀ = mean(V𝑀₀[j₀ - n_v:j₀ + n_v])
    𝑀₁ = mean(V𝑀₁[j₁ - n_v:j₁ + n_v])
    𝑀₂ = mean(V𝑀₂[j₂ - n_v:j₂ + n_v])
    𝑀₃ = mean(V𝑀₃[j₃ - n_v:j₃ + n_v])

    return MomentsData(𝑀₀, 𝑀₁, 𝑀₂, 𝑀₃)
end