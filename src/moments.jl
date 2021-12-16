#
# Project : Gardenia
# Source  : moments.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/16
#

function calc_moments(ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    niw = length(𝐺.value)

    N_FIT_MAX = 200
    N_FIT_FIN = 300
    
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
        push!(var𝑀₀, var(V𝑀₀[j - n_v + 1 : j + n_v + 1]))
        push!(var𝑀₁, var(V𝑀₁[j - n_v + 1 : j + n_v + 1]))
        push!(var𝑀₂, var(V𝑀₂[j - n_v + 1 : j + n_v + 1]))
        push!(var𝑀₃, var(V𝑀₃[j - n_v + 1 : j + n_v + 1]))
    end

    _, j₀ = findmin(var𝑀₀)
    _, j₁ = findmin(var𝑀₁)
    _, j₂ = findmin(var𝑀₂)
    _, j₃ = findmin(var𝑀₃)

    j₀ = j₀ + n_v
    j₁ = j₁ + n_v
    j₂ = j₂ + n_v
    j₃ = j₃ + n_v

    𝑀₀ = mean(V𝑀₀[j₀ - n_v : j₀ + n_v])
    𝑀₁ = mean(V𝑀₁[j₁ - n_v : j₁ + n_v])
    𝑀₂ = mean(V𝑀₂[j₂ - n_v : j₂ + n_v])
    𝑀₃ = mean(V𝑀₃[j₃ - n_v : j₃ + n_v])

    #@show j₀ - n_v, j₀ + n_v
    #@show 𝑀₀, 𝑀₁, 𝑀₂, 𝑀₃
    
    j = n_v + 1
    #@show abs(mean(V𝑀₀[j - n_v : j + n_v]))
    while j < j₀ && ( abs(mean(V𝑀₀[j - n_v : j + n_v]) - 𝑀₀) / 𝑀₀ > 0.002 || 
                      std(V𝑀₀[j - n_v : j + n_v]) / 𝑀₀ > 0.002 )
        #@show j, j₀, abs(mean(V𝑀₀[j - n_v : j + n_v]) - 𝑀₀) / 𝑀₀, std(V𝑀₀[j - n_v : j + n_v]) / 𝑀₀
        j = j + 1
    end
    #@show j
    ωc = j + n_min - 2
    #@show j, ω.grid[j], ω.grid[niw]



    n = j₀ - 1 + n_min - 1
    n_fit = N_FIT_FIN
    if n_fit > niw - n + 1
        n_fit = niw - n + 1
    end
    #@show n_fit, n_min, j₀ # 27 2 37
    #error()
    𝐶 = diagm(𝐺.covar)[2 * n - 1 : 2 * (n + n_fit - 1), 2 * n - 1 : 2 * (n + n_fit - 1)]
    𝑋 = zeros(F64, 2 * n_fit, 2 * n_c)
    for j = 1:n_c
        for i = n:(n + n_fit - 1)
            𝑋[2 * (i - n) + 1, 2 * j - 0] = (-1)^j / (ω.grid[i])^(2*j)
            𝑋[2 * (i - n) + 2, 2 * j - 1] = (-1)^j / (ω.grid[i])^(2*j-1)
        end
    end
    𝐴 = 𝑋' * inv(𝐶) * 𝑋
    𝐴 = (𝐴 + 𝐴') ./ 2.0
    𝐶𝑀 = (inv(𝐴))[1:4,1:4]
    #@show 𝐶𝑀

    𝐶𝑀[1,:] .= 0.0
    𝐶𝑀[:,1] .= 0.0
    𝐶𝑀[1,1] = 1.0E-4 ^ 2.0
    #@show 𝐶𝑀
    #error() 

    return ωc, MomentsData(𝑀₀, 𝑀₁, 𝑀₂, 𝑀₃, 𝐶𝑀)
end

function trunc_data!(ωc::I64, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    if ωc < 20
        ωc = 20
    end

    resize!(ω.grid, ωc)

    #@show ωc, ω.grid[ωc]
    #@show ω.grid

    resize!(𝐺.value, ωc)
    resize!(𝐺.error, ωc)
    resize!(𝐺.covar, ωc)
    #@show 𝐺.value
end

function diag_covar(𝑀::MomentsData)
    println("here")
    #@show 𝑀.𝐶𝑀
    #=
    𝑀.𝐶𝑀 = [1.0000e-08        0e+00        0e+00        0e+00;
            0e+00   3.8163e-11        0e+00   1.6561e-08;
            0e+00        0e+00   7.9080e-05        0e+00;
            0e+00   1.6561e-08        0e+00   7.3933e-06]
    =#
    𝐹 = eigen(𝑀.𝐶𝑀[2:4,2:4])
    #@show 𝑀.𝐶𝑀[2:4,2:4]
    #@show 𝐹.values
    #@show 𝐹.vectors
    #error()

    VM = zeros(F64,4,4)
    VM[1,1] = 1.0
    VM[2:4,2:4] .= 𝐹.vectors

    WM = diagm(1.0 ./ insert!(sqrt.(𝐹.values), 1, sqrt(𝑀.𝐶𝑀[1,1])))
    @show VM
    @show WM
end
