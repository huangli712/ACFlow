
function calc_moments(ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    niw = length(𝐺.value)

    N_FIT_MAX = 200
    
    n_c = 3
    n_min = 2
    n_max = niw - 2 * n_c  - 4
    #@show n_min, n_max

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
                #println(2 * (i - n) + 0, " ", 2 * j - 1, " ", 𝑋[2 * (i - n) + 1, 2 * j - 0])
                #println(2 * (i - n) + 1, " ", 2 * j - 2, " ", 𝑋[2 * (i - n) + 2, 2 * j - 1])
                #error()
            end 
        end
        #@show diag(𝐶)
        𝐴 = 𝑋' * inv(𝐶) * 𝑋

        𝐻 = F64[]
        for i = n:(n + n_fit - 1)
            v = 𝐺.value[i]
            push!(𝐻, real(v))
            push!(𝐻, imag(v))
            #@show i, real(v)
            #@show i, imag(v)
        end
        𝐵 = 𝑋' * inv(𝐶) * 𝐻

        𝑀 = 𝐴 \ 𝐵

        #@show n
        #@show 𝑀
        #@show 𝑋
        #error()

        push!(V𝑀₀, 𝑀[1])
        push!(V𝑀₁, 𝑀[2])
        push!(V𝑀₂, 𝑀[3])
        push!(V𝑀₃, 𝑀[4])
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
        # @show j
        push!(var𝑀₀, var(V𝑀₀[j - n_v + 1: j + n_v + 1]))
        push!(var𝑀₁, var(V𝑀₁[j - n_v + 1: j + n_v + 1]))
        push!(var𝑀₂, var(V𝑀₂[j - n_v + 1: j + n_v + 1]))
        push!(var𝑀₃, var(V𝑀₃[j - n_v + 1: j + n_v + 1]))
    end

    
end