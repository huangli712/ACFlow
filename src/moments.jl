
function calc_moments(ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    niw = length(𝐺.value)

    N_FIT_MAX = 200
    
    n_c = 3
    n_min = 2
    n_max = niw - 2 * n_c  - 4
    @show n_min, n_max

    for n = n_min:n_max
        n_fit = N_FIT_MAX
        if niw - n + 1 < N_FIT_MAX
            n_fit = niw - n + 1
        end

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

        𝐶 = diagm(𝐺.covar)[2 * n - 1 : 2 * (n + n_fit - 1), 2 * n - 1 : 2 * (n + n_fit - 1)]
        #@show diag(𝐶)

        𝐴 = 𝑋' * inv(𝐶) * 𝑋
        

        #@show 𝑋
        error()
    end
end