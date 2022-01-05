#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/03
#

function read_data!(::Type{ImaginaryTimeGrid})
    grid = F64[]
    value = F64[]
    error = F64[]
    covar = F64[]

    nbin = 1000
    ntau = 100

    bin_data = zeros(F64, ntau, nbin)

    open("tau.data", "r") do fin
        readline(fin)
        for i = 1:ntau
            tau = parse(F64, readline(fin))
            push!(grid, tau)
        end
    end
    #@show grid

    open("cor.data", "r") do fin
        readline(fin)
        for b = 1:nbin
            readline(fin)
            for i = 1:ntau
                g = parse(F64, readline(fin))
                bin_data[i,b] = g
            end
        end
    end
    #@show bin_data[:,end]

    # try to calculate mean value
    for i  = 1:ntau
        val = mean(bin_data[i,:])
        push!(value, val)
    end
    #@show value

    seed = rand(1:1000000)
    rng = MersenneTwister(seed)
    nbootstrap = 5000
    bootstrap = zeros(F64, ntau, nbootstrap)
    for i = 1:nbootstrap
        for b = 1:nbin
            ind = rand(rng, 1:nbin)
            @. bootstrap[:,i] = bootstrap[:,i] + bin_data[:,ind]
        end
    end
    bootstrap = bootstrap ./ nbin
    for i = 1:ntau
        err = sum((bootstrap[i,:] .- value[i]) .^ 2.0)
        err = sqrt((err / nbootstrap))
        push!(error, err)
    end
    #@show error

    unselected_tau = findall( x -> abs(x) ≥ 0.1, error[2:end] ./ value[2:end]) .+ 1
    push!(unselected_tau, 1)
    sort!(unselected_tau)
    #@show unselected_tau

    # Filter the data
    g0 = value[1]
    deleteat!(grid, unselected_tau)
    deleteat!(value, unselected_tau)
    deleteat!(error, unselected_tau)
    bootstrap = bootstrap[setdiff(1:end, unselected_tau), :]

    #@show grid
    #@show value
    #@show error

    value = value / g0
    error = error / g0
    bootstrap = bootstrap / g0
    #@show error
    #@show size(value), size(error), size(bootstrap)

    cov_mat_dim = length(value)
    cov_mat = zeros(F64, cov_mat_dim, cov_mat_dim)

    #@show cov_mat_dim

    for i = 1:cov_mat_dim
        for j = 1:cov_mat_dim
            cov_mat[j,i] = sum((bootstrap[i,:] .- value[i]) .* ( bootstrap[j,:] .- value[j]))
        end
    end


    open("test.data", "r") do fin
        for i = 1:cov_mat_dim
            for j = 1:cov_mat_dim
                cov_mat[i,j] = parse(F64, line_to_array(fin)[3])
            end
        end
    end


    #for i = 1:cov_mat_dim
    #    for j = 1:cov_mat_dim
    #        @show i, j, cov_mat[i,j]
    #    end
    #end

    #@show issymmetric(cov_mat)
    #F = eigen(Symmetric(cov_mat), 1:cov_mat_dim)
    #@show F.values
    #@show F.vectors

    eigs, evec = LAPACK.syev!('V', 'U', cov_mat)
    #for i in eachindex(eigs)
    #    @show eigs[i]
    #end

    #@show size(F.vectors), size(value)
    #@show F.vectors * value
    #@show value
    #@show evec' * value

    #value = evec' * value
    
    #@show evec'[1,:]
    #@show evec'[3,:]
    #@show size(evec)

    covar = sqrt(nbootstrap) ./ sqrt.(eigs)
    #@show covar

    return g0, GreenData(value, error, covar), ImaginaryTimeGrid(grid)
end

function read_data!(::Type{FermionicMatsubaraGrid})
    grid  = F64[] 
    value = C64[]
    error = C64[]
    covar = F64[]

    niw = 64
    #
    open("giw.data", "r") do fin
        for i = 1:niw
            arr = parse.(F64, line_to_array(fin))
            push!(grid, arr[1])
            push!(value, arr[2] + arr[3] * im)
        end
    end
    #
    open("err.data", "r") do fin
        for i = 1:niw
            arr = parse.(F64, line_to_array(fin))
            @assert grid[i] == arr[1]
            push!(error, arr[2] + arr[3] * im)
            push!(covar, arr[2]^2)
            push!(covar, arr[3]^2)
        end
    end

    return FermionicMatsubaraGrid(grid), GreenData(value, error, covar)
end

function read_data!(::Type{BosonicMatsubaraGrid})
    error()
end