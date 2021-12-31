#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/31
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
    @show value / g0
    #@show error
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