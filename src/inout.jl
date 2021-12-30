#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/30
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
        
    end
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