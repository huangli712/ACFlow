#
# Project : Gardenia
# Source  : inout.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/03
#

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