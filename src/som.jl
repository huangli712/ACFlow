#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/20
#

mutable struct Rectangle
    h :: F64
    w :: F64
    c :: F64
end

const P_SOM = Dict{String, Any}(
    "Lmax" => 400,
    "Ngrid" => 64,
    "Nf" => 2000,
    "Tmax" => 200,
    "Kmax" => 50,
    "nwout" => 100,
    "smin" => 0.005,
    "wmin" => 0.05,
    "gamma" => 2.0,
    "dmax" => 2.0,
    "ommax" => 10.0,
    "ommin" => -10.0,
    "alpha" => 2.0,
    "temp" => 0.05,
    "norm" => -1.0,
    "monitor" => false
)

function som_update()
end

function som_try_shift()
end

function som_try_add()
end

function som_try_remove()
end

function som_try_split()
end

function som_try_merge()
end
