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

mutable struct T_SOM
    dev :: Vector{F64}
    conf :: Vector{Vector{Rectangle}}

    att_conf :: Vector{Rectangle}
    tmp_conf :: Vector{Rectangle}
    new_conf :: Vector{Rectangle}

    att_elem_dev :: Vector{C64}
    elem_dev :: Vector{C64}
    new_elem_dev :: Vector{C64}

    trial_steps :: Vector{I64}
    accepted_steps :: Vector{I64}
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
    "monitor" => false,
)

function som_init()
    #println("here")
    Lmax = P_SOM["Lmax"]
    Kmax = P_SOM["Kmax"]
    Ngrid = P_SOM["Ngrid"]

    dev = zeros(F64, Lmax)

    conf = []
    for l = 1:Lmax
        _conf = Rectangle[]
        for k = 1:Kmax
            push!(_conf, Rectangle(0.0, 0.0, 0.0))
        end
        push!(conf, _conf)
    end

    att_conf = Rectangle[]
    tmp_conf = Rectangle[]
    new_conf = Rectangle[]
    for k = 1:Kmax
        push!(att_conf, Rectangle(0.0, 0.0, 0.0))
        push!(tmp_conf, Rectangle(0.0, 0.0, 0.0))
        push!(new_conf, Rectangle(0.0, 0.0, 0.0))
    end

    att_elem_dev = zeros(C64, Ngrid * Kmax)
    elem_dev = zeros(C64, Ngrid * Kmax)
    new_elem_dev = zeros(C64, Ngrid * Kmax)

    trial_steps = zeros(I64, 7)
    accepted_steps = zeros(I64, 7)
    #@show size(conf)
    #@show typeof(conf[1])

    return T_SOM(dev, 
                 conf, 
                 att_conf, 
                 tmp_conf, 
                 new_conf, 
                 att_elem_dev, 
                 elem_dev, 
                 new_elem_dev, 
                 trial_steps, 
                 accepted_steps)
end

function som_run(𝑆::T_SOM)
    Lmax = P_SOM["Lmax"]
    for l = 1:Lmax
        som_try()
    end
end

function som_try()
end

function som_random()
end

function som_update()
end

function _som_add()
end

function _som_remove()
end

function _som_split()
end

function _som_merge()
end

function _som_shift()
end

function _som_change()
end
