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
    rng :: AbstractRNG

    dev :: Vector{F64}
    conf :: Vector{Vector{Rectangle}}

    att_conf :: Vector{Rectangle}
    tmp_conf :: Vector{Rectangle}
    new_conf :: Vector{Rectangle}

    att_elem_dev :: Array{C64,2}
    elem_dev :: Array{C64,2}
    new_elem_dev :: Array{C64,2}

    trial_steps :: Vector{I64}
    accepted_steps :: Vector{I64}

    att_dev :: F64
    tmp_dev :: F64
    new_dev :: F64
    dacc    :: F64
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
    rng = MersenneTwister(2345)

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

    att_elem_dev = zeros(C64, Ngrid, Kmax)
    elem_dev = zeros(C64, Ngrid, Kmax)
    new_elem_dev = zeros(C64, Ngrid, Kmax)

    trial_steps = zeros(I64, 7)
    accepted_steps = zeros(I64, 7)
    #@show size(conf)
    #@show typeof(conf[1])

    return T_SOM(rng,
                 dev, 
                 conf, 
                 att_conf, 
                 tmp_conf, 
                 new_conf, 
                 att_elem_dev, 
                 elem_dev, 
                 new_elem_dev, 
                 trial_steps, 
                 accepted_steps, 0.0, 0.0, 0.0, 0.0)
end

function som_run(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Lmax = P_SOM["Lmax"]
    for l = 1:Lmax
        #@show rand(𝑆.rng, F64)
        som_try(𝑆, ω, 𝐺)
    end
end

function som_try(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Nf = P_SOM["Nf"]
    #println("here")
    som_random(𝑆, ω, 𝐺)
    #error()
    for f = 1:Nf
        som_update(𝑆::T_SOM)
    end
end

function som_random(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    smin = P_SOM["smin"]
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    _Know = 25
    _weight = zeros(F64, _Know)
    _weight = [
        0.139286,
        0.16858,
        0.265188,
        0.548449,
        0.350178,
        0.307434,
        0.34026,
        0.80599,
        0.715488,
        0.0402386,
        0.543467,
        0.71602,
        0.631526,
        0.580398,
        0.578101,
        0.917187,
        0.949232,
        0.74917,
        0.959639,
        0.245511,
        0.133757,
        0.0400198,
        0.308264,
        0.762588,
        1.0
    ]
    sort!(_weight)
    weight = diff(_weight)
    insert!(weight, 1, _weight[1])
    #@show weight
    sort!(weight)
    #@show weight

    plus_count = 1
    minus_count = _Know
    while weight[plus_count] < smin
        while weight[minus_count] < 2 * smin
            minus_count = minus_count - 1
        end
        weight[plus_count] = weight[plus_count] + smin
        weight[minus_count] = weight[minus_count] - smin
        plus_count = plus_count + 1
    end
    #@show weight

    empty!(𝑆.att_conf)
    fill!(𝑆.att_elem_dev, zero(C64))
    #@show size(𝑆.att_conf)

    #=
    c = [
        0.437433,
4.35723,
-4.86402,
-8.52167,
8.35216,
0.131443,
4.66806,
-5.48925,
2.18855,
-7.10065,
-8.64211,
-7.61222,
-1.01904,
-7.64489,
-6.92268,
-0.296024,
2.01579,
5.21624,
-0.430991,
1.5215,
-7.92062,
7.95186,
1.16635,
-4.30618,
5.27789
    ]
    =#
#=
    w = [
        3.80198,
9.33006,
6.98871,
0.282352,
0.58726,
6.56658,
8.78425,
2.47927,
1.9538,
2.65317,
1.9506,
4.59944,
16.0818,
4.70633,
2.38606,
16.4445,
2.18616,
2.21616,
8.25936,
13.7515,
1.32354,
1.37565,
13.4274,
2.7657,
4.68898
    ]
=#
    #=
    h = [
        0.00137267,
0.000592927,
0.000834196,
0.0258451,
0.0169977,
0.000841944,
0.00112916,
0.00419728,
0.00686794,
0.00741652,
0.0150182,
0.00644669,
0.00198957,
0.00680902,
0.0138929,
0.00243363,
0.0184623,
0.0190627,
0.00525481,
0.00371801,
0.0581247,
0.0610346,
0.00696474,
0.0402057,
0.0358903
    ]
    =#
    for k = 1:_Know
        c = ommin + wmin / 2.0 + (ommax - ommin - wmin) * rand(𝑆.rng, F64)
        w = wmin + (min(2 * (c - ommin), 2 * (ommax - c)) - wmin) * rand(𝑆.rng, F64)
        h = weight[k] / w
        push!(𝑆.att_conf, Rectangle(h, w, c))
        calc_dev_rec(Rectangle(h, w, c), k, 𝑆.att_elem_dev, ω)
    end
    𝑆.att_dev = calc_dev(𝑆.att_elem_dev, _Know, 𝐺)
    #@show att_dev
end

function som_update(𝑆::T_SOM)
    Tmax = P_SOM["Tmax"]
    dmax = P_SOM["dmax"]
    T1 = rand(𝑆.rng, 1:Tmax)
    #for i = 1:100
    #    @show rand(𝑆.rng, 1:Tmax)
    #end

    d1 = rand(𝑆.rng, F64)
    d2 = 1.0 + (dmax - 1.0) * rand(𝑆.rng, F64)

    𝑆.tmp_conf = copy(𝑆.att_conf)
    𝑆.tmp_dev = 𝑆.att_dev
    𝑆.elem_dev = copy(𝑆.att_elem_dev)

    for i = 1:T1
        𝑆.dacc = d1
        update_type = rand(𝑆.rng, 1:7)

        @cswitch update_type begin
            @case 1
                _som_add()
                break

            @case 2
                _som_remove()
                break

            @case 3
                _som_shift()
                break

            @case 4
                _som_change_width()
                break

            @case 5
                _som_change_weight()
                break

            @case 6
                _som_split()
                break

            @case 7
                _som_merge()
                break
        end
    end

    for j = T1+1:Tmax
        𝑆.dacc = d2
        update_type = rand(𝑆.rng, 1:7)

        @cswitch update_type begin
            @case 1
                _som_add()
                break

            @case 2
                _som_remove()
                break

            @case 3
                _som_shift()
                break

            @case 4
                _som_change_width()
                break

            @case 5
                _som_change_weight()
                break

            @case 6
                _som_split()
                break

            @case 7
                _som_merge()
                break
        end
    end
end

function _som_add()
    println("add Rectangle")
end

function _som_remove()
    println("remove Rectangle")
end

function _som_shift()
    println("shift Rectangle")
end

function _som_change_width()
    println("change width of Rectangle")
end

function _som_change_weight()
    println("change weight of Rectangle")
end

function _som_split()
    println("split Rectangle")
end

function _som_merge()
    println("Merge Rectangle")
end

function calc_dev_rec(r::Rectangle, k::I64, elem_dev::Array{C64,2}, ω::FermionicMatsubaraGrid)
    Ngrid = P_SOM["Ngrid"]

    #@show r.h, r.w, r.c
    for g = 1:Ngrid
        Gs = r.h * log((im * ω.grid[g] - r.c + 0.5 * r.w) / (im * ω.grid[g] - r.c - 0.5 * r.w))
        elem_dev[g,k] = Gs
        #@show g, Gs
    end
    #error()
end

function calc_dev(elem_dev::Array{C64,2}, nk::I64, 𝐺::GreenData)
    Ngrid = P_SOM["Ngrid"]
    
    res = 0.0
    for g = 1:Ngrid
        δ = 0.0
        for k = 1:nk
            δ = δ + elem_dev[g,k]
        end
        res = res + abs((δ - 𝐺.value[g]) / 𝐺.error[g])
    end

    return res
end
