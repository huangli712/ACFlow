#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/21
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
        som_update(𝑆, ω, 𝐺)
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

    for k = 1:_Know
        #c = ommin + wmin / 2.0 + (ommax - ommin - wmin) * rand(𝑆.rng, F64)
        #w = wmin + (min(2 * (c - ommin), 2 * (ommax - c)) - wmin) * rand(𝑆.rng, F64)
        #h = weight[k] / w
        push!(𝑆.att_conf, Rectangle(h[k], w[k], c[k]))
        calc_dev_rec(Rectangle(h[k], w[k], c[k]), k, 𝑆.att_elem_dev, ω)
    end
    𝑆.att_dev = calc_dev(𝑆.att_elem_dev, _Know, 𝐺)
    #@show att_dev
end

function som_update(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
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

    @show 𝑆.tmp_conf
    _som_change_weight(𝑆, ω, 𝐺)
    error()

    for i = 1:T1
        𝑆.dacc = d1
        update_type = rand(𝑆.rng, 1:7)

        @cswitch update_type begin
            @case 1
                _som_add(𝑆, ω, 𝐺)
                break

            @case 2
                _som_remove(𝑆, ω, 𝐺)
                break

            @case 3
                _som_shift(𝑆, ω, 𝐺)
                break

            @case 4
                _som_change_width(𝑆, ω, 𝐺)
                break

            @case 5
                _som_change_weight(𝑆, ω, 𝐺)
                break

            @case 6
                _som_split(𝑆, ω, 𝐺)
                break

            @case 7
                _som_merge(𝑆, ω, 𝐺)
                break
        end
    end

    for j = T1+1:Tmax
        𝑆.dacc = d2
        update_type = rand(𝑆.rng, 1:7)

        @cswitch update_type begin
            @case 1
                _som_add(𝑆, ω, 𝐺)
                break

            @case 2
                _som_remove(𝑆, ω, 𝐺)
                break

            @case 3
                _som_shift(𝑆, ω, 𝐺)
                break

            @case 4
                _som_change_width(𝑆, ω, 𝐺)
                break

            @case 5
                _som_change_weight(𝑆, ω, 𝐺)
                break

            @case 6
                _som_split(𝑆, ω, 𝐺)
                break

            @case 7
                _som_merge(𝑆, ω, 𝐺)
                break
        end
    end
end

function _som_add(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    println("add Rectangle")
    smin = P_SOM["smin"]
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]
    t = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t = 23
    if 𝑆.tmp_conf[t].h * 𝑆.tmp_conf[t].w ≤ 2.0 * smin
        return
    end
    #@show 𝑆.tmp_conf[t].h , 𝑆.tmp_conf[t].w

    dx_min = smin
    dx_max = 𝑆.tmp_conf[t].h * 𝑆.tmp_conf[t].w - smin
    #@show dx_min, dx_max
    if dx_max ≤ dx_min
        return
    end

    c = (ommin + wmin / 2.0) + (ommax - ommin - wmin) * rand(𝑆.rng, F64)
    c = -1.68255 # <----
    w_new_max = 2.0 * min(ommax - c, c - ommin)
    #@show c , w_new_max
    dx = Pdx(dx_min, dx_max, γ, 𝑆.rng)
    #@show dx

    r = rand(𝑆.rng, F64)
    r = 0.125254
    𝑆.new_conf = copy(𝑆.tmp_conf)
    𝑆.new_elem_dev = copy(𝑆.elem_dev)
    h = dx / w_new_max + (dx / wmin - dx / w_new_max) * r
    w = dx / h
    #@show c, h, w
    push!(𝑆.new_conf, Rectangle(h, w, c))
    𝑆.new_conf[t].h = 𝑆.new_conf[t].h - dx / 𝑆.new_conf[t].w
    #@show 𝑆.new_conf
    calc_dev_rec(𝑆.new_conf[t], t, 𝑆.new_elem_dev, ω)
    #@show 𝑆.new_conf[t]
    calc_dev_rec(𝑆.new_conf[end], length(𝑆.new_conf), 𝑆.new_elem_dev, ω)
    #@show 𝑆.new_conf[end]
    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)
    #@show 𝑆.new_dev

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = copy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = copy(𝑆.new_elem_dev)
        𝑆.accepted_steps[1] = 𝑆.accepted_steps[1] + 1
    end
    𝑆.trial_steps[1] = 𝑆.trial_steps[1] + 1

    #@show length(𝑆.tmp_conf)
end

function _som_remove(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    println("remove Rectangle")

    t1 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t2 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t1 = 23
    t2 = 25
    if t1 == t2
        t2 = (t1 + 1) % length(𝑆.tmp_conf)
    end

    _conf_size = length(𝑆.tmp_conf)
    dx = 𝑆.tmp_conf[t1].h * 𝑆.tmp_conf[t1].w
    #@show dx

    𝑆.new_conf = copy(𝑆.tmp_conf)
    𝑆.new_elem_dev = copy(𝑆.elem_dev)
    𝑆.new_conf[t2].h = 𝑆.new_conf[t2].h + dx / 𝑆.new_conf[t2].w
    𝑆.new_conf[t1] = 𝑆.new_conf[end]
    pop!(𝑆.new_conf)

    #@show 𝑆.new_conf
    if t1 < _conf_size
        calc_dev_rec(𝑆.new_conf[t1], t1, 𝑆.new_elem_dev, ω)
    end

    if t2 < _conf_size
        calc_dev_rec(𝑆.new_conf[t2], t2, 𝑆.new_elem_dev, ω)
    end

    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)
    #@show 𝑆.new_dev

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = copy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = copy(𝑆.new_elem_dev)
        𝑆.accepted_steps[2] = 𝑆.accepted_steps[2] + 1
    end
    𝑆.trial_steps[2] = 𝑆.trial_steps[2] + 1
    #@show length(𝑆.tmp_conf)
end

function _som_shift(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    println("shift Rectangle")
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t = 23

    dx_min = ommin + 𝑆.tmp_conf[t].w / 2.0 - 𝑆.tmp_conf[t].c
    dx_max = ommax - 𝑆.tmp_conf[t].w / 2.0 - 𝑆.tmp_conf[t].c
    if dx_max ≤ dx_min
        return
    end
    #@show dx_min, dx_max

    dc = Pdx(dx_min, dx_max, γ, 𝑆.rng)
    #@show dc

    _conf_size = length(𝑆.tmp_conf)
    𝑆.new_conf = copy(𝑆.tmp_conf)
    𝑆.new_elem_dev = copy(𝑆.elem_dev)
    𝑆.new_conf[t].c = 𝑆.new_conf[t].c + dc

    calc_dev_rec(𝑆.new_conf[t], t, 𝑆.new_elem_dev, ω)
    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)
    #@show 𝑆.new_dev

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = copy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = copy(𝑆.new_elem_dev)
        𝑆.accepted_steps[3] = 𝑆.accepted_steps[3] + 1
    end
    𝑆.trial_steps[3] = 𝑆.trial_steps[3] + 1
    #@show length(𝑆.tmp_conf)
end

function _som_change_width(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    println("change width of Rectangle")
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t = 23

    weight = 𝑆.tmp_conf[t].h * 𝑆.tmp_conf[t].w
    dx_min = wmin - 𝑆.tmp_conf[t].w
    dx_max = min(2 * (𝑆.tmp_conf[t].c - ommin), 2 * (ommax - 𝑆.tmp_conf[t].c)) - 𝑆.tmp_conf[t].w
    if dx_max ≤ dx_min
        return
    end
    dw = Pdx(dx_min, dx_max, γ, 𝑆.rng)
    #@show weight, dx_min, dx_max, dw

    _conf_size = length(𝑆.tmp_conf)
    𝑆.new_conf = copy(𝑆.tmp_conf)
    𝑆.new_elem_dev = copy(𝑆.elem_dev)
    𝑆.new_conf[t].w = 𝑆.new_conf[t].w + dw
    𝑆.new_conf[t].h = weight / 𝑆.new_conf[t].w
    calc_dev_rec(𝑆.new_conf[t], t, 𝑆.new_elem_dev, ω)

    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)
    #@show 𝑆.new_dev

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = copy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = copy(𝑆.new_elem_dev)
        𝑆.accepted_steps[4] = 𝑆.accepted_steps[4] + 1
    end
    𝑆.trial_steps[4] = 𝑆.trial_steps[4] + 1
    #@show length(𝑆.tmp_conf)
end

function _som_change_weight(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    println("change weight of Rectangle")
end

function _som_split(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    println("split Rectangle")
end

function _som_merge(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
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

function calc_kappa()
end

function Pdx(xmin::F64, xmax::F64, γ::F64, rng::AbstractRNG)
    _X = max(abs(xmin), abs(xmax))
    _lambda = γ / _X
    _elx = exp(-1 * _lambda * abs(xmin))
    _N = _lambda / ((xmin / abs(xmin)) * (exp(-1 * _lambda * abs(xmin)) - 1)
        + (xmax / abs(xmax)) * (1 - exp(-1 * _lambda * abs(xmax))))
 
    y = rand(rng, F64)
    y = 0.415661
    _lysn = _lambda * y / _N
    if xmin ≥ 0
        return -1 * log(_elx - _lysn) / _lambda
    elseif xmax ≤ 0
        return log(_lysn + _elx) / _lambda
    else
        _C1 = _N * (1 - _elx) / _lambda
        if y <= _C1
            return log(_lysn + _elx) / _lambda
        else
            return -1 * log(1 - _lysn + _lambda * _C1 / _N) / _lambda
        end
    end
end
