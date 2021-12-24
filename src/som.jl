#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/22
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
    "Lmax" => 100,
    "Ngrid" => 64,
    "Nf" => 1000,
    "Tmax" => 100,
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
    seed = rand(1:1000000)
#    seed = 571716
    rng = MersenneTwister(seed)
    @show "seed: ", seed

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
        println("try: $l")
        som_try(l, 𝑆, ω, 𝐺)
        som_output(l, 𝑆)
    end
end

function som_try(l, 𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Nf = P_SOM["Nf"]
    #println("here")
    som_random(𝑆, ω, 𝐺)
    #error()

    for f = 1:Nf
        #println("    update: $f")
        som_update(𝑆, ω, 𝐺)
        #error()
    end
    #error()

    𝑆.dev[l] = 𝑆.att_dev
    𝑆.conf[l] = deepcopy(𝑆.att_conf)
    #@show 𝑆.conf[l]
end

function som_random(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    smin = P_SOM["smin"]
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    Kmax = P_SOM["Kmax"]

    #_Know = 25
    _Know = rand(𝑆.rng, 2:Kmax)
    _weight = zeros(F64, _Know)
    for i = 1:_Know
        _weight[i] = rand(𝑆.rng, F64)
    end
    _weight[end] = 1.0

    #=
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
=#
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
=#

    for k = 1:_Know
        c = ommin + wmin / 2.0 + (ommax - ommin - wmin) * rand(𝑆.rng, F64)
        w = wmin + (min(2.0 * (c - ommin), 2.0 * (ommax - c)) - wmin) * rand(𝑆.rng, F64)
        h = weight[k] / w
        push!(𝑆.att_conf, Rectangle(h, w, c))
        calc_dev_rec(Rectangle(h, w, c), k, 𝑆.att_elem_dev, ω)
        #push!(𝑆.att_conf, Rectangle(h[k], w[k], c[k]))
        #calc_dev_rec(Rectangle(h[k], w[k], c[k]), k, 𝑆.att_elem_dev, ω)
    end
    𝑆.att_dev = calc_dev(𝑆.att_elem_dev, _Know, 𝐺)
    #@show 𝑆.att_dev
    #error()

    norm = 0.0
    for i = 1:length(𝑆.att_conf)
        norm = norm + 𝑆.att_conf[i].h * 𝑆.att_conf[i].w
    end
    @show norm
    #error()
end

function som_update(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Tmax = P_SOM["Tmax"]
    Kmax = P_SOM["Kmax"]
    dmax = P_SOM["dmax"]
    T1 = rand(𝑆.rng, 1:Tmax)

    #for i = 1:100
    #    @show rand(𝑆.rng, 1:Tmax)
    #end

    d1 = rand(𝑆.rng, F64)
    d2 = 1.0 + (dmax - 1.0) * rand(𝑆.rng, F64)

    𝑆.tmp_conf = deepcopy(𝑆.att_conf)
    𝑆.tmp_dev = 𝑆.att_dev
    𝑆.elem_dev = deepcopy(𝑆.att_elem_dev)

    #@show 𝑆.tmp_conf
    #_som_add(𝑆, ω, 𝐺)
    #_som_remove(𝑆, ω, 𝐺)
    #_som_shift(𝑆, ω, 𝐺)
    #_som_change_width(𝑆, ω, 𝐺)
    #_som_change_weight(𝑆, ω, 𝐺)
    #_som_split(𝑆, ω, 𝐺)
    #_som_merge(𝑆, ω, 𝐺)
    #@show 𝑆.tmp_conf
    #@show 𝑆.tmp_dev
    #error()

    for i = 1:T1
        𝑆.dacc = d1
        update_type = rand(𝑆.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(𝑆.tmp_conf) < Kmax - 1
                    _som_add(𝑆, ω, 𝐺)
                end
                break

            @case 2
                if length(𝑆.tmp_conf) > 1
                    _som_remove(𝑆, ω, 𝐺)
                end
                break

            @case 3
                _som_shift(𝑆, ω, 𝐺)
                break

            @case 4
                _som_change_width(𝑆, ω, 𝐺)
                break

            @case 5
                if length(𝑆.tmp_conf) > 1
                    _som_change_weight(𝑆, ω, 𝐺)
                end
                break

            @case 6
                if length(𝑆.tmp_conf) < Kmax - 1
                    _som_split(𝑆, ω, 𝐺)
                end
                break

            @case 7
                if length(𝑆.tmp_conf) > 1
                    _som_merge(𝑆, ω, 𝐺)
                end
                break
        end
    end

    for j = T1+1:Tmax
        𝑆.dacc = d2
        update_type = rand(𝑆.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(𝑆.tmp_conf) < Kmax - 1
                    _som_add(𝑆, ω, 𝐺)
                end
                break

            @case 2
                if length(𝑆.tmp_conf) > 1
                    _som_remove(𝑆, ω, 𝐺)
                end
                break

            @case 3
                _som_shift(𝑆, ω, 𝐺)
                break

            @case 4
                _som_change_width(𝑆, ω, 𝐺)
                break

            @case 5
                if length(𝑆.tmp_conf) > 1
                    _som_change_weight(𝑆, ω, 𝐺)
                end
                break

            @case 6
                if length(𝑆.tmp_conf) < Kmax - 1
                    _som_split(𝑆, ω, 𝐺)
                end
                break

            @case 7
                if length(𝑆.tmp_conf) > 1
                    _som_merge(𝑆, ω, 𝐺)
                end
                break
        end
    end

    #norm = 0.0
    #for i = 1:length(𝑆.tmp_conf)
    #    norm = norm + 𝑆.tmp_conf[i].h * 𝑆.tmp_conf[i].w
    #end

    #@show 𝑆.tmp_dev, 𝑆.att_dev, norm
    if 𝑆.tmp_dev < 𝑆.att_dev
        𝑆.att_conf = deepcopy(𝑆.tmp_conf)
        𝑆.att_dev = 𝑆.tmp_dev
        𝑆.att_elem_dev = deepcopy(𝑆.elem_dev)
    end
    #@show 𝑆.tmp_conf
    #error()
end

function som_output(count::I64, 𝑆::T_SOM)
    println("output")
    alpha = P_SOM["alpha"]
    #Lmax = P_SOM["Lmax"]
    Ngrid = P_SOM["Ngrid"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]

    dev_min = minimum(𝑆.dev[1:count])
    
    Lgood = 0
    Aom = zeros(F64, Ngrid)
    for l = 1:count
        if alpha * dev_min - 𝑆.dev[l] > 0
            Lgood = Lgood + 1
            for w = 1:Ngrid
                _omega = ommin + (w - 1) * (ommax - ommin) / (Ngrid - 1)
                for r = 1:length(𝑆.conf[l])
                    R = 𝑆.conf[l][r]
                    @show l, r, R
                    if R.c - 0.5 * R.w ≤ _omega ≤ R.c + 0.5 * R.w
                        Aom[w] = Aom[w] + R.h
                    end
                end
            end
        end
    end

    @show count, dev_min, Lgood

    if Lgood > 0
        @. Aom = Aom / Lgood
    end

    open("Aw.out", "w") do fout
        for w = 1:Ngrid
            _omega = ommin + (w - 1) * (ommax - ommin) / (Ngrid - 1)
            println(fout, _omega, " ", Aom[w])
        end
    end
end

function _som_add(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    #println("add Rectangle")
    smin = P_SOM["smin"]
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    if 𝑆.tmp_conf[t].h * 𝑆.tmp_conf[t].w ≤ 2.0 * smin
        return
    end

    dx_min = smin
    dx_max = 𝑆.tmp_conf[t].h * 𝑆.tmp_conf[t].w - smin
    if dx_max ≤ dx_min
        return
    end

    c = (ommin + wmin / 2.0) + (ommax - ommin - wmin) * rand(𝑆.rng, F64)
    w_new_max = 2.0 * min(ommax - c, c - ommin)
    dx = Pdx(dx_min, dx_max, γ, 𝑆.rng)

    r = rand(𝑆.rng, F64)
    𝑆.new_conf = deepcopy(𝑆.tmp_conf)
    𝑆.new_elem_dev = deepcopy(𝑆.elem_dev)
    h = dx / w_new_max + (dx / wmin - dx / w_new_max) * r
    w = dx / h

    push!(𝑆.new_conf, Rectangle(h, w, c))
    𝑆.new_conf[t].h = 𝑆.new_conf[t].h - dx / 𝑆.new_conf[t].w
    calc_dev_rec(𝑆.new_conf[t], t, 𝑆.new_elem_dev, ω)
    calc_dev_rec(𝑆.new_conf[end], length(𝑆.new_conf), 𝑆.new_elem_dev, ω)
    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = deepcopy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = deepcopy(𝑆.new_elem_dev)
        𝑆.accepted_steps[1] = 𝑆.accepted_steps[1] + 1
    end
    𝑆.trial_steps[1] = 𝑆.trial_steps[1] + 1
end

function _som_remove(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    t1 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t2 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    while t1 == t2
        t2 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    end

    if t1 < t2
        t1, t2 = t2, t1
    end

    _conf_size = length(𝑆.tmp_conf)
    dx = 𝑆.tmp_conf[t1].h * 𝑆.tmp_conf[t1].w

    𝑆.new_conf = deepcopy(𝑆.tmp_conf)
    𝑆.new_elem_dev = deepcopy(𝑆.elem_dev)


    𝑆.new_conf[t2].h = 𝑆.new_conf[t2].h + dx / 𝑆.new_conf[t2].w
    if t1 < _conf_size
        𝑆.new_conf[t1] = deepcopy(𝑆.new_conf[end])
    else
        @assert t1 == _conf_size
    end
    pop!(𝑆.new_conf)

    if t1 < _conf_size
        calc_dev_rec(𝑆.new_conf[t1], t1, 𝑆.new_elem_dev, ω)
    end
    calc_dev_rec(𝑆.new_conf[t2], t2, 𝑆.new_elem_dev, ω)

    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = deepcopy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = deepcopy(𝑆.new_elem_dev)
        𝑆.accepted_steps[2] = 𝑆.accepted_steps[2] + 1
    end
    𝑆.trial_steps[2] = 𝑆.trial_steps[2] + 1
end

function _som_shift(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(𝑆.rng, 1:length(𝑆.tmp_conf))

    dx_min = ommin + 𝑆.tmp_conf[t].w / 2.0 - 𝑆.tmp_conf[t].c
    dx_max = ommax - 𝑆.tmp_conf[t].w / 2.0 - 𝑆.tmp_conf[t].c
    if dx_max ≤ dx_min
        return
    end

    dc = Pdx(dx_min, dx_max, γ, 𝑆.rng)

    _conf_size = length(𝑆.tmp_conf)
    𝑆.new_conf = deepcopy(𝑆.tmp_conf)
    𝑆.new_elem_dev = deepcopy(𝑆.elem_dev)
    𝑆.new_conf[t].c = 𝑆.new_conf[t].c + dc

    calc_dev_rec(𝑆.new_conf[t], t, 𝑆.new_elem_dev, ω)
    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = deepcopy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = deepcopy(𝑆.new_elem_dev)
        𝑆.accepted_steps[3] = 𝑆.accepted_steps[3] + 1
    end
    𝑆.trial_steps[3] = 𝑆.trial_steps[3] + 1
end

function _som_change_width(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(𝑆.rng, 1:length(𝑆.tmp_conf))

    weight = 𝑆.tmp_conf[t].h * 𝑆.tmp_conf[t].w
    dx_min = wmin - 𝑆.tmp_conf[t].w
    dx_max = min(2.0 * (𝑆.tmp_conf[t].c - ommin), 2.0 * (ommax - 𝑆.tmp_conf[t].c)) - 𝑆.tmp_conf[t].w
    if dx_max ≤ dx_min
        return
    end
    dw = Pdx(dx_min, dx_max, γ, 𝑆.rng)

    _conf_size = length(𝑆.tmp_conf)
    𝑆.new_conf = deepcopy(𝑆.tmp_conf)
    𝑆.new_elem_dev = deepcopy(𝑆.elem_dev)
    𝑆.new_conf[t].w = 𝑆.new_conf[t].w + dw
    𝑆.new_conf[t].h = weight / 𝑆.new_conf[t].w
    calc_dev_rec(𝑆.new_conf[t], t, 𝑆.new_elem_dev, ω)

    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = deepcopy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = deepcopy(𝑆.new_elem_dev)
        𝑆.accepted_steps[4] = 𝑆.accepted_steps[4] + 1
    end
    𝑆.trial_steps[4] = 𝑆.trial_steps[4] + 1
end

function _som_change_weight(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    smin = P_SOM["smin"]
    γ = P_SOM["gamma"]

    t1 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t2 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    while t1 == t2
        t2 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    end
    w1 = 𝑆.tmp_conf[t1].w
    w2 = 𝑆.tmp_conf[t2].w
    h1 = 𝑆.tmp_conf[t1].h
    h2 = 𝑆.tmp_conf[t2].h
    dx_min = smin / w1 - h1
    dx_max = (h2 - smin / w2) * w2 / w1
    if dx_max ≤ dx_min
        return
    end
    dh = Pdx(dx_min, dx_max, γ, 𝑆.rng)

    _conf_size = length(𝑆.tmp_conf)
    𝑆.new_conf = deepcopy(𝑆.tmp_conf)
    𝑆.new_elem_dev = deepcopy(𝑆.elem_dev)
    𝑆.new_conf[t1].h = 𝑆.new_conf[t1].h + dh
    𝑆.new_conf[t2].h = 𝑆.new_conf[t2].h - dh * w1 / w2
    calc_dev_rec(𝑆.new_conf[t1], t1, 𝑆.new_elem_dev, ω)
    calc_dev_rec(𝑆.new_conf[t2], t2, 𝑆.new_elem_dev, ω)
    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = deepcopy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = deepcopy(𝑆.new_elem_dev)
        𝑆.accepted_steps[5] = 𝑆.accepted_steps[5] + 1
    end
    𝑆.trial_steps[5] = 𝑆.trial_steps[5] + 1
end

function _som_split(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    wmin = P_SOM["wmin"]
    smin = P_SOM["smin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(𝑆.rng, 1:length(𝑆.tmp_conf))

    old_conf = 𝑆.tmp_conf[t]
    if old_conf.w ≤ 2 * wmin || old_conf.w * old_conf.h ≤ 2.0 * smin
        return
    end

    h = old_conf.h
    w1 = wmin + (old_conf.w - 2.0 * wmin) * rand(𝑆.rng, F64)
    w2 = old_conf.w - w1
    if w1 > w2
        w1, w2 = w2, w1
    end

    c1 = old_conf.c - old_conf.w / 2.0 + w1 / 2.0
    c2 = old_conf.c + old_conf.w / 2.0 - w2 / 2.0
    dx_min = ommin + w1 / 2.0 - c1
    dx_max = ommax - w1 / 2.0 - c1
    if dx_max ≤ dx_min
        return
    end
    dc1 = Pdx(dx_min, dx_max, γ, 𝑆.rng)

    _conf_size = length(𝑆.tmp_conf)
    𝑆.new_conf = deepcopy(𝑆.tmp_conf)
    𝑆.new_elem_dev = deepcopy(𝑆.elem_dev)
    dc2 = -1.0 * w1 * dc1 / w2

    if (c1 + dc1 ≥ ommin + w1 / 2.0) &&
       (c1 + dc1 ≤ ommax - w1 / 2.0) &&
       (c2 + dc2 ≥ ommin + w2 / 2.0) &&
       (c2 + dc2 ≤ ommax - w2 / 2.0)

        𝑆.new_conf[t] = deepcopy(𝑆.new_conf[end])
        pop!(𝑆.new_conf)
        push!(𝑆.new_conf, Rectangle(h, w1, c1 + dc1))
        push!(𝑆.new_conf, Rectangle(h, w2, c2 + dc2))

        if t < _conf_size
            calc_dev_rec(𝑆.new_conf[t], t, 𝑆.new_elem_dev, ω)
        end
        calc_dev_rec(𝑆.new_conf[_conf_size], _conf_size, 𝑆.new_elem_dev, ω)
        calc_dev_rec(𝑆.new_conf[_conf_size+1], _conf_size+1, 𝑆.new_elem_dev, ω)
        𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)
        if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
            𝑆.tmp_conf = deepcopy(𝑆.new_conf)
            𝑆.tmp_dev = 𝑆.new_dev
            𝑆.elem_dev = deepcopy(𝑆.new_elem_dev)
            𝑆.accepted_steps[6] = 𝑆.accepted_steps[6] + 1
        end
    end
    𝑆.trial_steps[6] = 𝑆.trial_steps[6] + 1
end

function _som_merge(𝑆::T_SOM, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t1 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    t2 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    while t1 == t2
        t2 = rand(𝑆.rng, 1:length(𝑆.tmp_conf))
    end
    @assert t1 != t2

    old_conf1 = 𝑆.tmp_conf[t1]
    old_conf2 = 𝑆.tmp_conf[t2]

    weight = old_conf1.h * old_conf1.w + old_conf2.h * old_conf2.w
    w_new = 0.5 * (old_conf1.w + old_conf2.w)
    h_new = weight / w_new
    c_new = old_conf1.c + (old_conf2.c - old_conf1.c) * old_conf2.h * old_conf2.w / weight
    dx_min = ommin + w_new / 2.0 - c_new
    dx_max = ommax - w_new / 2.0 - c_new
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, γ, 𝑆.rng)

    _conf_size = length(𝑆.tmp_conf)
    𝑆.new_conf = deepcopy(𝑆.tmp_conf)
    𝑆.new_elem_dev = deepcopy(𝑆.elem_dev)

    if t1 > t2
        t1, t2 = t2, t1
    end

    𝑆.new_conf[t1] = deepcopy(Rectangle(h_new, w_new, c_new + dc))
    if t2 < _conf_size
        𝑆.new_conf[t2] = deepcopy(𝑆.new_conf[end])
    else
        @assert t2 == _conf_size
    end
    pop!(𝑆.new_conf)

    calc_dev_rec(𝑆.new_conf[t1], t1, 𝑆.new_elem_dev, ω)
    if t2 < _conf_size
        calc_dev_rec(𝑆.new_conf[t2], t2, 𝑆.new_elem_dev, ω)
    end

    calc_dev_rec(𝑆.new_conf[_conf_size - 1], _conf_size - 1, 𝑆.new_elem_dev, ω)
    𝑆.new_dev = calc_dev(𝑆.new_elem_dev, length(𝑆.new_conf), 𝐺)

    if rand(𝑆.rng, F64) < ((𝑆.tmp_dev / 𝑆.new_dev) ^ (1.0 + 𝑆.dacc))
        𝑆.tmp_conf = deepcopy(𝑆.new_conf)
        𝑆.tmp_dev = 𝑆.new_dev
        𝑆.elem_dev = deepcopy(𝑆.new_elem_dev)
        𝑆.accepted_steps[7] = 𝑆.accepted_steps[7] + 1
        #@show "test"
    end
    𝑆.trial_steps[7] = 𝑆.trial_steps[7] + 1
end

function calc_dev_rec(r::Rectangle, k::I64, elem_dev::Array{C64,2}, ω::FermionicMatsubaraGrid)
    Ngrid = P_SOM["Ngrid"]
    for g = 1:Ngrid
        Gs = r.h * log((im * ω.grid[g] - r.c + 0.5 * r.w) / (im * ω.grid[g] - r.c - 0.5 * r.w))
        elem_dev[g,k] = Gs
    end
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

function calc_norm(𝑆::T_SOM)
    norm1 = 0.0
    for i = 1:length(𝑆.tmp_conf)
        norm1 = norm1 + 𝑆.tmp_conf[i].h * 𝑆.tmp_conf[i].w
    end
    println("tmp norm is: $norm1")

    norm2 = 0.0
    for i = 1:length(𝑆.new_conf)
        norm2 = norm2 + 𝑆.new_conf[i].h * 𝑆.new_conf[i].w
    end
    println("new norm is: $norm2")

    if abs(norm1 - norm2) > 0.0001
        error()
    end
end

function Pdx(xmin::F64, xmax::F64, γ::F64, rng::AbstractRNG)
    _X = max(abs(xmin), abs(xmax))
    _λ = γ / _X
    _elx = exp(-1.0 * _λ * abs(xmin))
    _N = _λ / ( (xmin / abs(xmin)) * (exp(-1.0 * _λ * abs(xmin)) - 1.0)
              + (xmax / abs(xmax)) * (1.0 - exp(-1.0 * _λ * abs(xmax))) )

    y = rand(rng, F64)
    _lysn = _λ * y / _N
    if xmin ≥ 0
        return -1.0 * log(_elx - _lysn) / _λ
    elseif xmax ≤ 0
        return log(_lysn + _elx) / _λ
    else
        _C1 = _N * (1.0 - _elx) / _λ
        if y ≤ _C1
            return log(_lysn + _elx) / _λ
        else
            return -1.0 * log(1.0 - _lysn + _λ * _C1 / _N) / _λ
        end
    end
end
