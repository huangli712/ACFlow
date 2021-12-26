#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/26
#

mutable struct Rectangle
    h :: F64
    w :: F64
    c :: F64
end

abstract type AbstractMonteCarlo end
mutable struct SOMMonteCarlo <: AbstractMonteCarlo
    rng :: AbstractRNG
    tri :: Vector{I64}
    acc :: Vector{I64}
end

mutable struct SOMElement
    C :: Vector{Rectangle}
    Λ :: Array{C64,2}
    Δ :: F64
end

mutable struct SOMContext
    dev :: Vector{F64}
    conf :: Vector{Vector{Rectangle}}

    att_conf :: Vector{Rectangle}
    tmp_conf :: Vector{Rectangle}

    att_elem_dev :: Array{C64,2}
    tmp_elem_dev :: Array{C64,2}

    att_dev :: F64
    tmp_dev :: F64
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
    for k = 1:Kmax
        push!(att_conf, Rectangle(0.0, 0.0, 0.0))
        push!(tmp_conf, Rectangle(0.0, 0.0, 0.0))
    end

    att_elem_dev = zeros(C64, Ngrid, Kmax)
    tmp_elem_dev = zeros(C64, Ngrid, Kmax)

    trial_steps = zeros(I64, 7)
    accepted_steps = zeros(I64, 7)

    return SOMContext(dev,
                 conf,
                 att_conf,
                 tmp_conf,
                 att_elem_dev,
                 tmp_elem_dev,
                 0.0, 0.0), SOMMonteCarlo(rng, trial_steps, accepted_steps)
end

function som_run(𝑆::SOMContext, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Lmax = P_SOM["Lmax"]
    for l = 1:Lmax
        println("try: $l")
        som_try(l, 𝑆, MC, ω, 𝐺)
        som_output(l, 𝑆)
    end
end

function som_try(l, 𝑆::SOMContext, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Nf = P_SOM["Nf"]
    som_random(𝑆, MC, ω, 𝐺)

    for f = 1:Nf
        som_update(𝑆, MC, ω, 𝐺)
    end

    𝑆.dev[l] = 𝑆.att_dev
    𝑆.conf[l] = deepcopy(𝑆.att_conf)
end

function som_random(𝑆::SOMContext, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    smin = P_SOM["smin"]
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    Kmax = P_SOM["Kmax"]

    _Know = rand(MC.rng, 2:Kmax)
    _weight = zeros(F64, _Know)
    for i = 1:_Know
        _weight[i] = rand(MC.rng, F64)
    end
    _weight[end] = 1.0

    sort!(_weight)
    weight = diff(_weight)
    insert!(weight, 1, _weight[1])
    sort!(weight)

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

    empty!(𝑆.att_conf)
    fill!(𝑆.att_elem_dev, zero(C64))

    for k = 1:_Know
        c = ommin + wmin / 2.0 + (ommax - ommin - wmin) * rand(MC.rng, F64)
        w = wmin + (min(2.0 * (c - ommin), 2.0 * (ommax - c)) - wmin) * rand(MC.rng, F64)
        h = weight[k] / w
        push!(𝑆.att_conf, Rectangle(h, w, c))
        calc_dev_rec(Rectangle(h, w, c), k, 𝑆.att_elem_dev, ω)
    end
    𝑆.att_dev = calc_dev(𝑆.att_elem_dev, _Know, 𝐺)
end

function som_update(SA::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Tmax = P_SOM["Tmax"]
    Kmax = P_SOM["Kmax"]
    dmax = P_SOM["dmax"]
    T1 = rand(MC.rng, 1:Tmax)

    d1 = rand(MC.rng, F64)
    d2 = 1.0 + (dmax - 1.0) * rand(MC.rng, F64)

    ST = deepcopy(SA)

    for i = 1:T1
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < Kmax - 1
                    _som_add(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 2
                if length(ST.C) > 1
                    _som_remove(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 3
                _som_shift(ST, MC, ω, 𝐺, d1)
                break

            @case 4
                _som_change_width(ST, MC, ω, 𝐺, d1)
                break

            @case 5
                if length(ST.C) > 1
                    _som_change_weight(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 6
                if length(ST.C) < Kmax - 1
                    _som_split(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 7
                if length(ST.C) > 1
                    _som_merge(ST, MC, ω, 𝐺, d1)
                end
                break
        end
    end

    for j = T1+1:Tmax
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < Kmax - 1
                    _som_add(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 2
                if length(ST.C) > 1
                    _som_remove(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 3
                _som_shift(ST, MC, ω, 𝐺, d2)
                break

            @case 4
                _som_change_width(ST, MC, ω, 𝐺, d2)
                break

            @case 5
                if length(ST.C) > 1
                    _som_change_weight(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 6
                if length(ST.C) < Kmax - 1
                    _som_split(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 7
                if length(ST.C) > 1
                    _som_merge(ST, MC, ω, 𝐺, d2)
                end
                break
        end
    end

    if ST.Δ < SA.Δ
        SA = deepcopy(ST)
    end
end

function som_output(count::I64, 𝑆::SOMContext)
    println("output")
    alpha = P_SOM["alpha"]
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

function _som_add(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    smin = P_SOM["smin"]
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(MC.rng, 1:length(𝑆.C))
    if 𝑆.C[t].h * 𝑆.C[t].w ≤ 2.0 * smin
        return
    end

    dx_min = smin
    dx_max = 𝑆.C[t].h * 𝑆.C[t].w - smin
    if dx_max ≤ dx_min
        return
    end

    c = (ommin + wmin / 2.0) + (ommax - ommin - wmin) * rand(MC.rng, F64)
    w_new_max = 2.0 * min(ommax - c, c - ommin)
    dx = Pdx(dx_min, dx_max, γ, MC.rng)

    r = rand(MC.rng, F64)
    new_conf = deepcopy(𝑆.C)
    new_elem_dev = deepcopy(𝑆.Λ)
    h = dx / w_new_max + (dx / wmin - dx / w_new_max) * r
    w = dx / h

    push!(new_conf, Rectangle(h, w, c))
    new_conf[t].h = new_conf[t].h - dx / new_conf[t].w
    calc_dev_rec(new_conf[t], t, new_elem_dev, ω)
    calc_dev_rec(new_conf[end], length(new_conf), new_elem_dev, ω)
    new_dev = calc_dev(new_elem_dev, length(new_conf), 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/ new_dev) ^ (1.0 + dacc))
        𝑆.C = deepcopy(new_conf)
        𝑆.Δ = new_dev
        𝑆.Λ = deepcopy(new_elem_dev)
        MC.acc[1] = MC.acc[1] + 1
    end
    MC.tri[1] = MC.tri[1] + 1
end

function _som_remove(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    t1 = rand(MC.rng, 1:length(𝑆.C))
    t2 = rand(MC.rng, 1:length(𝑆.C))
    while t1 == t2
        t2 = rand(MC.rng, 1:length(𝑆.C))
    end

    if t1 < t2
        t1, t2 = t2, t1
    end

    _conf_size = length(𝑆.C)
    dx = 𝑆.C[t1].h * 𝑆.C[t1].w

    new_conf = deepcopy(𝑆.C)
    new_elem_dev = deepcopy(𝑆.Λ)


    new_conf[t2].h = new_conf[t2].h + dx / new_conf[t2].w
    if t1 < _conf_size
        new_conf[t1] = deepcopy(new_conf[end])
    else
        @assert t1 == _conf_size
    end
    pop!(new_conf)

    if t1 < _conf_size
        calc_dev_rec(new_conf[t1], t1, new_elem_dev, ω)
    end
    calc_dev_rec(new_conf[t2], t2, new_elem_dev, ω)

    new_dev = calc_dev(new_elem_dev, length(new_conf), 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/ new_dev) ^ (1.0 + dacc))
        𝑆.C = deepcopy(new_conf)
        𝑆.Δ = new_dev
        𝑆.Λ = deepcopy(new_elem_dev)
        MC.acc[2] = MC.acc[2] + 1
    end
    MC.tri[2] = MC.tri[2] + 1
end

function _som_shift(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(MC.rng, 1:length(𝑆.C))

    dx_min = ommin + 𝑆.C[t].w / 2.0 - 𝑆.C[t].c
    dx_max = ommax - 𝑆.C[t].w / 2.0 - 𝑆.C[t].c
    if dx_max ≤ dx_min
        return
    end

    dc = Pdx(dx_min, dx_max, γ, MC.rng)

    _conf_size = length(𝑆.C)
    new_conf = deepcopy(𝑆.C)
    new_elem_dev = deepcopy(𝑆.Λ)
    new_conf[t].c = new_conf[t].c + dc

    calc_dev_rec(new_conf[t], t, new_elem_dev, ω)
    new_dev = calc_dev(new_elem_dev, length(new_conf), 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ / new_dev) ^ (1.0 + dacc))
        𝑆.C = deepcopy(new_conf)
        𝑆.Δ = new_dev
        𝑆.Λ = deepcopy(new_elem_dev)
        MC.acc[3] = MC.acc[3] + 1
    end
    MC.tri[3] = MC.tri[3] + 1
end

function _som_change_width(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(MC.rng, 1:length(𝑆.C))

    weight = 𝑆.C[t].h * 𝑆.C[t].w
    dx_min = wmin - 𝑆.C[t].w
    dx_max = min(2.0 * (𝑆.C[t].c - ommin), 2.0 * (ommax - 𝑆.C[t].c)) - 𝑆.C[t].w
    if dx_max ≤ dx_min
        return
    end
    dw = Pdx(dx_min, dx_max, γ, MC.rng)

    _conf_size = length(𝑆.C)
    new_conf = deepcopy(𝑆.C)
    new_elem_dev = deepcopy(𝑆.Λ)
    new_conf[t].w = new_conf[t].w + dw
    new_conf[t].h = weight / new_conf[t].w
    calc_dev_rec(new_conf[t], t, new_elem_dev, ω)

    new_dev = calc_dev(new_elem_dev, length(new_conf), 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/ new_dev) ^ (1.0 + dacc))
        𝑆.C = deepcopy(new_conf)
        𝑆.Δ = new_dev
        𝑆.Λ = deepcopy(new_elem_dev)
        MC.acc[4] = MC.acc[4] + 1
    end
    MC.tri[4] = MC.tri[4] + 1
end

function _som_change_weight(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    smin = P_SOM["smin"]
    γ = P_SOM["gamma"]

    t1 = rand(MC.rng, 1:length(𝑆.C))
    t2 = rand(MC.rng, 1:length(𝑆.C))
    while t1 == t2
        t2 = rand(MC.rng, 1:length(𝑆.C))
    end
    w1 = 𝑆.C[t1].w
    w2 = 𝑆.C[t2].w
    h1 = 𝑆.C[t1].h
    h2 = 𝑆.C[t2].h
    dx_min = smin / w1 - h1
    dx_max = (h2 - smin / w2) * w2 / w1
    if dx_max ≤ dx_min
        return
    end
    dh = Pdx(dx_min, dx_max, γ, MC.rng)

    _conf_size = length(𝑆.C)
    new_conf = deepcopy(𝑆.C)
    new_elem_dev = deepcopy(𝑆.Λ)
    new_conf[t1].h = new_conf[t1].h + dh
    new_conf[t2].h = new_conf[t2].h - dh * w1 / w2
    calc_dev_rec(new_conf[t1], t1, new_elem_dev, ω)
    calc_dev_rec(new_conf[t2], t2, new_elem_dev, ω)
    new_dev = calc_dev(new_elem_dev, length(new_conf), 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
        𝑆.C = deepcopy(new_conf)
        𝑆.Δ = new_dev
        𝑆.Λ = deepcopy(new_elem_dev)
        MC.acc[5] = MC.acc[5] + 1
    end
    MC.tri[5] = MC.tri[5] + 1
end

function _som_split(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    wmin = P_SOM["wmin"]
    smin = P_SOM["smin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t = rand(MC.rng, 1:length(𝑆.C))

    old_conf = 𝑆.C[t]
    if old_conf.w ≤ 2 * wmin || old_conf.w * old_conf.h ≤ 2.0 * smin
        return
    end

    h = old_conf.h
    w1 = wmin + (old_conf.w - 2.0 * wmin) * rand(MC.rng, F64)
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
    dc1 = Pdx(dx_min, dx_max, γ, MC.rng)

    _conf_size = length(𝑆.C)
    new_conf = deepcopy(𝑆.C)
    new_elem_dev = deepcopy(𝑆.Λ)
    dc2 = -1.0 * w1 * dc1 / w2

    if (c1 + dc1 ≥ ommin + w1 / 2.0) &&
       (c1 + dc1 ≤ ommax - w1 / 2.0) &&
       (c2 + dc2 ≥ ommin + w2 / 2.0) &&
       (c2 + dc2 ≤ ommax - w2 / 2.0)

        new_conf[t] = deepcopy(new_conf[end])
        pop!(new_conf)
        push!(new_conf, Rectangle(h, w1, c1 + dc1))
        push!(new_conf, Rectangle(h, w2, c2 + dc2))

        if t < _conf_size
            calc_dev_rec(new_conf[t], t, new_elem_dev, ω)
        end
        calc_dev_rec(new_conf[_conf_size], _conf_size, new_elem_dev, ω)
        calc_dev_rec(new_conf[_conf_size+1], _conf_size+1, new_elem_dev, ω)
        new_dev = calc_dev(new_elem_dev, length(new_conf), 𝐺)
        if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
            𝑆.C = deepcopy(new_conf)
            𝑆.Δ = new_dev
            𝑆.Λ = deepcopy(new_elem_dev)
            MC.acc[6] = MC.acc[6] + 1
        end
    end
    MC.tri[6] = MC.tri[6] + 1
end

function _som_merge(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ = P_SOM["gamma"]

    t1 = rand(MC.rng, 1:length(𝑆.C))
    t2 = rand(MC.rng, 1:length(𝑆.C))
    while t1 == t2
        t2 = rand(MC.rng, 1:length(𝑆.C))
    end
    @assert t1 != t2

    old_conf1 = 𝑆.C[t1]
    old_conf2 = 𝑆.C[t2]

    weight = old_conf1.h * old_conf1.w + old_conf2.h * old_conf2.w
    w_new = 0.5 * (old_conf1.w + old_conf2.w)
    h_new = weight / w_new
    c_new = old_conf1.c + (old_conf2.c - old_conf1.c) * old_conf2.h * old_conf2.w / weight
    dx_min = ommin + w_new / 2.0 - c_new
    dx_max = ommax - w_new / 2.0 - c_new
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, γ, MC.rng)

    _conf_size = length(𝑆.C)
    new_conf = deepcopy(𝑆.C)
    new_elem_dev = deepcopy(𝑆.Λ)

    if t1 > t2
        t1, t2 = t2, t1
    end

    new_conf[t1] = deepcopy(Rectangle(h_new, w_new, c_new + dc))
    if t2 < _conf_size
        new_conf[t2] = deepcopy(new_conf[end])
    else
        @assert t2 == _conf_size
    end
    pop!(new_conf)

    calc_dev_rec(new_conf[t1], t1, new_elem_dev, ω)
    if t2 < _conf_size
        calc_dev_rec(new_conf[t2], t2, new_elem_dev, ω)
    end

    calc_dev_rec(new_conf[_conf_size - 1], _conf_size - 1, new_elem_dev, ω)
    new_dev = calc_dev(new_elem_dev, length(new_conf), 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
        𝑆.C = deepcopy(new_conf)
        𝑆.Δ = new_dev
        𝑆.Λ = deepcopy(new_elem_dev)
        MC.acc[7] = MC.acc[7] + 1
    end
    MC.tri[7] = MC.tri[7] + 1
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

function calc_norm(𝑆::SOMContext)
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
