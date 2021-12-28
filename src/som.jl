#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/28
#

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
    G :: Vector{C64}
    Δ :: F64
end

mutable struct SOMContext
    Cv :: Vector{Vector{Rectangle}}
    Δv :: Vector{F64}
end

function som_init()
    Lmax = P_SOM["Lmax"]
    Kmax = P_SOM["Kmax"]

    Δv = zeros(F64, Lmax)

    Cv = []
    for _ = 1:Lmax
        C = Rectangle[]
        for _ = 1:Kmax
            push!(C, Rectangle(0.0, 0.0, 0.0))
        end
        push!(Cv, C)
    end

    seed = rand(1:1000000)#;  seed = 840959
    rng = MersenneTwister(seed)
    @show "seed: ", seed
    tri = zeros(I64, 7)
    acc = zeros(I64, 7)

    return SOMContext(Cv, Δv), SOMMonteCarlo(rng, tri, acc)
end

function som_run(SC::SOMContext, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Lmax = P_SOM["Lmax"]
    for l = 1:Lmax
        println("try: $l")
        som_try(l, SC, MC, ω, 𝐺)
        som_output(l, SC)
    end
end

function som_try(l::I64, SC::SOMContext, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Nf = P_SOM["Nf"]

    SE = som_random(MC, ω, 𝐺)

    #@timev 
    for _ = 1:Nf
        som_update(SE, MC, ω, 𝐺)

        G = calc_gf(SE.Λ, length(SE.C))
        if sum( abs.(G - SE.G) ) / 64.0 > 0.00001
            error()
        end    
     end
     #error()

    SC.Δv[l] = SE.Δ
    SC.Cv[l] = deepcopy(SE.C)
end

function som_random(MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    smin = P_SOM["smin"]
    wmin = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    Kmax = P_SOM["Kmax"]
    Ngrid = P_SOM["Ngrid"]

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

    C = Rectangle[]
    Λ = zeros(C64, Ngrid, Kmax)
    Δ = 0.0

    for k = 1:_Know
        c = ommin + wmin / 2.0 + (ommax - ommin - wmin) * rand(MC.rng, F64)
        w = wmin + (min(2.0 * (c - ommin), 2.0 * (ommax - c)) - wmin) * rand(MC.rng, F64)
        h = weight[k] / w
        push!(C, Rectangle(h, w, c))
        calc_dev_rec(Rectangle(h, w, c), k, Λ, ω)
    end
    Δ = calc_dev(Λ, _Know, 𝐺)
    G = calc_gf(Λ, _Know)
    
    return SOMElement(C, Λ, G, Δ)
end

function som_update(SE::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Tmax = P_SOM["Tmax"]
    Kmax = P_SOM["Kmax"]
    dmax = P_SOM["dmax"]

    T1 = rand(MC.rng, 1:Tmax)
    d1 = rand(MC.rng, F64)
    d2 = 1.0 + (dmax - 1.0) * rand(MC.rng, F64)

    ST = deepcopy(SE)

    for _ = 1:T1
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

    for _ = T1+1:Tmax
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

    if ST.Δ < SE.Δ
        SE = deepcopy(ST)
        #SE.C = deepcopy(ST.C)
        #SE.Λ = deepcopy(ST.Λ)
        #SE.G = deepcopy(ST.G)
        #SE.Δ = ST.Δ
    end
end

function som_output(count::I64, 𝑆::SOMContext)
    alpha = P_SOM["alpha"]
    Ngrid = P_SOM["Ngrid"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]

    dev_min = minimum(𝑆.Δv[1:count])

    Lgood = 0
    Aom = zeros(F64, Ngrid)
    for l = 1:count
        if alpha * dev_min - 𝑆.Δv[l] > 0
            Lgood = Lgood + 1
            for w = 1:Ngrid
                _omega = ommin + (w - 1) * (ommax - ommin) / (Ngrid - 1)
                for r = 1:length(𝑆.Cv[l])
                    R = 𝑆.Cv[l][r]
                    if R.c - 0.5 * R.w ≤ _omega ≤ R.c + 0.5 * R.w
                        Aom[w] = Aom[w] + R.h
                    end
                end
            end
        end
    end

    @show count, 𝑆.Δv[1:count], dev_min, Lgood

    if Lgood > 0
        @. Aom = Aom / Lgood
    end

    open("Aw.data", "w") do fout
        for w = 1:Ngrid
            _omega = ommin + (w - 1) * (ommax - ommin) / (Ngrid - 1)
            println(fout, _omega, " ", Aom[w])
        end
    end
end

function _som_add(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    smin  = P_SOM["smin"]
    wmin  = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ     = P_SOM["gamma"]

    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]
    if R.h * R.w ≤ 2.0 * smin
        return
    end

    dx_min = smin
    dx_max = R.h * R.w - smin
    if dx_max ≤ dx_min
        return
    end

    r1 = rand(MC.rng, F64)
    r2 = rand(MC.rng, F64)
    c = (ommin + wmin / 2.0) + (ommax - ommin - wmin) * r1
    w_new_max = 2.0 * min(ommax - c, c - ommin)
    dx = Pdx(dx_min, dx_max, γ, MC.rng)
    h = dx / w_new_max + (dx / wmin - dx / w_new_max) * r2
    w = dx / h

    Rnew = Rectangle(R.h - dx / R.w, R.w, R.c)
    Radd = Rectangle(h, w, c)

    G1 = calc_dev_rec(R, ω)
    G2 = calc_dev_rec(Rnew, ω)
    G3 = calc_dev_rec(Radd, ω)
    new_dev = calc_dev(𝑆.G - G1 + G2 + G3, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
        𝑆.C[t] = Rnew
        push!(𝑆.C, Radd)
        𝑆.Δ = new_dev
        @. 𝑆.G = 𝑆.G - G1 + G2 + G3
        @. 𝑆.Λ[:,t] = G2
        @. 𝑆.Λ[:,csize+1] = G3
        MC.acc[1] = MC.acc[1] + 1
    end
    MC.tri[1] = MC.tri[1] + 1
end

function _som_remove(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    csize = length(𝑆.C)

    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end
    if t1 < t2
        t1, t2 = t2, t1
    end

    R1 = 𝑆.C[t1]
    R2 = 𝑆.C[t2]
    Re = 𝑆.C[end]

    dx = R1.h * R1.w

    G1 = calc_dev_rec(R1, ω)
    G2 = calc_dev_rec(R2, ω)
    Ge = calc_dev_rec(Re, ω)

    R2n = Rectangle(R2.h + dx / R2.w, R2.w, R2.c)
    G2n = calc_dev_rec(R2n, ω)

    new_dev = calc_dev(𝑆.G - G1 - G2 + G2n, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
        𝑆.C[t2] = R2n
        if t1 < csize
            𝑆.C[t1] = Re
        end
        pop!(𝑆.C)
        𝑆.Δ = new_dev
        @. 𝑆.G = 𝑆.G - G1 - G2 + G2n
        @. 𝑆.Λ[:,t2] = G2n
        if t1 < csize
            @. 𝑆.Λ[:,t1] = Ge
        end
        MC.acc[2] = MC.acc[2] + 1
    end
    MC.tri[2] = MC.tri[2] + 1
end

function _som_shift(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ     = P_SOM["gamma"]

    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]
    dx_min = ommin + R.w / 2.0 - R.c
    dx_max = ommax - R.w / 2.0 - R.c
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, γ, MC.rng)
    
    Rn = Rectangle(R.h, R.w, R.c + dc)
    G1 = calc_dev_rec(R, ω)
    G2 = calc_dev_rec(Rn, ω)
    new_dev = calc_dev(𝑆.G - G1 + G2, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ / new_dev) ^ (1.0 + dacc))
        𝑆.C[t] = Rn
        𝑆.Δ = new_dev
        @. 𝑆.G = 𝑆.G - G1 + G2
        @. 𝑆.Λ[:,t] = G2
        MC.acc[3] = MC.acc[3] + 1
    end
    MC.tri[3] = MC.tri[3] + 1
end

function _som_change_width(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    wmin  = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ     = P_SOM["gamma"]

    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]
    weight = R.h * R.w
    dx_min = wmin - R.w
    dx_max = min(2.0 * (R.c - ommin), 2.0 * (ommax - R.c)) - R.w
    if dx_max ≤ dx_min
        return
    end
    dw = Pdx(dx_min, dx_max, γ, MC.rng)

    w = R.w + dw
    h = weight / w
    c = R.c
    Rn = Rectangle(h, w, c)
    G1 = calc_dev_rec(R, ω)
    G2 = calc_dev_rec(Rn, ω)
    new_dev = calc_dev(𝑆.G - G1 + G2, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/ new_dev) ^ (1.0 + dacc))
        𝑆.C[t] = Rn
        𝑆.Δ = new_dev
        @. 𝑆.G = 𝑆.G - G1 + G2
        @. 𝑆.Λ[:,t] = G2
        MC.acc[4] = MC.acc[4] + 1
    end
    MC.tri[4] = MC.tri[4] + 1
end

function _som_change_weight(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    smin = P_SOM["smin"]
    γ    = P_SOM["gamma"]

    csize = length(𝑆.C)

    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end

    R1 = 𝑆.C[t1]
    R2 = 𝑆.C[t2]

    w1 = R1.w
    w2 = R2.w
    h1 = R1.h
    h2 = R2.h
    dx_min = smin / w1 - h1
    dx_max = (h2 - smin / w2) * w2 / w1
    if dx_max ≤ dx_min
        return
    end
    dh = Pdx(dx_min, dx_max, γ, MC.rng)

    R1n = Rectangle(R1.h + dh, R1.w, R1.c)
    G1A = calc_dev_rec(R1, ω)
    G1B = calc_dev_rec(R1n, ω)
    
    R2n = Rectangle(R2.h - dh * w1 / w2, R2.w, R2.c)
    G2A = calc_dev_rec(R2, ω)
    G2B = calc_dev_rec(R2n, ω)
    new_dev = calc_dev(𝑆.G - G1A + G1B - G2A + G2B, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
        𝑆.C[t1] = R1n
        𝑆.C[t2] = R2n
        𝑆.Δ = new_dev
        @. 𝑆.G = 𝑆.G - G1A + G1B - G2A + G2B
        @. 𝑆.Λ[:,t1] = G1B
        @. 𝑆.Λ[:,t2] = G2B
        MC.acc[5] = MC.acc[5] + 1
    end
    MC.tri[5] = MC.tri[5] + 1
end

function _som_split(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    wmin  = P_SOM["wmin"]
    smin  = P_SOM["smin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ     = P_SOM["gamma"]

    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R1 = 𝑆.C[t]
    Re = 𝑆.C[end]

    if R1.w ≤ 2 * wmin || R1.w * R1.h ≤ 2.0 * smin
        return
    end

    h = R1.h
    w1 = wmin + (R1.w - 2.0 * wmin) * rand(MC.rng, F64)
    w2 = R1.w - w1
    if w1 > w2
        w1, w2 = w2, w1
    end
    c1 = R1.c - R1.w / 2.0 + w1 / 2.0
    c2 = R1.c + R1.w / 2.0 - w2 / 2.0
    dx_min = ommin + w1 / 2.0 - c1
    dx_max = ommax - w1 / 2.0 - c1
    if dx_max ≤ dx_min
        return
    end
    dc1 = Pdx(dx_min, dx_max, γ, MC.rng)
    dc2 = -1.0 * w1 * dc1 / w2

    if (c1 + dc1 ≥ ommin + w1 / 2.0) &&
       (c1 + dc1 ≤ ommax - w1 / 2.0) &&
       (c2 + dc2 ≥ ommin + w2 / 2.0) &&
       (c2 + dc2 ≤ ommax - w2 / 2.0)

        G1 = calc_dev_rec(R1, ω)
        Ge = calc_dev_rec(Re, ω)

        R2 = Rectangle(h, w1, c1 + dc1)
        G2 = calc_dev_rec(R2, ω)

        R3 = Rectangle(h, w2, c2 + dc2)
        G3 = calc_dev_rec(R3, ω)
        new_dev = calc_dev(𝑆.G - G1 + G2 + G3, 𝐺)

        if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
            𝑆.C[t] = deepcopy(𝑆.C[end])
            pop!(𝑆.C)
            push!(𝑆.C, R2)
            push!(𝑆.C, R3)
            𝑆.Δ = new_dev
            @. 𝑆.G = 𝑆.G - G1 + G2 + G3
            
            if t < csize
                @. 𝑆.Λ[:,t] = Ge
            end
            @. 𝑆.Λ[:,csize] = G2
            @. 𝑆.Λ[:,csize+1] = G3
            MC.acc[6] = MC.acc[6] + 1
        end
    end
    MC.tri[6] = MC.tri[6] + 1
end

function _som_merge(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    γ     = P_SOM["gamma"]

    csize = length(𝑆.C)

    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end
    if t1 > t2
        t1, t2 = t2, t1
    end

    R1 = 𝑆.C[t1]
    R2 = 𝑆.C[t2]

    weight = R1.h * R1.w + R2.h * R2.w
    w_new = 0.5 * (R1.w + R2.w)
    h_new = weight / w_new
    c_new = R1.c + (R2.c - R1.c) * R2.h * R2.w / weight
    dx_min = ommin + w_new / 2.0 - c_new
    dx_max = ommax - w_new / 2.0 - c_new
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, γ, MC.rng)

    G1 = calc_dev_rec(R1, ω)
    G2 = calc_dev_rec(R2, ω)

    Rn = Rectangle(h_new, w_new, c_new + dc)
    Re = 𝑆.C[end]
    Gn = calc_dev_rec(Rn, ω)
    Ge = calc_dev_rec(Re, ω)

    new_dev = calc_dev(𝑆.G - G1 - G2 + Gn, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/new_dev) ^ (1.0 + dacc))
        𝑆.C[t1] = Rn
        if t2 < csize
            𝑆.C[t2] = deepcopy(𝑆.C[end])
        end
        pop!(𝑆.C)
        𝑆.Δ = new_dev
        @. 𝑆.G = 𝑆.G - G1 - G2 + Gn

        @. 𝑆.Λ[:,t1] = Gn
        if t2 < csize
            @. 𝑆.Λ[:,t2] = Ge
        end

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

function calc_dev_rec(r::Rectangle, ω::FermionicMatsubaraGrid)
    Ngrid = P_SOM["Ngrid"]
    elem_dev = zeros(C64, Ngrid)
    for g = 1:Ngrid
        elem_dev[g] = r.h * log((im * ω.grid[g] - r.c + 0.5 * r.w) / (im * ω.grid[g] - r.c - 0.5 * r.w))
    end
    return elem_dev
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

function calc_dev(Gc::Vector{C64}, 𝐺::GreenData)
    return sum( abs.((Gc .- 𝐺.value) ./ (𝐺.error)) )
end

function calc_gf(elem_dev::Array{C64,2}, nk::I64)
    Ngrid = P_SOM["Ngrid"]

    G = zeros(C64, Ngrid)
    for g = 1:Ngrid
        for k = 1:nk
            G[g] = G[g] + elem_dev[g,k]
        end
    end
    return G
end

function calc_norm(V1::Vector{Rectangle}, V2::Vector{Rectangle})
    norm1 = 0.0
    for i = 1:length(V1)
        norm1 = norm1 + V1[i].h * V1[i].w
    end
    println("tmp norm is: $norm1")

    norm2 = 0.0
    for i = 1:length(V2)
        norm2 = norm2 + V2[i].h * V2[i].w
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
