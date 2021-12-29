#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/29
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

function som_run(ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    Lmax = P_SOM["Lmax"]
    Nf = P_SOM["Nf"]

    SC, MC = som_init()

    for l = 1:Lmax
        println("try: $l")

        SE = som_random(MC, ω, 𝐺)
    
        for _ = 1:Nf
            som_update(SE, MC, ω, 𝐺)
        end
    
        SC.Δv[l] = SE.Δ
        SC.Cv[l] = deepcopy(SE.C)    
    end

    return som_spectra(SC)
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

    seed = rand(1:1000000)#;  seed = 112414
    rng = MersenneTwister(seed)
    @show "seed: ", seed
    tri = zeros(I64, 7)
    acc = zeros(I64, 7)

    return SOMContext(Cv, Δv), SOMMonteCarlo(rng, tri, acc)
end

function som_random(MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData)
    smin  = P_SOM["smin"]
    wmin  = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    Kmax  = P_SOM["Kmax"]
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
        R = Rectangle(h, w, c)
        push!(C, R)
        Λ[:,k] .= _calc_lambda(R, ω)
    end
    Δ = _calc_err(Λ, _Know, 𝐺)
    G = _calc_gf(Λ, _Know)
    
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
                    _try_insert(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 2
                if length(ST.C) > 1
                    _try_remove(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 3
                _try_position(ST, MC, ω, 𝐺, d1)
                break

            @case 4
                _try_width(ST, MC, ω, 𝐺, d1)
                break

            @case 5
                if length(ST.C) > 1
                    _try_height(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 6
                if length(ST.C) < Kmax - 1
                    _try_split(ST, MC, ω, 𝐺, d1)
                end
                break

            @case 7
                if length(ST.C) > 1
                    _try_merge(ST, MC, ω, 𝐺, d1)
                end
                break
        end

    end

    for _ = T1+1:Tmax
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < Kmax - 1
                    _try_insert(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 2
                if length(ST.C) > 1
                    _try_remove(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 3
                _try_position(ST, MC, ω, 𝐺, d2)
                break

            @case 4
                _try_width(ST, MC, ω, 𝐺, d2)
                break

            @case 5
                if length(ST.C) > 1
                    _try_height(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 6
                if length(ST.C) < Kmax - 1
                    _try_split(ST, MC, ω, 𝐺, d2)
                end
                break

            @case 7
                if length(ST.C) > 1
                    _try_merge(ST, MC, ω, 𝐺, d2)
                end
                break
        end
    end

    if ST.Δ < SE.Δ
        SE.C = deepcopy(ST.C)
        SE.Λ .= ST.Λ
        SE.G .= ST.G
        SE.Δ  = ST.Δ
    end
end

function som_spectra(𝑆::SOMContext)
    alpha = P_SOM["alpha"]
    Ngrid = P_SOM["Ngrid"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    Lmax  = P_SOM["Lmax"]

    dev_min = minimum(𝑆.Δv)

    Lgood = 0
    Aom = zeros(F64, Ngrid)
    for l = 1:Lmax
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

    @show 𝑆.Δv, dev_min, Lgood

    if Lgood > 0
        @. Aom = Aom / Lgood
    end

    return Aom
end

function som_output(Aom::Vector{F64})
    Ngrid = P_SOM["Ngrid"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]

    open("Aw.data", "w") do fout
        for w = 1:Ngrid
            _omega = ommin + (w - 1) * (ommax - ommin) / (Ngrid - 1)
            println(fout, _omega, " ", Aom[w])
        end
    end
end

function _try_insert(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    smin  = P_SOM["smin"]
    wmin  = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
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
    dx = Pdx(dx_min, dx_max, MC.rng)
    h = dx / w_new_max + (dx / wmin - dx / w_new_max) * r2
    w = dx / h

    Rnew = Rectangle(R.h - dx / R.w, R.w, R.c)
    Radd = Rectangle(h, w, c)

    G1 = 𝑆.Λ[:,t]
    G2 = _calc_lambda(Rnew, ω)
    G3 = _calc_lambda(Radd, ω)

    Δ = _calc_err(𝑆.G - G1 + G2 + G3, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t] = Rnew
        push!(𝑆.C, Radd)
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1 + G2 + G3
        @. 𝑆.Λ[:,t] = G2
        @. 𝑆.Λ[:,csize+1] = G3
        MC.acc[1] = MC.acc[1] + 1
    end

    MC.tri[1] = MC.tri[1] + 1
end

function _try_remove(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
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

    G1 = 𝑆.Λ[:,t1]
    G2 = 𝑆.Λ[:,t2]
    Ge = 𝑆.Λ[:,csize]

    R2n = Rectangle(R2.h + dx / R2.w, R2.w, R2.c)
    G2n = _calc_lambda(R2n, ω)

    Δ = _calc_err(𝑆.G - G1 - G2 + G2n, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t2] = R2n
        if t1 < csize
            𝑆.C[t1] = Re
        end
        pop!(𝑆.C)
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1 - G2 + G2n
        @. 𝑆.Λ[:,t2] = G2n
        if t1 < csize
            @. 𝑆.Λ[:,t1] = Ge
        end
        MC.acc[2] = MC.acc[2] + 1
    end

    MC.tri[2] = MC.tri[2] + 1
end

function _try_position(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]

    dx_min = ommin + R.w / 2.0 - R.c
    dx_max = ommax - R.w / 2.0 - R.c
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, MC.rng)
    
    Rn = Rectangle(R.h, R.w, R.c + dc)
    G1 = 𝑆.Λ[:,t]
    G2 = _calc_lambda(Rn, ω)

    Δ = _calc_err(𝑆.G - G1 + G2, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t] = Rn
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1 + G2
        @. 𝑆.Λ[:,t] = G2
        MC.acc[3] = MC.acc[3] + 1
    end

    MC.tri[3] = MC.tri[3] + 1
end

function _try_width(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    wmin  = P_SOM["wmin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]

    weight = R.h * R.w
    dx_min = wmin - R.w
    dx_max = min(2.0 * (R.c - ommin), 2.0 * (ommax - R.c)) - R.w
    if dx_max ≤ dx_min
        return
    end
    dw = Pdx(dx_min, dx_max, MC.rng)
    w = R.w + dw
    h = weight / w
    c = R.c

    Rn = Rectangle(h, w, c)
    G1 = 𝑆.Λ[:,t]
    G2 = _calc_lambda(Rn, ω)

    Δ = _calc_err(𝑆.G - G1 + G2, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t] = Rn
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1 + G2
        @. 𝑆.Λ[:,t] = G2
        MC.acc[4] = MC.acc[4] + 1
    end

    MC.tri[4] = MC.tri[4] + 1
end

function _try_height(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    smin  = P_SOM["smin"]
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
    dh = Pdx(dx_min, dx_max, MC.rng)

    R1n = Rectangle(R1.h + dh, R1.w, R1.c)
    G1A = 𝑆.Λ[:,t1]
    G1B = _calc_lambda(R1n, ω)
    R2n = Rectangle(R2.h - dh * w1 / w2, R2.w, R2.c)
    G2A = 𝑆.Λ[:,t2]
    G2B = _calc_lambda(R2n, ω)

    Δ = _calc_err(𝑆.G - G1A + G1B - G2A + G2B, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t1] = R1n
        𝑆.C[t2] = R2n
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1A + G1B - G2A + G2B
        @. 𝑆.Λ[:,t1] = G1B
        @. 𝑆.Λ[:,t2] = G2B
        MC.acc[5] = MC.acc[5] + 1
    end

    MC.tri[5] = MC.tri[5] + 1
end

function _try_split(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    wmin  = P_SOM["wmin"]
    smin  = P_SOM["smin"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R1 = 𝑆.C[t]
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
    dc1 = Pdx(dx_min, dx_max, MC.rng)
    dc2 = -1.0 * w1 * dc1 / w2

    if (c1 + dc1 ≥ ommin + w1 / 2.0) &&
       (c1 + dc1 ≤ ommax - w1 / 2.0) &&
       (c2 + dc2 ≥ ommin + w2 / 2.0) &&
       (c2 + dc2 ≤ ommax - w2 / 2.0)

        G1 = 𝑆.Λ[:,t]
        Ge = 𝑆.Λ[:,csize]

        R2 = Rectangle(h, w1, c1 + dc1)
        G2 = _calc_lambda(R2, ω)

        R3 = Rectangle(h, w2, c2 + dc2)
        G3 = _calc_lambda(R3, ω)
        Δ = _calc_err(𝑆.G - G1 + G2 + G3, 𝐺)

        if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
            𝑆.C[t] = 𝑆.C[end]
            pop!(𝑆.C)
            push!(𝑆.C, R2)
            push!(𝑆.C, R3)
            𝑆.Δ = Δ
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

function _try_merge(𝑆::SOMElement, MC::SOMMonteCarlo, ω::FermionicMatsubaraGrid, 𝐺::GreenData, dacc)
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]
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
    dc = Pdx(dx_min, dx_max, MC.rng)

    G1 = 𝑆.Λ[:,t1]
    G2 = 𝑆.Λ[:,t2]
    Ge = 𝑆.Λ[:,csize]

    Rn = Rectangle(h_new, w_new, c_new + dc)
    Gn = _calc_lambda(Rn, ω)

    Δ = _calc_err(𝑆.G - G1 - G2 + Gn, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t1] = Rn
        if t2 < csize
            𝑆.C[t2] = 𝑆.C[end]
        end
        pop!(𝑆.C)
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1 - G2 + Gn
        @. 𝑆.Λ[:,t1] = Gn
        if t2 < csize
            @. 𝑆.Λ[:,t2] = Ge
        end
        MC.acc[7] = MC.acc[7] + 1
    end

    MC.tri[7] = MC.tri[7] + 1
end

function _calc_lambda(r::Rectangle, ω::FermionicMatsubaraGrid)
    Λ = @. r.h * log((im * ω.grid - r.c + 0.5 * r.w) / (im * ω.grid - r.c - 0.5 * r.w))
    return Λ
end

function _calc_err(Λ::Array{C64,2}, nk::I64, 𝐺::GreenData)
    Ngrid, Kmax = size(Λ)
    @assert nk ≤ Kmax

    res = 0.0
    for w = 1:Ngrid
        g = sum(Λ[w,1:nk])
        res = res + abs((g - 𝐺.value[w]) / 𝐺.error[w])
    end

    return res
end

function _calc_err(Gc::Vector{C64}, 𝐺::GreenData)
    return sum( @. abs((Gc - 𝐺.value) / 𝐺.error) )
end

function _calc_gf(Λ::Array{C64,2}, nk::I64)
    Ngrid, Kmax = size(Λ)
    @assert nk ≤ Kmax

    G = zeros(C64, Ngrid)
    for k = 1:nk
        for g = 1:Ngrid
            G[g] = G[g] + Λ[g,k]
        end
    end

    return G
end

function _calc_norm(C::Vector{Rectangle})
    norm = sum(map(x -> x.h * x.w, C))
    return norm
end

function Pdx(xmin::F64, xmax::F64, rng::AbstractRNG)
    γ = 2.0
    y = rand(rng, F64)

    _X = max(abs(xmin), abs(xmax))
    _λ = γ / _X
    _elx = exp(-1.0 * _λ * abs(xmin))
    _N = _λ / ( (xmin / abs(xmin)) * (exp(-1.0 * _λ * abs(xmin)) - 1.0)
              + (xmax / abs(xmax)) * (1.0 - exp(-1.0 * _λ * abs(xmax))) )
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
