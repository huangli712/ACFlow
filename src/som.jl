#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/04/03
#

function read_data!()
    grid  = F64[]
    value = C64[]
    error = C64[]

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
        end
    end

    β = 2.0 * π / (grid[2] - grid[1])
    return FermionicMatsubaraGrid(niw, β, grid), SOMData(value, error)
end

struct SOMData <: AbstractData
    value :: Vector{N64}
    error :: Vector{N64}
end

const P_SOM = Dict{String, Any}(
    "ngrid" => 64,
    "nmesh" => 501,
    "nstep" => 4000,
    "ntry"  => 1000,
    "nbox"  => 100,
    "sbox"  => 0.005,
    "wbox"  => 0.02,
    "dmax"  => 2.0,
    "wmax"  => 10.0,
    "wmin"  => -10.0,
    "alpha" => 10.0,
    "norm"  => -1.0,
)

mutable struct Box
    h :: F64
    w :: F64
    c :: F64
end

mutable struct SOMElement
    C :: Vector{Box}
    Λ :: Array{C64,2}
    G :: Vector{C64}
    Δ :: F64
end

mutable struct SOMContext
    Cv :: Vector{Vector{Box}}
    Δv :: Vector{F64}
end

function som_run(ω::FermionicMatsubaraGrid, 𝐺::SOMData)
    nstep = P_SOM["nstep"]
    ntry = P_SOM["ntry"]

    SC, MC = som_init()

    for l = 1:nstep
        println("try: $l")

        SE = som_random(MC, ω, 𝐺)

        for _ = 1:ntry
            som_update(SE, MC, ω, 𝐺)
        end

        SC.Δv[l] = SE.Δ
        SC.Cv[l] = deepcopy(SE.C)
    end

    return som_spectra(SC)
end

function som_init()
    nstep = P_SOM["nstep"]
    nbox = P_SOM["nbox"]

    Δv = zeros(F64, nstep)

    Cv = []
    for _ = 1:nstep
        C = Box[]
        for _ = 1:nbox
            push!(C, Box(0.0, 0.0, 0.0))
        end
        push!(Cv, C)
    end

    seed = rand(1:1000000)#;  seed = 112414
    rng = MersenneTwister(seed)
    @show "seed: ", seed
    Mtry = zeros(I64, 7)
    Macc = zeros(I64, 7)

    return SOMContext(Cv, Δv), StochOMMC(rng, Mtry, Macc)
end

function som_random(MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData)
    sbox  = P_SOM["sbox"]
    wbox  = P_SOM["wbox"]
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]
    nbox  = P_SOM["nbox"]
    ngrid = P_SOM["ngrid"]

    _Know = rand(MC.rng, 2:nbox)
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
    while weight[plus_count] < sbox
        while weight[minus_count] < 2 * sbox
            minus_count = minus_count - 1
        end
        weight[plus_count] = weight[plus_count] + sbox
        weight[minus_count] = weight[minus_count] - sbox
        plus_count = plus_count + 1
    end

    C = Box[]
    Λ = zeros(C64, ngrid, nbox)
    Δ = 0.0

    for k = 1:_Know
        c = wmin + wbox / 2.0 + (wmax - wmin - wbox) * rand(MC.rng, F64)
        w = wbox + (min(2.0 * (c - wmin), 2.0 * (wmax - c)) - wbox) * rand(MC.rng, F64)
        h = weight[k] / w
        R = Box(h, w, c)
        push!(C, R)
        Λ[:,k] .= _calc_lambda(R, ω)
    end
    Δ = _calc_err(Λ, _Know, 𝐺)
    G = _calc_gf(Λ, _Know)

    return SOMElement(C, Λ, G, Δ)
end

function som_update(SE::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData)
    Tmax = 100
    nbox = P_SOM["nbox"]
    dmax = P_SOM["dmax"]

    T1 = rand(MC.rng, 1:Tmax)
    d1 = rand(MC.rng, F64)
    d2 = 1.0 + (dmax - 1.0) * rand(MC.rng, F64)

    ST = deepcopy(SE)

    for _ = 1:T1
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < nbox - 1
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
                if length(ST.C) < nbox - 1
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
                if length(ST.C) < nbox - 1
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
                if length(ST.C) < nbox - 1
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
    nmesh = P_SOM["nmesh"]
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]
    nstep  = P_SOM["nstep"]

    dev_min = minimum(𝑆.Δv)

    Lgood = 0
    Aom = zeros(F64, nmesh)
    for l = 1:nstep
        if alpha * dev_min - 𝑆.Δv[l] > 0
            Lgood = Lgood + 1
            for w = 1:nmesh
                _omega = wmin + (w - 1) * (wmax - wmin) / (nmesh - 1)
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
    nmesh = P_SOM["nmesh"]
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]

    open("Aw.data", "w") do fout
        for w = 1:nmesh
            _omega = wmin + (w - 1) * (wmax - wmin) / (nmesh - 1)
            println(fout, _omega, " ", Aom[w])
        end
    end
end

function _try_insert(𝑆::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData, dacc)
    sbox  = P_SOM["sbox"]
    wbox  = P_SOM["wbox"]
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]
    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]
    if R.h * R.w ≤ 2.0 * sbox
        return
    end

    dx_min = sbox
    dx_max = R.h * R.w - sbox
    if dx_max ≤ dx_min
        return
    end
    r1 = rand(MC.rng, F64)
    r2 = rand(MC.rng, F64)
    c = (wmin + wbox / 2.0) + (wmax - wmin - wbox) * r1
    w_new_max = 2.0 * min(wmax - c, c - wmin)
    dx = Pdx(dx_min, dx_max, MC.rng)
    h = dx / w_new_max + (dx / wbox - dx / w_new_max) * r2
    w = dx / h

    Rnew = Box(R.h - dx / R.w, R.w, R.c)
    Radd = Box(h, w, c)

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
        MC.Macc[1] = MC.Macc[1] + 1
    end

    MC.Mtry[1] = MC.Mtry[1] + 1
end

function _try_remove(𝑆::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData, dacc)
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

    R2n = Box(R2.h + dx / R2.w, R2.w, R2.c)
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
        MC.Macc[2] = MC.Macc[2] + 1
    end

    MC.Mtry[2] = MC.Mtry[2] + 1
end

function _try_position(𝑆::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData, dacc)
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]
    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]

    dx_min = wmin + R.w / 2.0 - R.c
    dx_max = wmax - R.w / 2.0 - R.c
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, MC.rng)

    Rn = Box(R.h, R.w, R.c + dc)
    G1 = 𝑆.Λ[:,t]
    G2 = _calc_lambda(Rn, ω)

    Δ = _calc_err(𝑆.G - G1 + G2, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t] = Rn
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1 + G2
        @. 𝑆.Λ[:,t] = G2
        MC.Macc[3] = MC.Macc[3] + 1
    end

    MC.Mtry[3] = MC.Mtry[3] + 1
end

function _try_width(𝑆::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData, dacc)
    wbox  = P_SOM["wbox"]
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]
    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R = 𝑆.C[t]

    weight = R.h * R.w
    dx_min = wbox - R.w
    dx_max = min(2.0 * (R.c - wmin), 2.0 * (wmax - R.c)) - R.w
    if dx_max ≤ dx_min
        return
    end
    dw = Pdx(dx_min, dx_max, MC.rng)
    w = R.w + dw
    h = weight / w
    c = R.c

    Rn = Box(h, w, c)
    G1 = 𝑆.Λ[:,t]
    G2 = _calc_lambda(Rn, ω)

    Δ = _calc_err(𝑆.G - G1 + G2, 𝐺)

    if rand(MC.rng, F64) < ((𝑆.Δ/Δ) ^ (1.0 + dacc))
        𝑆.C[t] = Rn
        𝑆.Δ = Δ
        @. 𝑆.G = 𝑆.G - G1 + G2
        @. 𝑆.Λ[:,t] = G2
        MC.Macc[4] = MC.Macc[4] + 1
    end

    MC.Mtry[4] = MC.Mtry[4] + 1
end

function _try_height(𝑆::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData, dacc)
    sbox  = P_SOM["sbox"]
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
    dx_min = sbox / w1 - h1
    dx_max = (h2 - sbox / w2) * w2 / w1
    if dx_max ≤ dx_min
        return
    end
    dh = Pdx(dx_min, dx_max, MC.rng)

    R1n = Box(R1.h + dh, R1.w, R1.c)
    G1A = 𝑆.Λ[:,t1]
    G1B = _calc_lambda(R1n, ω)
    R2n = Box(R2.h - dh * w1 / w2, R2.w, R2.c)
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
        MC.Macc[5] = MC.Macc[5] + 1
    end

    MC.Mtry[5] = MC.Mtry[5] + 1
end

function _try_split(𝑆::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData, dacc)
    wbox  = P_SOM["wbox"]
    sbox  = P_SOM["sbox"]
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]
    csize = length(𝑆.C)

    t = rand(MC.rng, 1:csize)

    R1 = 𝑆.C[t]
    if R1.w ≤ 2 * wbox || R1.w * R1.h ≤ 2.0 * sbox
        return
    end

    h = R1.h
    w1 = wbox + (R1.w - 2.0 * wbox) * rand(MC.rng, F64)
    w2 = R1.w - w1
    if w1 > w2
        w1, w2 = w2, w1
    end
    c1 = R1.c - R1.w / 2.0 + w1 / 2.0
    c2 = R1.c + R1.w / 2.0 - w2 / 2.0
    dx_min = wmin + w1 / 2.0 - c1
    dx_max = wmax - w1 / 2.0 - c1
    if dx_max ≤ dx_min
        return
    end
    dc1 = Pdx(dx_min, dx_max, MC.rng)
    dc2 = -1.0 * w1 * dc1 / w2

    if (c1 + dc1 ≥ wmin + w1 / 2.0) &&
       (c1 + dc1 ≤ wmax - w1 / 2.0) &&
       (c2 + dc2 ≥ wmin + w2 / 2.0) &&
       (c2 + dc2 ≤ wmax - w2 / 2.0)

        G1 = 𝑆.Λ[:,t]
        Ge = 𝑆.Λ[:,csize]

        R2 = Box(h, w1, c1 + dc1)
        G2 = _calc_lambda(R2, ω)

        R3 = Box(h, w2, c2 + dc2)
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
            MC.Macc[6] = MC.Macc[6] + 1
        end
    end

    MC.Mtry[6] = MC.Mtry[6] + 1
end

function _try_merge(𝑆::SOMElement, MC::StochOMMC, ω::FermionicMatsubaraGrid, 𝐺::SOMData, dacc)
    wmin = P_SOM["wmin"]
    wmax = P_SOM["wmax"]
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
    dx_min = wmin + w_new / 2.0 - c_new
    dx_max = wmax - w_new / 2.0 - c_new
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, MC.rng)

    G1 = 𝑆.Λ[:,t1]
    G2 = 𝑆.Λ[:,t2]
    Ge = 𝑆.Λ[:,csize]

    Rn = Box(h_new, w_new, c_new + dc)
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
        MC.Macc[7] = MC.Macc[7] + 1
    end

    MC.Mtry[7] = MC.Mtry[7] + 1
end

function _calc_lambda(r::Box, ω::FermionicMatsubaraGrid)
    Λ = @. r.h * log((im * ω.ω - r.c + 0.5 * r.w) / (im * ω.ω - r.c - 0.5 * r.w))
    return Λ
end

function _calc_err(Λ::Array{C64,2}, nk::I64, 𝐺::SOMData)
    ngrid, nbox = size(Λ)
    @assert nk ≤ nbox

    res = 0.0
    for w = 1:ngrid
        g = sum(Λ[w,1:nk])
        res = res + abs((g - 𝐺.value[w]) / 𝐺.error[w])
    end

    return res
end

function _calc_err(Gc::Vector{C64}, 𝐺::SOMData)
    return sum( @. abs((Gc - 𝐺.value) / 𝐺.error) )
end

function _calc_gf(Λ::Array{C64,2}, nk::I64)
    ngrid, nbox = size(Λ)
    @assert nk ≤ nbox

    G = zeros(C64, ngrid)
    for k = 1:nk
        for g = 1:ngrid
            G[g] = G[g] + Λ[g,k]
        end
    end

    return G
end

function _calc_norm(C::Vector{Box})
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
