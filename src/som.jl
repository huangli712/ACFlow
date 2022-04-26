#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/04/26
#

#=
### *Customized Structs* : *StochOM Solver*
=#

"""
    Box

Rectangle. The field configuration consists of many boxes. They exhibit
various areas (width × height). We used the Metropolis important sampling
algorithm to sample them and evaluate their contributions to the spectrum.

### Members

* h -> Height of the box.
* w -> Width of the box.
* c -> Position of the box.
"""
mutable struct Box
    h :: F64
    w :: F64
    c :: F64
end

"""
    StochOMElement

Mutable struct. It is used to record the field configurations, which will
be sampled by monte carlo procedure.

### Members

* C -> Field configuration.
* Λ -> Contributions of the field configuration to the correlator.
* G -> Reproduced correlator.
* Δ -> Difference between reproduced and raw correlators.
"""
mutable struct StochOMElement
    C :: Vector{Box}
    Λ :: Array{F64,2}
    G :: Vector{F64}
    Δ :: F64
end

"""
    StochOMContext

Mutable struct. It is used within the StochOM solver only.

### Members

* Gᵥ    -> Input data for correlator.
* σ¹    -> Actually 1.0 / σ¹.
* grid  -> Grid for input data.
* mesh  -> Mesh for output spectrum.
* Cᵥ    -> It is used to record the field configurations for all attempts.
* Δᵥ    -> It is used to record the errors for all attempts.
"""
mutable struct StochOMContext
    Gᵥ   :: Vector{F64}
    σ¹   :: Vector{F64}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    Cᵥ   :: Vector{Vector{Box}}
    Δᵥ   :: Vector{F64}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochOMSolver, rd::RawData)

Solve the analytical continuation problem by the stochastic optimization
method.
"""
function solve(S::StochOMSolver, rd::RawData)
    println("[ StochOM ]")
    MC, SC = init(S, rd)

    if nworkers() > 1
        println("Using $(nworkers()) workers")
        #
        p1 = deepcopy(PCOMM)
        p2 = deepcopy(PStochOM)
        #
        sol = pmap((x) -> prun(S, p1, p2, MC, SC), 1:nworkers())
        @assert length(sol) == nworkers()
        #
        Aout = similar(sol[end])
        fill!(Aout, 0.0)
        for i in eachindex(sol)
            @. Aout = Aout + sol[i] / nworkers()
        end
        #
        Gout = last(SC, Aout)
    else
        Aout = run(S, MC, SC)
        Gout = last(SC, Aout)
    end
    return SC.mesh, Aout, Gout
end

"""
    init(S::StochOMSolver, rd::RawData)

Initialize the StochOM solver and return the StochOMMC and StochOMContext
structs.
"""
function init(S::StochOMSolver, rd::RawData)
    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    Gᵥ, σ¹ = init_iodata(S, rd)
    println("Postprocess input data: ", length(σ¹), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    Cᵥ, Δᵥ = init_context(S)
    SC = StochOMContext(Gᵥ, σ¹, grid, mesh, Cᵥ, Δᵥ)

    return MC, SC
end

"""
    run(S::StochOMSolver, MC::StochOMMC, SC::StochOMContext)

Perform stochastic optimization simulation, sequential version.
"""
function run(S::StochOMSolver, MC::StochOMMC, SC::StochOMContext)
    ntry = get_s("ntry")
    nstep = get_s("nstep")

    for l = 1:ntry
        SE = init_element(MC, SC)

        for _ = 1:nstep
            update(MC, SE, SC)
        end

        SC.Δᵥ[l] = SE.Δ
        SC.Cᵥ[l] = deepcopy(SE.C)
        @printf("try -> %5i (%5i) Δ -> %8.4e \n", l, ntry, SE.Δ)
        (l % 10 == 0) && write_statistics(MC)
    end

    return average(SC)
end

"""
    prun(S::StochOMSolver,
         p1::Dict{String,Vector{Any}},
         p2::Dict{String,Vector{Any}},
         MC::StochOMMC, SC::StochOMContext)

Perform stochastic optimization simulation, parallel version.
The arguments `p1` and `p2` are copies of PCOMM and PStochOM, respectively.
"""
function prun(S::StochOMSolver,
              p1::Dict{String,Vector{Any}},
              p2::Dict{String,Vector{Any}},
              MC::StochOMMC, SC::StochOMContext)
    rev_dict(p1)
    rev_dict(S, p2)

    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    ntry = get_s("ntry")
    nstep = get_s("nstep")

    for l = 1:ntry
        SE = init_element(MC, SC)

        for _ = 1:nstep
            update(MC, SE, SC)
        end

        SC.Δᵥ[l] = SE.Δ
        SC.Cᵥ[l] = deepcopy(SE.C)
        @printf("try -> %5i (%5i) Δ -> %8.4e \n", l, ntry, SE.Δ)
        (myid() == 2) && (l % 10 == 0) && write_statistics(MC)
    end

    return average(SC)
end

"""
    average(SC::StochOMContext)

Postprocess the collected results after the stochastic optimization
simulations. It will calculate real spectral functions.
"""
function average(SC::StochOMContext)
    nmesh = get_c("nmesh")
    ntry  = get_s("ntry")

    # Calculate the median of SC.Δᵥ
    dev_ave = median(SC.Δᵥ)

    # Determine the αgood parameter, which is used to filter the
    # calculated spectra.
    αgood = 1.2
    if count(x -> x < dev_ave / αgood, SC.Δᵥ) == 0
        αgood = 1.0
    end

    # Accumulate the final spectrum
    Aom = zeros(F64, nmesh)
    for l = 1:ntry
        if SC.Δᵥ[l] < dev_ave / αgood
            for w = 1:nmesh
                _omega = SC.mesh[w]
                for r = 1:length(SC.Cᵥ[l])
                    R = SC.Cᵥ[l][r]
                    if R.c - 0.5 * R.w ≤ _omega ≤ R.c + 0.5 * R.w
                        Aom[w] = Aom[w] + R.h
                    end
                end
            end
        end
    end

    # Normalize the spectrum
    Lgood = count(x -> x < dev_ave / αgood, SC.Δᵥ)
    @. Aom = Aom / Lgood

    @printf("Median χ² : %16.12e Accepted configurations : %5i \n", dev_ave, Lgood)

    return Aom
end

"""
    last(SC::StochOMContext, Aout::Vector{F64})

It will process and write the calculated results by the StochOM solver,
including final spectral function and reproduced correlator.
"""
function last(SC::StochOMContext, Aout::Vector{F64})
    write_spectrum(SC.mesh, Aout)

    # Reproduce input data
    kernel = make_kernel(SC.mesh, SC.grid)
    G = reprod(kernel, SC.mesh, Aout)
    write_backward(SC.grid, G)

    _G = kramers(SC.mesh, Aout)
    write_complete(SC.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

"""
    update(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext)

Using the Metropolis algorithm to update the field configuration, i.e, a
collection of hundreds of boxes.
"""
function update(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext)
    Tmax = 100
    nbox = get_s("nbox")

    T1 = rand(MC.rng, 1:Tmax)
    d1 = rand(MC.rng, F64)
    d2 = 1.0 + rand(MC.rng, F64)

    ST = deepcopy(SE)

    for _ = 1:T1
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < nbox - 1
                    try_insert(MC, ST, SC, d1)
                end
                break

            @case 2
                if length(ST.C) > 1
                    try_remove(MC, ST, SC, d1)
                end
                break

            @case 3
                try_shift(MC, ST, SC, d1)
                break

            @case 4
                try_width(MC, ST, SC, d1)
                break

            @case 5
                if length(ST.C) > 1
                    try_height(MC, ST, SC, d1)
                end
                break

            @case 6
                if length(ST.C) < nbox - 1
                    try_split(MC, ST, SC, d1)
                end
                break

            @case 7
                if length(ST.C) > 1
                    try_merge(MC, ST, SC, d1)
                end
                break
        end

    end

    for _ = T1+1:Tmax
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < nbox - 1
                    try_insert(MC, ST, SC, d2)
                end
                break

            @case 2
                if length(ST.C) > 1
                    try_remove(MC, ST, SC, d2)
                end
                break

            @case 3
                try_shift(MC, ST, SC, d2)
                break

            @case 4
                try_width(MC, ST, SC, d2)
                break

            @case 5
                if length(ST.C) > 1
                    try_height(MC, ST, SC, d2)
                end
                break

            @case 6
                if length(ST.C) < nbox - 1
                    try_split(MC, ST, SC, d2)
                end
                break

            @case 7
                if length(ST.C) > 1
                    try_merge(MC, ST, SC, d2)
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

#=
### *Service Functions*
=#

"""
    init_mc(S::StochOMSolver)

Try to create a StochOMMC struct.

See also: [`StochOMMC`](@ref).
"""
function init_mc(S::StochOMSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    Macc = zeros(I64, 7)
    Mtry = zeros(I64, 7)

    MC = StochOMMC(rng, Macc, Mtry)

    return MC
end

"""
    init_element(MC::StochOMMC, SC::StochOMContext)

Try to initialize a StochOMElement struct.

See also: [`StochOMElement`](@ref).
"""
function init_element(MC::StochOMMC, SC::StochOMContext)
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    nbox = get_s("nbox")
    sbox = get_s("sbox")
    wbox = get_s("wbox")

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
    Λ = zeros(F64, length(SC.Gᵥ), nbox)
    Δ = 0.0

    for k = 1:_Know
        c = wmin + wbox / 2.0 + (wmax - wmin - wbox) * rand(MC.rng, F64)
        w = wbox + (min(2.0 * (c - wmin), 2.0 * (wmax - c)) - wbox) * rand(MC.rng, F64)
        h = weight[k] / w
        R = Box(h, w, c)
        push!(C, R)
        Λ[:,k] .= calc_lambda(R, SC.grid)
    end
    G = calc_green(Λ, _Know)
    Δ = calc_error(G, SC.Gᵥ, SC.σ¹)

    return StochOMElement(C, Λ, G, Δ)
end

"""
    init_context(S::StochOMSolver)

Try to initialize the key members of a StochOMContext struct.

See also: [`StochOMContext`](@ref).
"""
function init_context(S::StochOMSolver)
    ntry = get_s("ntry")
    nbox = get_s("nbox")

    Δv = zeros(F64, ntry)

    Cv = []
    for _ = 1:ntry
        C = Box[]
        for _ = 1:nbox
            push!(C, Box(0.0, 0.0, 0.0))
        end
        push!(Cv, C)
    end

    return Cv, Δv
end

"""
    init_iodata(S::StochACSolver, rd::RawData)

Preprocess the input data (`rd`).

See also: [`RawData`](@ref), [`GreenData`](@ref).
"""
function init_iodata(S::StochOMSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value
    σ¹ = 1.0 ./ G.error

    return Gᵥ, σ¹
end

"""
    calc_lambda(r::Box, grid::FermionicMatsubaraGrid)

Try to calculate the kernel-related function Λ. This function works for
FermionicMatsubaraGrid only.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function calc_lambda(r::Box, grid::FermionicMatsubaraGrid)
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w
    iw = im * grid.ω
    Λ = @. r.h * log((iw - e₁) / (iw - e₂))
    return vcat(real(Λ), imag(Λ))
end

"""
    calc_lambda(r::Box, grid::BosonicMatsubaraGrid)

Try to calculate the kernel-related function Λ. This function works for
BosonicMatsubaraGrid only.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function calc_lambda(r::Box, grid::BosonicMatsubaraGrid)
    ktype = get_c("ktype")

    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    if ktype == "bsymm"
        Λ = @. atan( e₁ / grid.ω ) - atan( e₂ / grid.ω )
        Λ = r.h * (r.w .+ grid.ω .* Λ)
        return Λ
    else
        iw = im * grid.ω
        Λ = @. r.h * log((iw - e₁) / (iw - e₂))
        return vcat(real(Λ), imag(Λ))
    end
end

"""
    calc_error(G::Vector{F64}, Gᵥ::Vector{F64}, σ¹::Vector{F64})

Try to calculate χ². Here `Gᵥ` and `σ¹` denote the raw correlator and
related standard deviation. `G` means the reproduced correlator.
"""
function calc_error(G::Vector{F64}, Gᵥ::Vector{F64}, σ¹::Vector{F64})
    return sum( abs.((G .- Gᵥ) .* σ¹) )
end

"""
    calc_green(Λ::Array{F64,2}, nk::I64)

Try to reconstruct the correlator via the field configuration.
"""
function calc_green(Λ::Array{F64,2}, nk::I64)
    ngrid, nbox = size(Λ)
    @assert nk ≤ nbox

    G = zeros(F64, ngrid)
    for k = 1:nk
        for g = 1:ngrid
            G[g] = G[g] + Λ[g,k]
        end
    end

    return G
end

"""
    calc_norm(C::Vector{Box})

Calculate the total area of all boxes.
"""
function calc_norm(C::Vector{Box})
    norm = sum(map(x -> x.h * x.w, C))
    return norm
end

"""
    try_insert(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Insert a new box into the field configuration.
"""
function try_insert(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    sbox  = get_s("sbox")
    wbox  = get_s("wbox")
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    csize = length(SE.C)

    t = rand(MC.rng, 1:csize)

    R = SE.C[t]
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

    G1 = SE.Λ[:,t]
    G2 = calc_lambda(Rnew, SC.grid)
    G3 = calc_lambda(Radd, SC.grid)

    Δ = calc_error(SE.G - G1 + G2 + G3, SC.Gᵥ, SC.σ¹)

    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        SE.C[t] = Rnew
        push!(SE.C, Radd)
        SE.Δ = Δ
        @. SE.G = SE.G - G1 + G2 + G3
        @. SE.Λ[:,t] = G2
        @. SE.Λ[:,csize+1] = G3
        MC.Macc[1] = MC.Macc[1] + 1
    end

    MC.Mtry[1] = MC.Mtry[1] + 1
end

"""
    try_remove(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Remove an old box from the field configuration.
"""
function try_remove(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    csize = length(SE.C)

    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end
    if t1 < t2
        t1, t2 = t2, t1
    end

    R1 = SE.C[t1]
    R2 = SE.C[t2]
    Re = SE.C[end]

    dx = R1.h * R1.w

    G1 = SE.Λ[:,t1]
    G2 = SE.Λ[:,t2]
    Ge = SE.Λ[:,csize]

    R2n = Box(R2.h + dx / R2.w, R2.w, R2.c)
    G2n = calc_lambda(R2n, SC.grid)

    Δ = calc_error(SE.G - G1 - G2 + G2n, SC.Gᵥ, SC.σ¹)

    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        SE.C[t2] = R2n
        if t1 < csize
            SE.C[t1] = Re
        end
        pop!(SE.C)
        SE.Δ = Δ
        @. SE.G = SE.G - G1 - G2 + G2n
        @. SE.Λ[:,t2] = G2n
        if t1 < csize
            @. SE.Λ[:,t1] = Ge
        end
        MC.Macc[2] = MC.Macc[2] + 1
    end

    MC.Mtry[2] = MC.Mtry[2] + 1
end

"""
    try_shift(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Change the position of given box in the field configuration.
"""
function try_shift(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    csize = length(SE.C)

    t = rand(MC.rng, 1:csize)

    R = SE.C[t]

    dx_min = wmin + R.w / 2.0 - R.c
    dx_max = wmax - R.w / 2.0 - R.c
    if dx_max ≤ dx_min
        return
    end
    dc = Pdx(dx_min, dx_max, MC.rng)

    Rn = Box(R.h, R.w, R.c + dc)
    G1 = SE.Λ[:,t]
    G2 = calc_lambda(Rn, SC.grid)

    Δ = calc_error(SE.G - G1 + G2, SC.Gᵥ, SC.σ¹)

    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        SE.C[t] = Rn
        SE.Δ = Δ
        @. SE.G = SE.G - G1 + G2
        @. SE.Λ[:,t] = G2
        MC.Macc[3] = MC.Macc[3] + 1
    end

    MC.Mtry[3] = MC.Mtry[3] + 1
end

"""
    try_width(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Change the width of given box in the field configuration.
"""
function try_width(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wbox  = get_s("wbox")
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    csize = length(SE.C)

    t = rand(MC.rng, 1:csize)

    R = SE.C[t]

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
    G1 = SE.Λ[:,t]
    G2 = calc_lambda(Rn, SC.grid)

    Δ = calc_error(SE.G - G1 + G2, SC.Gᵥ, SC.σ¹)

    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        SE.C[t] = Rn
        SE.Δ = Δ
        @. SE.G = SE.G - G1 + G2
        @. SE.Λ[:,t] = G2
        MC.Macc[4] = MC.Macc[4] + 1
    end

    MC.Mtry[4] = MC.Mtry[4] + 1
end

"""
    try_height(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Change the height of given box in the field configuration.
"""
function try_height(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    sbox  = get_s("sbox")
    csize = length(SE.C)

    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end

    R1 = SE.C[t1]
    R2 = SE.C[t2]

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
    G1A = SE.Λ[:,t1]
    G1B = calc_lambda(R1n, SC.grid)
    R2n = Box(R2.h - dh * w1 / w2, R2.w, R2.c)
    G2A = SE.Λ[:,t2]
    G2B = calc_lambda(R2n, SC.grid)

    Δ = calc_error(SE.G - G1A + G1B - G2A + G2B, SC.Gᵥ, SC.σ¹)

    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        SE.C[t1] = R1n
        SE.C[t2] = R2n
        SE.Δ = Δ
        @. SE.G = SE.G - G1A + G1B - G2A + G2B
        @. SE.Λ[:,t1] = G1B
        @. SE.Λ[:,t2] = G2B
        MC.Macc[5] = MC.Macc[5] + 1
    end

    MC.Mtry[5] = MC.Mtry[5] + 1
end

"""
    try_split(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Split a given box into two boxes in the field configuration.
"""
function try_split(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wbox  = get_s("wbox")
    sbox  = get_s("sbox")
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    csize = length(SE.C)

    t = rand(MC.rng, 1:csize)

    R1 = SE.C[t]
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

        G1 = SE.Λ[:,t]
        Ge = SE.Λ[:,csize]

        R2 = Box(h, w1, c1 + dc1)
        G2 = calc_lambda(R2, SC.grid)

        R3 = Box(h, w2, c2 + dc2)
        G3 = calc_lambda(R3, SC.grid)
        Δ = calc_error(SE.G - G1 + G2 + G3, SC.Gᵥ, SC.σ¹)

        if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
            SE.C[t] = SE.C[end]
            pop!(SE.C)
            push!(SE.C, R2)
            push!(SE.C, R3)
            SE.Δ = Δ
            @. SE.G = SE.G - G1 + G2 + G3
            if t < csize
                @. SE.Λ[:,t] = Ge
            end
            @. SE.Λ[:,csize] = G2
            @. SE.Λ[:,csize+1] = G3
            MC.Macc[6] = MC.Macc[6] + 1
        end
    end

    MC.Mtry[6] = MC.Mtry[6] + 1
end

"""
    try_merge(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Merge two given boxes into one box in the field configuration.
"""
function try_merge(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    csize = length(SE.C)

    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end
    if t1 > t2
        t1, t2 = t2, t1
    end

    R1 = SE.C[t1]
    R2 = SE.C[t2]

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

    G1 = SE.Λ[:,t1]
    G2 = SE.Λ[:,t2]
    Ge = SE.Λ[:,csize]

    Rn = Box(h_new, w_new, c_new + dc)
    Gn = calc_lambda(Rn, SC.grid)

    Δ = calc_error(SE.G - G1 - G2 + Gn, SC.Gᵥ, SC.σ¹)

    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        SE.C[t1] = Rn
        if t2 < csize
            SE.C[t2] = SE.C[end]
        end
        pop!(SE.C)
        SE.Δ = Δ
        @. SE.G = SE.G - G1 - G2 + Gn
        @. SE.Λ[:,t1] = Gn
        if t2 < csize
            @. SE.Λ[:,t2] = Ge
        end
        MC.Macc[7] = MC.Macc[7] + 1
    end

    MC.Mtry[7] = MC.Mtry[7] + 1
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
