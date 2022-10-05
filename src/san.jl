#
# Project : Gardenia
# Source  : san.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/10/05
#

#=
### *Customized Structs* : *StochSK Solver*
=#

"""
    StochSKElement
"""
mutable struct StochSKElement
    P :: Vector{I64}
    A :: F64
    W :: I64
end

"""
    StochACContext
"""
mutable struct StochSKContext
    Gᵥ     :: Vector{F64}
    Gᵧ     :: Vector{F64}
    σ¹     :: Vector{F64}
    grid   :: AbstractGrid
    mesh   :: AbstractMesh
    kernel :: Array{F64,2}
    Aout   :: Vector{F64}
    χ²     :: F64
    χ²min  :: F64
    χ²vec  :: Vector{F64}
    Θ      :: F64
    Θvec   :: Vector{F64}
    𝒞ᵧ     :: Vector{StochSKElement}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochSKSolver, rd::RawData)
"""
function solve(S::StochSKSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    mc = init_mc()
    fmesh = LinearMesh(get_k("nfine"), get_b("wmin"), get_b("wmax"))
    grid = make_grid(rd)
    kernel = init_kernel(grid.τ, fmesh)
    SE = init_element(mc.rng)

    Gᵧ = calc_correlator(SE, kernel)

    #
    mesh = make_mesh()
    Aout = zeros(F64, get_b("nmesh"))
    #
    χ = calc_goodness(Gᵧ, Gᵥ, σ¹)
    χ² = χ
    χ²min = χ
    χ²vec = zeros(F64, get_k("nwarm"))
    #
    Θ = get_k("theta")
    Θvec = zeros(F64, get_k("nwarm"))
    #
    𝒞ᵧ = StochSKElement[]
    #
    SC = StochSKContext(Gᵥ, Gᵧ, σ¹, grid, mesh, kernel, Aout, χ², χ²min, χ²vec, Θ, Θvec, 𝒞ᵧ)

    warmup(mc, SE, SC)
    SE = analyze(SC)
    sample(mc, SE, SC)
end

function init()
end

function run()
end

function prun()
end

function average()
end

function last()
end

function warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    anneal_length = get_k("nwarm")

    for i = 1:anneal_length
        try_update(MC, SE, SC)

        push!(SC.𝒞ᵧ, deepcopy(SE))
        SC.χ²vec[i] = SC.χ²
        SC.Θvec[i] = SC.Θ

        @show i, SC.χ², SC.χ²min, SC.χ² - SC.χ²min
        if SC.χ² - SC.χ²min < 1e-3
            break
        end

        SC.Θ = SC.Θ * get_k("ratio")
    end
end

function analyze(SC::StochSKContext)
    num_anneal = length(SC.𝒞ᵧ)
    @assert num_anneal ≤ get_k("nwarm")

    c = num_anneal
    while c ≥ 1
        if SC.χ²vec[c] > SC.χ²min + 2.0 * sqrt(SC.χ²min)
            break
        end
        c = c - 1
    end
    @assert 1 ≤ c ≤ num_anneal

    SE = deepcopy(SC.𝒞ᵧ[c])
    SC.Θ = SC.Θvec[c]
    SC.Gᵧ = calc_correlator(SE, SC.kernel)
    SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
    @show SC.Θ, SC.χ²

    return SE
end

function sample(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    nmesh = get_b("nmesh")
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")
    nstep = get_k("nstep")
    retry = get_k("retry")

    try_update(MC, SE, SC)

    for i = 1:nstep
        if i % retry == 0
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
            @show i, SC.χ²
        end

        try_update_s(MC, SE, SC)

        for j = 1:ngamm
            d_pos = SE.P[j]
            s_pos = ceil(I64, d_pos / nfine * nmesh)
            SC.Aout[s_pos] = SC.Aout[s_pos] + SE.A
        end
    end

    SC.Aout = SC.Aout / (nstep * (SC.mesh[2] - SC.mesh[1]))

    open("Aout.data", "w") do fout
        for i in eachindex(SC.mesh)
            println(fout, SC.mesh[i], " ", SC.Aout[i])
        end
    end
end

function init_mc()
    seed = rand(1:1000000); seed = 840443
    rng = MersenneTwister(seed)

    Sacc = 0
    Stry = 0
    Pacc = 0
    Ptry = 0
    MC = StochSKMC(rng, Sacc, Stry, Pacc, Ptry)

    return MC
end

function init_iodata()
end

function init_kernel(tmesh, fmesh::AbstractMesh)
    beta = get_b("beta")
    nfine = get_k("nfine")

    ntau = length(tmesh)
    kernel = zeros(F64, ntau, nfine)

    for f = 1:nfine
        ω = fmesh[f]
        de = 1.0 + exp(-beta * ω)
        kernel[:,f] = exp.(-ω * tmesh) / de
    end

    return kernel
end

function init_element(rng)
    β = get_b("beta")
    wmax = get_b("wmax")
    wmin = get_b("wmin")
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    position = zeros(I64, ngamm)
    rand!(rng, position, 1:nfine)
    #
    amplitude = 1.0 / ngamm
    #
    δf = (wmax - wmin) / (nfine - 1)
    average_freq = abs(log(2.0) / β)
    window_width = ceil(I64, 0.1 * average_freq / δf)

    return StochSKElement(position, amplitude, window_width)
end

function calc_correlator(SE::StochSKElement, kernel::Array{F64,2})
    ngamm = length(SE.P)
    𝐴 = fill(SE.A, ngamm)
    𝐾 = kernel[:, SE.P]
    return 𝐾 * 𝐴
end

function calc_goodness(Gₙ::Vector{F64,}, Gᵥ::Vector{F64}, σ¹::Vector{F64})
    χ = sum( ( (Gₙ .- Gᵥ) .* σ¹ ) .^ 2.0 )
    return χ
end

function try_update(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    nfine = get_k("nfine")
    retry = get_k("retry")
    max_bin_size = 100

    bin_χ²  = zeros(F64, max_bin_size)
    bin_acc = zeros(I64, max_bin_size)
    bin_try = zeros(I64, max_bin_size)

    for s = 1:max_bin_size
        if s % retry == 0
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
        end

        try_update_s(MC, SE, SC)

        bin_χ²[s]  = SC.χ²
        bin_acc[s] = MC.Sacc + MC.Pacc
        bin_try[s] = MC.Stry + MC.Ptry
    end

    𝑝 = sum(bin_acc) / sum(bin_try)
    #
    if 𝑝 > 0.5
        r = SE.W * 1.5
        if ceil(I64, r) < nfine
            SE.W = ceil(I64, r)
        else
            SE.W = nfine
        end
    end
    #
    if 𝑝 < 0.4
        SE.W = ceil(I64, SE.W / 1.5)
    end

    SC.χ² = mean(bin_χ²)
end

function try_update_s(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    MC.Sacc = 0
    MC.Stry = ngamm
    @assert 1 < SE.W ≤ nfine

    for i = 1:ngamm
        s = rand(MC.rng, 1:ngamm)
        pcurr = SE.P[s]

        if 1 < SE.W < nfine
            move_width = rand(MC.rng, 1:SE.W)

            if rand(MC.rng) > 0.5
                pnext = pcurr + move_width
            else
                pnext = pcurr - move_width
            end

            pnext < 1     && (pnext = pnext + nfine)
            pnext > nfine && (pnext = pnext - nfine)
        else
            pnext = rand(MC.rng, 1:nfine)
        end

        Knext = view(SC.kernel, :, pnext)
        Kcurr = view(SC.kernel, :, pcurr)
        Gₙ = SC.Gᵧ + SE.A * (Knext - Kcurr)
        χ²new = calc_goodness(Gₙ, SC.Gᵥ, SC.σ¹)
        prob = exp( 0.5 * (SC.χ² - χ²new) / SC.Θ )

        if rand(MC.rng) < min(prob, 1.0)
            SE.P[s] = pnext
            SC.Gᵧ = Gₙ

            SC.χ² = χ²new
            if χ²new < SC.χ²min
                SC.χ²min = χ²new
            end

            MC.Sacc = MC.Sacc + 1
        end
    end
end

function try_update_p()
end
