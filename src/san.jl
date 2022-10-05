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

Mutable struct. It is used to record the field configurations, which will
be sampled within monte carlo procedure.

### Members

* P -> It means the positions of the δ functions.
* A -> It means the weights / amplitudes of the δ functions.
* W -> It denotes the window that is used to restrict the δ functions.
"""
mutable struct StochSKElement
    P :: Vector{I64}
    A :: F64
    W :: I64
end

"""
    StochSKContext

Mutable struct. It is used within the StochSK solver only.

### Members

* Gᵥ     -> Input data for correlator.
* σ¹     -> Actually 1.0 / σ¹.
* allow  -> Allowable indices.
* grid   -> Grid for input data.
* mesh   -> Mesh for output spectrum.
* model  -> Default model function.
* kernel -> Default kernel function.
* Aout   -> Calculate spectrum, it is actually ⟨n(x)⟩.
* Δ      -> Precomputed δ functions.
* hτ     -> α-resolved h(τ).
* Hα     -> α-resolved Hc.
* Uα     -> α-resolved internal energy, it is actually ⟨Hα⟩.
* αₗ     -> Vector of the α parameters.
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
    nmesh = get_b("nmesh")

    println("[ StochSK ]")
    MC, SE, SC = init(S, rd)

    if nworkers() > 1
        println("Using $(nworkers()) workers")
        #
        # Copy configuration dicts
        p1 = deepcopy(PBASE)
        p2 = deepcopy(PStochSK)
        #
            #sol = pmap((x) -> prun(S, p1, p2, MC, SE, SC), 1:nworkers())
            #@assert length(sol) == nworkers()
        #
        # Launch the task
        𝐹 = Future[]
        for i = 1:nworkers()
            𝑓 = @spawnat i + 1 prun(S, p1, p2, MC, SE, SC)
            push!(𝐹, 𝑓)
        end
        #
        # Wait and collect the solutions
        sol = []
        for i = 1:nworkers()
            wait(𝐹[i])
            push!(sol, fetch(𝐹[i]))
        end
        #
        # Average the solutions
        Aout = zeros(F64, nmesh)
        for i in eachindex(sol)
            a = sol[i]
            @. Aout = Aout + a / nworkers()
        end
        #
        # Postprocess the solutions
        Gout = last(SC, Aout)
    else
        Aout = run(MC, SE, SC)
        Gout = last(SC, Aout)    
    end

    return SC.mesh.mesh, Aout, Gout
end

"""
    init(S::StochSKSolver, rd::RawData)

Initialize the StochSK solver and return the StochSKMC, StochSKElement,
and StochSKContext structs.
"""
function init(S::StochSKSolver, rd::RawData)
    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    SE = init_element(S, MC.rng)
    println("Randomize Monte Carlo configurations")

    Gᵥ, σ¹, Aout = init_iodata(S, rd)
    println("Postprocess input data: ", length(σ¹), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    fmesh = calc_fmesh(S)
    kernel = make_kernel(fmesh, grid)
    println("Build default kernel: ", get_b("ktype"))

    Gᵧ = calc_correlator(SE, kernel)

    χ = calc_goodness(Gᵧ, Gᵥ, σ¹)
    χ², χ²min = χ, χ
    χ²vec = zeros(F64, get_k("nwarm"))

    Θ = get_k("theta")
    Θvec = zeros(F64, get_k("nwarm"))

    𝒞ᵧ = StochSKElement[]
    #
    SC = StochSKContext(Gᵥ, Gᵧ, σ¹, grid, mesh, kernel, Aout, χ², χ²min, χ²vec, Θ, Θvec, 𝒞ᵧ)

    return MC, SE, SC
end

"""
    run(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Perform stochastic analytical continuation simulation, sequential version.
"""
function run(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    nstep = get_k("nstep")
    retry = get_k("retry")
    output_per_steps = get_k("ndump")
    measure_per_steps = 10

    println("Start thermalization...")
    SE = warmup(MC, SE, SC)

    shuffle(MC, SE, SC)

    step = 0.0
    for iter = 1:nstep
        if iter % retry == 0
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
        end

        sample(MC, SE, SC)

        if iter % measure_per_steps == 0
            step = step + 1.0
            measure(SE, SC)
        end

        if iter % output_per_steps == 0
            prog = round(I64, iter / nstep * 100)
            println("Start stochastic sampling (prog: $prog)")
            flush(stdout)
            write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    prun(S::StochSKSolver,
         p1::Dict{String,Vector{Any}},
         p2::Dict{String,Vector{Any}},
         MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Perform stochastic analytical continuation simulation, parallel version.
The arguments `p1` and `p2` are copies of PBASE and PStochSK, respectively.
"""
function prun(S::StochSKSolver,
              p1::Dict{String,Vector{Any}},
              p2::Dict{String,Vector{Any}},
              MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    rev_dict(p1)
    rev_dict(S, p2)

    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    nstep = get_k("nstep")
    retry = get_k("retry")
    output_per_steps = get_k("ndump")
    measure_per_steps = 10

    println("Start thermalization...")
    SE = warmup(MC, SE, SC)

    shuffle(MC, SE, SC)

    step = 0.0
    for iter = 1:nstep
        if iter % retry == 0
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
        end

        sample(MC, SE, SC)

        if iter % measure_per_steps == 0
            step = step + 1.0
            measure(SE, SC)
        end

        if iter % output_per_steps == 0
            prog = round(I64, iter / nstep * 100)
            println("Start stochastic sampling (prog: $prog)")
            flush(stdout)
            myid() == 2 && write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    average(step::F64, SC::StochSKContext)

Postprocess the results generated during the stochastic analytical
continuation simulations. It will generate real spectral functions.
"""
function average(step::F64, SC::StochSKContext)
    SC.Aout = SC.Aout / (step * (SC.mesh[2] - SC.mesh[1]))

    return SC.Aout
end

"""
    last(SC::StochSKContext, Asum::Vector{F64})

It will process and write the calculated results by the StochSK solver,
including effective hamiltonian, final spectral function, reproduced
correlator.
"""
function last(SC::StochSKContext, Asum::Vector{F64})
    # Calculate final spectral function
    write_spectrum(SC.mesh, Asum)

    # Reproduce input data
    kernel = make_kernel(SC.mesh, SC.grid)
    G = reprod(SC.mesh, kernel, Asum)
    write_backward(SC.grid, G)

    # Calculate full response function at real frequency
    _G = kramers(SC.mesh, Asum)
    write_complete(SC.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

function warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    nwarm = get_k("nwarm")
    ratio = get_k("ratio")

    for i = 1:nwarm
        shuffle(MC, SE, SC)

        push!(SC.𝒞ᵧ, deepcopy(SE))
        SC.χ²vec[i] = SC.χ²
        SC.Θvec[i] = SC.Θ

        @show i, SC.χ², SC.χ²min, SC.χ² - SC.χ²min
        if SC.χ² - SC.χ²min < 1e-3
            break
        end

        SC.Θ = SC.Θ * ratio
    end

    num_anneal = length(SC.𝒞ᵧ)
    @assert num_anneal ≤ nwarm

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
    if rand(MC.rng) > 0.95
        try_move_s(MC, SE, SC)
    else
        try_move_p(MC, SE, SC)
    end
end

function measure(SE::StochSKElement, SC::StochSKContext)
    nmesh = get_b("nmesh")
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    for j = 1:ngamm
        d_pos = SE.P[j]
        s_pos = ceil(I64, d_pos / nfine * nmesh)
        SC.Aout[s_pos] = SC.Aout[s_pos] + SE.A
    end
end

function shuffle(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
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

        sample(MC, SE, SC)

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

#=
### *Service Functions*
=#

"""
    init_mc(S::StochSKSolver)

Try to create a StochSKMC struct.

See also: [`StochSK`](@ref).
"""
function init_mc(S::StochSKSolver)
    seed = rand(1:1000000); seed = 840443
    rng = MersenneTwister(seed)
    Sacc = 0
    Stry = 0
    Pacc = 0
    Ptry = 0

    MC = StochSKMC(rng, Sacc, Stry, Pacc, Ptry)

    return MC
end

"""
    init_element(S::StochSKSolver, rng::AbstractRNG)

Randomize the configurations for future monte carlo sampling. It will
return a StochSKElement object.

See also: [`StochSKElement`](@ref).
"""
function init_element(S::StochSKSolver, rng::AbstractRNG)
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

"""
    init_iodata(S::StochSKSolver, rd::RawData)

Preprocess the input data (`rd`), then allocate memory for the calculated
spectral functions.

See also: [`RawData`](@ref).
"""
function init_iodata(S::StochSKSolver, rd::RawData)
    nmesh = get_b("nmesh")

    Aout = zeros(F64, nmesh)

    G = make_data(rd)
    Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹, Aout
end

"""
    calc_fmesh(S::StochSKSolver)

Try to calculate very fine (dense) linear mesh in [wmin, wmax], which
is used internally to build the kernel function.

See also: [`LinearMesh`](@ref).
"""
function calc_fmesh(S::StochSKSolver)
    nfine = get_k("nfine")
    wmin = get_b("wmin")
    wmax = get_b("wmax")

    fmesh = LinearMesh(nfine, wmin, wmax)

    return fmesh
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

function try_move_s(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
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

function try_move_p(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    MC.Pacc = 0
    MC.Ptry = ngamm
    @assert 1 < SE.W ≤ nfine

    for i = 1:ngamm
        s₁ = rand(MC.rng, 1:ngamm)
        s₂ = s₁
        while s₁ == s₂
            s₂ = rand(MC.rng, 1:ngamm)
        end

        pcurr₁ = SE.P[s₁]
        pcurr₂ = SE.P[s₂]

        if 1 < SE.W < nfine
            move_width_1 = rand(MC.rng, 1:SE.W)
            move_width_2 = rand(MC.rng, 1:SE.W)

            if rand(MC.rng) > 0.5
                pnext₁ = pcurr₁ + move_width_1
                pnext₂ = pcurr₂ - move_width_2
            else
                pnext₁ = pcurr₁ - move_width_1
                pnext₂ = pcurr₂ - move_width_2
            end

            pnext₁ < 1     && (pnext₁ = pnext₁ + nfine)
            pnext₁ > nfine && (pnext₁ = pnext₁ - nfine)
            pnext₂ < 1     && (pnext₂ = pnext₂ + nfine)
            pnext₂ > nfine && (pnext₂ = pnext₂ - nfine)
        else
            pnext₁ = rand(MC.rng, 1:nfine)
            pnext₂ = rand(MC.rng, 1:nfine)
        end

        Knext₁ = view(SC.kernel, :, pnext₁)
        Kcurr₁ = view(SC.kernel, :, pcurr₁)
        Knext₂ = view(SC.kernel, :, pnext₂)
        Kcurr₂ = view(SC.kernel, :, pcurr₂)
        Gₙ = SC.Gᵧ + SE.A * (Knext₁ - Kcurr₁ + Knext₂ - Kcurr₂)
        χ²new = calc_goodness(Gₙ, SC.Gᵥ, SC.σ¹)
        prob = exp( 0.5 * (SC.χ² - χ²new) / SC.Θ )

        if rand(MC.rng) < min(prob, 1.0)
            SE.P[s₁] = pnext₁
            SE.P[s₂] = pnext₂
            SC.Gᵧ = Gₙ

            SC.χ² = χ²new
            if χ²new < SC.χ²min
                SC.χ²min = χ²new
            end

            MC.Pacc = MC.Pacc + 1
        end
    end
end
