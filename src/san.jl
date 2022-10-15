#
# Project : Gardenia
# Source  : san.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/10/15
#

#=
### *Customized Structs* : *StochSK Solver*
=#

"""
    StochSKElement

Mutable struct. It is used to record the field configurations, which will
be sampled by Monte Carlo sweeping procedure.

### Members

* P -> It means the positions of the δ functions.
* A -> It means the weights / amplitudes of the δ functions.
* W -> It denotes the window that is used to shift the δ functions.
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
* Gᵧ     -> Generated correlator.
* σ¹     -> Actually 1.0 / σ¹.
* allow  -> Allowable indices.
* grid   -> Grid for input data.
* mesh   -> Mesh for output spectrum.
* kernel -> Default kernel function.
* Aout   -> Calculated spectral function.
* χ²     -> Current goodness function.
* χ²min  -> Mininum goodness function.
* χ²vec  -> Vector of goodness function.
* Θ      -> Current Θ parameter.
* Θvec   -> Vector of Θ parameter.
* 𝒞ᵧ     -> Historical field configuration.
"""
mutable struct StochSKContext
    Gᵥ     :: Vector{F64}
    Gᵧ     :: Vector{F64}
    σ¹     :: Vector{F64}
    allow  :: Vector{I64}
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

Solve the analytical continuation problem by the stochastic analytical
continuation algorithm (A. W. Sandvik's version).
"""
function solve(S::StochSKSolver, rd::RawData)
    nmesh = get_b("nmesh")

    println("[ StochSK ]")
    MC, SE, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
        println("Using $(nworkers()) workers")
        #
        # Copy configuration dicts
        p1 = deepcopy(PBASE)
        p2 = deepcopy(PStochSK)
        #
        # Launch the tasks one by one
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

    # Sequential version
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
    allow = constraints(S)

    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    SE = init_element(S, MC.rng, allow)
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
    println("Precompute correlator")

    χ = calc_goodness(Gᵧ, Gᵥ, σ¹)
    χ², χ²min = χ, χ
    χ²vec = zeros(F64, get_k("nwarm"))
    println("Precompute goodness function")

    Θ = get_k("theta")
    Θvec = zeros(F64, get_k("nwarm"))
    println("Setup Θ parameter")

    𝒞ᵧ = StochSKElement[]
    println("Setup historical Monte Carlo configurations")

    SC = StochSKContext(Gᵥ, Gᵧ, σ¹, allow, grid, mesh, kernel, Aout,
                        χ², χ²min, χ²vec, Θ, Θvec, 𝒞ᵧ)

    return MC, SE, SC
end

"""
    run(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Perform stochastic analytical continuation simulation, sequential version.
"""
function run(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Setup essential parameters
    nstep = get_k("nstep")
    retry = get_k("retry")
    output_per_steps = get_k("ndump")
    measure_per_steps = 10

    # Warmup the Monte Carlo engine 
    println("Start thermalization...")
    warmup(MC, SE, SC)

    # Shuffle the Monte Carlo configuration again
    shuffle(MC, SE, SC)

    # Sample and collect data
    step = 0.0
    println("Start stochastic sampling...")
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
            @printf("step = %6i ", iter)
            @printf("(progress = %3i)\n", prog)
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
    # Revise parameteric dicts
    rev_dict(p1)
    rev_dict(S, p2)

    # Initialize random number generator again
    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    # Setup essential parameters
    nstep = get_k("nstep")
    retry = get_k("retry")
    output_per_steps = get_k("ndump")
    measure_per_steps = 10

    # Warmup the Monte Carlo engine 
    println("Start thermalization...")
    warmup(MC, SE, SC)

    # Shuffle the Monte Carlo configuration again
    shuffle(MC, SE, SC)

    # Sample and collect data
    step = 0.0
    println("Start stochastic sampling...")
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
            @printf("step = %6i ", iter)
            @printf("(progress = %3i)\n", prog)
            flush(stdout)
            myid() == 2 && write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    average(step::F64, SC::StochSKContext)

Postprocess the results generated during the stochastic analytical
continuation simulations. It will generate the spectral functions.
"""
function average(step::F64, SC::StochSKContext)
    SC.Aout = SC.Aout ./ (step * SC.mesh.weight)

    return SC.Aout
end

"""
    last(SC::StochSKContext, Asum::Vector{F64})

It will process and write the calculated results by the StochSK solver,
including final spectral function and reproduced correlator.
"""
function last(SC::StochSKContext, Asum::Vector{F64})
    # Write final spectral function
    write_spectrum(SC.mesh, Asum)

    # Write Θ-dependent goodness function
    write_goodness(SC.Θvec, SC.χ²vec)

    # Reproduce input data and write them
    kernel = make_kernel(SC.mesh, SC.grid)
    G = reprod(SC.mesh, kernel, Asum)
    write_backward(SC.grid, G)

    # Calculate full response function on real axis and write them
    _G = kramers(SC.mesh, Asum)
    write_complete(SC.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

"""
    warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Warmup the Monte Carlo engine to acheieve thermalized equilibrium. Then
it will try to figure out the optimized Θ and the corresponding Monte
Carlo field configuration.
"""
function warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get essential parameters
    nwarm = get_k("nwarm")
    ratio = get_k("ratio")

    # Change the Θ parameter and approch the equilibrium state
    for i = 1:nwarm
        # Shuffle the Monte Carlo configurations
        shuffle(MC, SE, SC)

        # Backup key parameters and Monte Carlo field configurations
        SC.χ²vec[i] = SC.χ²
        SC.Θvec[i] = SC.Θ
        push!(SC.𝒞ᵧ, deepcopy(SE))

        # Check whether the equilibrium state is reached 
        δχ² = SC.χ² - SC.χ²min
        @printf("step : %5i ", i)
        @printf("χ² - χ²min -> %12.6e\n", δχ²)
        if δχ² < 1e-3
            println("Reach equilibrium state")
            break
        end

        # Adjust the Θ parameter
        SC.Θ = SC.Θ * ratio
    end

    # Well, we have vectors for Θ and χ². We have to figure out the
    # optimized Θ and χ², and then extract the corresponding Monte
    # Carlo field configuration.
    c = length(SC.𝒞ᵧ)
    while c ≥ 1
        if SC.χ²vec[c] > SC.χ²min + 2.0 * sqrt(SC.χ²min)
            break
        end
        c = c - 1
    end
    @assert 1 ≤ c ≤ length(SC.𝒞ᵧ)

    # Retrieve the Monte Carlo field configuration
    @. SE.P = SC.𝒞ᵧ[c].P
    SE.A = SC.𝒞ᵧ[c].A
    SE.W = SC.𝒞ᵧ[c].W

    # Reset Θ
    SC.Θ = SC.Θvec[c]

    # Update Gᵧ and χ²
    SC.Gᵧ = calc_correlator(SE, SC.kernel)
    SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
    println("Θ = ", SC.Θ, " χ² = ", SC.χ², " (step = $c)")
end

"""
    sample(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Perform Monte Carlo sweeps and sample the field configurations.
"""
function sample(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    if rand(MC.rng) > 0.95
        try_move_s(MC, SE, SC)
    else
        try_move_p(MC, SE, SC)
    end
end

"""
    measure(SE::StochSKElement, SC::StochSKContext)

Measure the final spectral functions.

See also: [`nearest`](@ref).
"""
function measure(SE::StochSKElement, SC::StochSKContext)
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    for j = 1:ngamm
        d_pos = SE.P[j]
        # d_pos / nfine denotes the position of the selected δ function
        # in the fine linear mesh.
        #
        # The nearest() function is used to extract the approximated
        # position (index) of the selected δ function in the spectral
        # mesh, which could be linear or non-linear.
        s_pos = nearest(SC.mesh, d_pos / nfine)
        SC.Aout[s_pos] = SC.Aout[s_pos] + SE.A
    end
end

"""
    shuffle(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Try to shuffle the Monte Carlo field configuration via the Metropolis
algorithm. Then the window for shifting the δ functions is adjusted.
"""
function shuffle(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get/set essential parameters
    nfine = get_k("nfine")
    retry = get_k("retry")
    max_bin_size = 100 # You can increase it to improve the accuracy

    # Allocate memory
    bin_χ²  = zeros(F64, max_bin_size)
    bin_acc = zeros(I64, max_bin_size)
    bin_try = zeros(I64, max_bin_size)

    # Perform Monte Carlo sweeping
    for s = 1:max_bin_size
        # Recalculate the goodness function
        if s % retry == 0
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ, SC.σ¹)
        end

        sample(MC, SE, SC)

        # Update the counters
        bin_χ²[s]  = SC.χ²
        bin_acc[s] = MC.Sacc + MC.Pacc
        bin_try[s] = MC.Stry + MC.Ptry
    end

    # Calculate the transition probability, and then adjust the window,
    # which restricts the movement of the δ functions.
    #
    # The transition probability will be kept around 0.5.
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

    # Update χ² with averaged χ² 
    SC.χ² = mean(bin_χ²)
end

#=
### *Service Functions*
=#

"""
    init_mc(S::StochSKSolver)

Try to create a StochSKMC struct.

See also: [`StochSKMC`](@ref).
"""
function init_mc(S::StochSKSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    Sacc = 0
    Stry = 0
    Pacc = 0
    Ptry = 0

    MC = StochSKMC(rng, Sacc, Stry, Pacc, Ptry)

    return MC
end

"""
    init_element(S::StochSKSolver, rng::AbstractRNG, allow::Vector{I64})

Randomize the configurations for future Monte Carlo sampling. It will
return a StochSKElement object.

See also: [`StochSKElement`](@ref).
"""
function init_element(S::StochSKSolver, rng::AbstractRNG, allow::Vector{I64})
    β = get_b("beta")
    wmax = get_b("wmax")
    wmin = get_b("wmin")
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    position = rand(rng, allow, ngamm)
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

"""
    calc_correlator(SE::StochSKElement, kernel::Array{F64,2})

Try to calculate correlator with the kernel function and the Monte Carlo
field configuration. This correlator will then be used to evaluate the
goodness function.

See also: [`calc_goodness`](@ref).
"""
function calc_correlator(SE::StochSKElement, kernel::Array{F64,2})
    ngamm = length(SE.P)
    𝐴 = fill(SE.A, ngamm)
    𝐾 = kernel[:, SE.P]
    return 𝐾 * 𝐴
end

"""
    calc_goodness(Gₙ::Vector{F64,}, Gᵥ::Vector{F64}, σ¹::Vector{F64})

Try to calculate the goodness function (i.e, χ²), which measures the
distance between input and regenerated correlators.

See also: [`calc_correlator`](@ref).
"""
function calc_goodness(Gₙ::Vector{F64,}, Gᵥ::Vector{F64}, σ¹::Vector{F64})
    χ = sum( ( (Gₙ .- Gᵥ) .* σ¹ ) .^ 2.0 )
    return χ
end

"""
    constraints(S::StochSKSolver)

Try to implement the constrained stochastic analytical continuation
method. This function will return a collection. It contains all the
allowable indices.

See also: [`StochSKSolver`](@ref).
"""
function constraints(S::StochSKSolver)
    exclude = get_b("exclude")
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nfine = get_k("nfine")

    allow = I64[]
    mesh = collect(LinRange(wmin, wmax, nfine))

    # Go through the fine linear mesh and check each mesh point.
    # Is is excluded ?
    for i in eachindex(mesh)
        is_excluded = false
        #
        if !isa(exclude, Missing)
            for j in eachindex(exclude)
                if exclude[j][1] ≤ mesh[i] ≤ exclude[j][2]
                    is_excluded = true
                    break
                end
            end
        end
        #
        if !is_excluded
            push!(allow, i)
        end
    end

    return allow
end

"""
    try_move_s(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Try to update the Monte Carlo field configurations via the Metropolis
algorithm. In each update, only single δ function is shifted.

See also: [`try_move_p`](@ref).
"""
function try_move_s(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get parameters
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    # Reset counters
    MC.Sacc = 0
    MC.Stry = ngamm
    @assert 1 < SE.W ≤ nfine

    # Allocate memory for new correlator
    Gₙ = zeros(F64, size(SC.Gᵧ))
    ΔG = zeros(F64, size(SC.Gᵧ))

    for i = 1:ngamm
        # Choose single δ function
        s = rand(MC.rng, 1:ngamm)

        # Evaluate new position for the δ function
        pcurr = SE.P[s]
        #
        if 1 < SE.W < nfine
            δW = rand(MC.rng, 1:SE.W)
            #
            if rand(MC.rng) > 0.5
                pnext = pcurr + δW
            else
                pnext = pcurr - δW
            end
            #
            pnext < 1     && (pnext = pnext + nfine)
            pnext > nfine && (pnext = pnext - nfine)
        else
            pnext = rand(MC.rng, 1:nfine)
        end

        # Apply the constraints
        !(pnext in SC.allow) && continue

        # Calculate the transition probability
        Knext = view(SC.kernel, :, pnext)
        Kcurr = view(SC.kernel, :, pcurr)
        #
        @. Gₙ = SC.Gᵧ + SE.A * (Knext - Kcurr)
        @. ΔG = (Gₙ - SC.Gᵥ) * SC.σ¹
        χ²new = dot(ΔG, ΔG)
        #
        prob = exp( 0.5 * (SC.χ² - χ²new) / SC.Θ )

        # Important sampling, if true, the δ function is shifted and the
        # corresponding objects are updated.
        if rand(MC.rng) < min(prob, 1.0)
            SE.P[s] = pnext
            @. SC.Gᵧ = Gₙ
            #
            SC.χ² = χ²new
            if χ²new < SC.χ²min
                SC.χ²min = χ²new
            end
            #
            MC.Sacc = MC.Sacc + 1
        end
    end
end

"""
    try_move_p(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Try to update the Monte Carlo field configurations via the Metropolis
algorithm. In each update, only a pair of δ functions is shifted.

See also: [`try_move_s`](@ref).
"""
function try_move_p(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get parameters
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    # Reset counters
    MC.Pacc = 0
    MC.Ptry = ngamm
    @assert 1 < SE.W ≤ nfine

    # Allocate memory for new correlator
    Gₙ = zeros(F64, size(SC.Gᵧ))
    ΔG = zeros(F64, size(SC.Gᵧ))

    for i = 1:ngamm
        # Choose a pair of δ functions
        s₁ = rand(MC.rng, 1:ngamm)
        s₂ = s₁
        while s₁ == s₂
            s₂ = rand(MC.rng, 1:ngamm)
        end

        # Evaluate new positions for the two δ functions
        pcurr₁ = SE.P[s₁]
        pcurr₂ = SE.P[s₂]
        #
        if 1 < SE.W < nfine
            δW₁ = rand(MC.rng, 1:SE.W)
            δW₂ = rand(MC.rng, 1:SE.W)
            #
            if rand(MC.rng) > 0.5
                pnext₁ = pcurr₁ + δW₁
                pnext₂ = pcurr₂ - δW₂
            else
                pnext₁ = pcurr₁ - δW₁
                pnext₂ = pcurr₂ + δW₂
            end
            #
            pnext₁ < 1     && (pnext₁ = pnext₁ + nfine)
            pnext₁ > nfine && (pnext₁ = pnext₁ - nfine)
            pnext₂ < 1     && (pnext₂ = pnext₂ + nfine)
            pnext₂ > nfine && (pnext₂ = pnext₂ - nfine)
        else
            pnext₁ = rand(MC.rng, 1:nfine)
            pnext₂ = rand(MC.rng, 1:nfine)
        end

        # Apply the constraints
        !(pnext₁ in SC.allow) && continue
        !(pnext₂ in SC.allow) && continue

        # Calculate the transition probability
        Knext₁ = view(SC.kernel, :, pnext₁)
        Kcurr₁ = view(SC.kernel, :, pcurr₁)
        Knext₂ = view(SC.kernel, :, pnext₂)
        Kcurr₂ = view(SC.kernel, :, pcurr₂)
        #
        @. Gₙ = SC.Gᵧ + SE.A * (Knext₁ - Kcurr₁ + Knext₂ - Kcurr₂)
        @. ΔG = (Gₙ - SC.Gᵥ) * SC.σ¹
        χ²new = dot(ΔG, ΔG)
        #
        prob = exp( 0.5 * (SC.χ² - χ²new) / SC.Θ )

        # Important sampling, if true, the δ function is shifted and the
        # corresponding objects are updated.
        if rand(MC.rng) < min(prob, 1.0)
            SE.P[s₁] = pnext₁
            SE.P[s₂] = pnext₂
            @. SC.Gᵧ = Gₙ
            #
            SC.χ² = χ²new
            if χ²new < SC.χ²min
                SC.χ²min = χ²new
            end
            #
            MC.Pacc = MC.Pacc + 1
        end
    end
end
