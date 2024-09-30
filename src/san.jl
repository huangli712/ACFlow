#
# Project : Gardenia
# Source  : san.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/09/30
#

#=
### *Customized Structs* : *StochSK Solver*
=#

"""
    StochSKElement

Mutable struct. It is used to record the field configurations, which will
be sampled by Monte Carlo sweeping procedure.

In the present implementation of StochSK solver, the amplitudes of the δ
functions are fixed. But in principles, they could be sampled in the Monte
Carlo procedure.

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
* grid   -> Imaginary axis grid for input data.
* mesh   -> Real frequency mesh for output spectrum.
* kernel -> Default kernel function.
* Aout   -> Calculated spectral function.
* χ²     -> Current goodness-of-fit function.
* χ²min  -> Mininum goodness-of-fit function.
* χ²vec  -> Vector of goodness-of-fit function.
* Θ      -> Current Θ parameter.
* Θvec   -> Vector of Θ parameter.
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
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochSKSolver, rd::RawData)

Solve the analytic continuation problem by using the stochastic analytic
continuation algorithm (A. W. Sandvik's version). It is the driver for the
StochSK solver.

If the input correlators are bosonic, this solver will return A(ω) / ω
via `Aout`, instead of A(ω). At this time, `Aout` is not compatible with
`Gout`. If the input correlators are fermionic, this solver will return
A(ω) in `Aout`. Now it is compatible with `Gout`. These behaviors are just
similar to the MaxEnt, StochAC, and StochOM solvers.

Now the StochSK solver supports both continuous and δ-like spectra.

### Arguments
* S -> A StochSKSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* mesh -> Real frequency mesh, ω.
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
"""
function solve(S::StochSKSolver, rd::RawData)
    nmesh = get_b("nmesh")
    nwarm = get_k("nwarm")

    println("[ StochSK ]")
    MC, SE, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
        #
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
        χ²out = zeros(F64, nwarm)
        Θout = zeros(F64, nwarm)
        for i in eachindex(sol)
            a, b, c = sol[i]
            @. Aout = Aout + a / nworkers()
            @. χ²out = χ²out + b / nworkers()
            @. Θout = Θout + c / nworkers()
        end
        #
        # Postprocess the solutions
        Gout = last(SC, Aout, χ²out, Θout)
        #
    # Sequential version
    else
        #
        Aout, χ²out, Θout = run(MC, SE, SC)
        Gout = last(SC, Aout, χ²out, Θout)
        #
    end

    return SC.mesh.mesh, Aout, Gout
end

"""
    init(S::StochSKSolver, rd::RawData)

Initialize the StochSK solver and return the StochSKMC, StochSKElement,
and StochSKContext structs. Please don't call this function directly.

### Arguments
* S -> A StochSKSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.
"""
function init(S::StochSKSolver, rd::RawData)
    # Initialize possible constraints.
    # The array allow contains all the possible indices for δ functions.
    fmesh = calc_fmesh(S)
    allow = constraints(S, fmesh)

    # Prepare input data
    Gᵥ, σ¹ = init_iodata(S, rd)
    println("Postprocess input data: ", length(σ¹), " points")

    # Prepare grid for input data
    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    # Prepare mesh for output spectrum
    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    # Initialize counters for Monte Carlo engine
    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    # Initialize Monte Carlo configurations
    SE = init_element(S, MC.rng, allow)
    println("Randomize Monte Carlo configurations")

    # Prepare some key variables
    SC = init_context(SE, Gᵥ, σ¹, allow, grid, mesh, fmesh)
    println("Initialize context for the StochSK solver")

    return MC, SE, SC
end

"""
    run(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Perform stochastic analytic continuation simulation, sequential version.

### Arguments
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* χ²vec -> Θ-dependent χ², χ²(Θ).
* Θvec -> List of Θ parameters.
"""
function run(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

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
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ)
        end

        sample(MC, SE, SC)

        if iter % measure_per_steps == 0
            step = step + 1.0
            measure(SE, SC)
        end

        if iter % output_per_steps == 0
            prog = round(I64, iter / nstep * 100)
            @printf("step = %6i ", iter)
            @printf("[progress = %3i]\n", prog)
            flush(stdout)
            fwrite && write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    prun(
        S::StochSKSolver,
        p1::Dict{String,Vector{Any}},
        p2::Dict{String,Vector{Any}},
        MC::StochSKMC,
        SE::StochSKElement,
        SC::StochSKContext
    )

Perform stochastic analytic continuation simulation, parallel version.
The arguments `p1` and `p2` are copies of PBASE and PStochSK, respectively.

### Arguments
* S -> A StochSKSolver struct.
* p1 -> A copy of PBASE.
* p2 -> A copy of PStochSK.
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* χ²vec -> Θ-dependent χ², χ²(Θ).
* Θvec -> List of Θ parameters.
"""
function prun(
    S::StochSKSolver,
    p1::Dict{String,Vector{Any}},
    p2::Dict{String,Vector{Any}},
    MC::StochSKMC,
    SE::StochSKElement,
    SC::StochSKContext
    )
    # Revise parameteric dicts
    rev_dict_b(p1)
    rev_dict_k(S, p2)

    # Initialize random number generator again
    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

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
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ)
        end

        sample(MC, SE, SC)

        if iter % measure_per_steps == 0
            step = step + 1.0
            measure(SE, SC)
        end

        if iter % output_per_steps == 0
            prog = round(I64, iter / nstep * 100)
            @printf("step = %6i ", iter)
            @printf("[progress = %3i]\n", prog)
            flush(stdout)
            myid() == 2 && fwrite && write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    average(step::F64, SC::StochSKContext)

Postprocess the results generated during the stochastic analytic
continuation simulations. It will generate the spectral functions.

### Arguments
* step -> Number of Monte Carlo samplings.
* SC   -> A StochSKContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* χ²vec -> Θ-dependent χ², χ²(Θ).
* Θvec -> List of Θ parameters.
"""
function average(step::F64, SC::StochSKContext)
    #
    # Here, the factor SC.mesh.weight in denominator is used to make sure
    # that the sum-rule
    #
    # ∫ A(ω) dω = 1
    #
    # is obeyed by the obtained spectral function A(ω).
    #
    SC.Aout = SC.Aout ./ (step * SC.mesh.weight)

    return SC.Aout, SC.χ²vec, SC.Θvec
end

"""
    last(
        SC::StochSKContext,
        Asum::Vector{F64},
        χ²vec::Vector{F64},
        Θvec::Vector{F64}
    )

It will process and write the calculated results by the StochSK solver,
including the final spectral function and reproduced correlator.

### Arguments
* SC    -> A StochSKContext struct.
* Asum  -> Spectral function, A(ω).
* χ²vec -> Θ-dependent χ².
* Θvec  -> List of Θ parameters.

### Returns
* G -> Retarded Green's function, G(ω).
"""
function last(
    SC::StochSKContext,
    Asum::Vector{F64},
    χ²vec::Vector{F64},
    Θvec::Vector{F64}
    )
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Write Θ-dependent goodness-of-fit function
    fwrite && write_goodness(Θvec, χ²vec)

    # Write final spectral function
    fwrite && write_spectrum(SC.mesh, Asum)

    # Reproduce input data and write them
    kernel = make_kernel(SC.mesh, SC.grid)
    G = reprod(SC.mesh, kernel, Asum)
    fwrite && write_backward(SC.grid, G)

    # Calculate full response function on real axis and write them
    if get_b("ktype") == "fermi"
        _G = kramers(SC.mesh, Asum)
    else
        _G = kramers(SC.mesh, Asum .* SC.mesh)
    end
    fwrite && write_complete(SC.mesh, _G)

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

### Arguments
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
N/A
"""
function warmup(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get essential parameters
    nwarm = get_k("nwarm")
    ratio = get_k("ratio")
    threshold = 1e-3

    # To store the historic Monte Carlo field configurations
    𝒞ᵧ = StochSKElement[]

    # Change the Θ parameter and approch the equilibrium state
    for i = 1:nwarm
        # Shuffle the Monte Carlo configurations
        shuffle(MC, SE, SC)

        # Backup key parameters and Monte Carlo field configurations
        SC.χ²vec[i] = SC.χ²
        SC.Θvec[i] = SC.Θ
        push!(𝒞ᵧ, deepcopy(SE))

        # Check whether the equilibrium state is reached
        δχ² = SC.χ² - SC.χ²min
        @printf("step : %5i ", i)
        @printf("χ² - χ²min -> %12.6e\n", δχ²)
        if δχ² < threshold
            println("Reach equilibrium state")
            break
        else
            if i == nwarm
                error("Fail to reach equilibrium state")
            end
        end

        # Adjust the Θ parameter
        SC.Θ = SC.Θ * ratio
    end

    # Well, we have vectors for Θ and χ². We have to figure out the
    # optimized Θ and χ², and then extract the corresponding Monte
    # Carlo field configuration.
    c = calc_theta(length(𝒞ᵧ), SC)
    @assert 1 ≤ c ≤ length(𝒞ᵧ)

    # Retrieve the Monte Carlo field configuration
    @. SE.P = 𝒞ᵧ[c].P
    SE.A = 𝒞ᵧ[c].A
    SE.W = 𝒞ᵧ[c].W

    # Reset Θ
    SC.Θ = SC.Θvec[c]

    # Update Gᵧ and χ²
    SC.Gᵧ = calc_correlator(SE, SC.kernel)
    SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ)
    println("Θ = ", SC.Θ, " χ² = ", SC.χ², " (step = $c)")
end

"""
    sample(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Perform Monte Carlo sweeps and sample the field configurations.

### Arguments
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
N/A
"""
function sample(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    if rand(MC.rng) < 0.80
        try_move_s(MC, SE, SC)
    else
        if rand(MC.rng) < 0.50
            try_move_p(MC, SE, SC)
        else
            try_move_q(MC, SE, SC)
        end
    end
end

"""
    measure(SE::StochSKElement, SC::StochSKContext)

Accumulate the final spectral functions A(ω).

### Arguments
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
N/A

See also: [`nearest`](@ref).
"""
function measure(SE::StochSKElement, SC::StochSKContext)
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    for j = 1:ngamm
        d_pos = SE.P[j]
        # d_pos / nfine denotes the position of the selected δ-like peak
        # in the fine linear mesh.
        #
        # The nearest() function is used to extract the approximated
        # position (index) of the selected δ function in the spectral
        # mesh, which could be linear or non-linear.
        #
        # Note that nearest() is defined in mesh.jl.
        s_pos = nearest(SC.mesh, d_pos / nfine)
        SC.Aout[s_pos] = SC.Aout[s_pos] + SE.A
    end
end

"""
    shuffle(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Try to shuffle the Monte Carlo field configuration via the Metropolis
algorithm. Then the window for shifting the δ functions is adjusted.

### Arguments
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
N/A
"""
function shuffle(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get/set essential parameters
    nfine = get_k("nfine")
    retry = get_k("retry")
    max_bin_size = 100 # You can increase it to improve the accuracy

    # Announce counters
    bin_χ²  = 0.0
    bin_acc = 0.0
    bin_try = 0.0

    # Perform Monte Carlo sweeping
    for s = 1:max_bin_size
        # Recalculate the goodness-of-fit function
        if s % retry == 0
            SC.χ² = calc_goodness(SC.Gᵧ, SC.Gᵥ)
        end

        sample(MC, SE, SC)

        # Update the counters
        bin_χ² = bin_χ² + SC.χ²
        bin_acc = bin_acc + (MC.Sacc + MC.Pacc)
        bin_try = bin_try + (MC.Stry + MC.Ptry)
    end

    # Calculate the transition probability, and then adjust the window,
    # which restricts the movement of the δ functions.
    #
    # The transition probability will be kept around 0.5.
    𝑝 = bin_acc / bin_try
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
    SC.χ² = bin_χ² / max_bin_size
end

#=
### *Service Functions*
=#

"""
    init_iodata(S::StochSKSolver, rd::RawData)

Preprocess the input data (`rd`).

### Arguments
* S -> A StochSKSolver struct.
* rd -> A RawData struct, which contains essential input data.

### Returns
* Gᵥ -> Input correlator.
* σ¹ -> 1.0 / σ¹.

See also: [`RawData`](@ref).
"""
function init_iodata(S::StochSKSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value # Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹
end

"""
    init_mc(S::StochSKSolver)

Try to create a StochSKMC struct. Some counters for Monte Carlo updates
are initialized here.

### Arguments
* S -> A StochSKSolver struct.

### Returns
* MC -> A StochSKMC struct.

See also: [`StochSKMC`](@ref).
"""
function init_mc(S::StochSKSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    #
    Sacc = 0
    Stry = 0
    Pacc = 0
    Ptry = 0
    Qacc = 0
    Qtry = 0
    #
    MC = StochSKMC(rng, Sacc, Stry, Pacc, Ptry, Qacc, Qtry)

    return MC
end

"""
    init_element(
        S::StochSKSolver,
        rng::AbstractRNG,
        allow::Vector{I64}
    )

Randomize the configurations for future Monte Carlo sampling. It will
return a StochSKElement struct.

### Arguments
* S     -> A StochSKSolver struct.
* rng   -> Random number generator.
* allow -> Allowed positions for the δ peaks.

### Returns
* SE -> A StochSKElement struct.

See also: [`StochSKElement`](@ref).
"""
function init_element(
    S::StochSKSolver,
    rng::AbstractRNG,
    allow::Vector{I64}
    )
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
    init_context(
        SE::StochSKElement,
        Gᵥ::Vector{F64},
        σ¹::Vector{F64},
        allow::Vector{I64},
        grid::AbstractGrid,
        mesh::AbstractMesh,
        fmesh::AbstractMesh
    )

Try to create a StochSKContext struct, which contains some key variables,
including grid, mesh, input correlator and the corresponding standard
deviation, kernel matrix, spectral function, and goodness-of-fit function.

### Arguments
* SE -> A StochSKElement struct.
* Gᵥ -> Input correlator. It will be changed in this function.
* σ¹ -> Standard deviation for input correlator.
* allow -> Allowable indices for δ-like peaks.
* grid -> Imaginary axis grid for input data.
* mesh -> Real frequency mesh for output spectrum.
* fmesh -> Very fine mesh in [wmin, wmax].

### Returns
* SC -> A StochSKContext struct.
"""
function init_context(
    SE::StochSKElement,
    Gᵥ::Vector{F64},
    σ¹::Vector{F64},
    allow::Vector{I64},
    grid::AbstractGrid,
    mesh::AbstractMesh,
    fmesh::AbstractMesh
    )
    # Get parameters
    nmesh = get_b("nmesh")
    nwarm = get_k("nwarm")
    Θ = get_k("theta")

    # Allocate memory for spectral function, A(ω)
    Aout = zeros(F64, nmesh)

    # Allocate memory for χ² and Θ
    χ²vec = zeros(F64, nwarm)
    Θvec = zeros(F64, nwarm)

    # Build kernel matrix
    kernel = make_kernel(fmesh, grid)

    # In order to accelerate the calculations, the singular space of the
    # kernel function is used. At first, we preform singular value
    # decomposition for K/σ:
    #     K/σ = U S Vᵀ
    # Then
    #     (G - KA)/σ = G/σ - K/σA
    #                = UU'(G/σ - USVᵀA)
    #                = U(U'G/σ - U'USVᵀA)
    #                = U(U'G/σ - SVᵀA)
    #                = U(G' - K'A)
    # In the StochAC solver, let Gᵥ → G', kernel → K'. Then new χ² is
    # calculated by
    #     |G' - K'A|²
    # instead of
    #     |G - KA|²/σ²

    # Singular value decomposition of K/σ
    U, V, S = make_singular_space(Diagonal(σ¹) * kernel)

    # Get new kernel matrix
    kernel = Diagonal(S) * V'

    # Get new (input) correlator
    Gᵥ = U' *  (Gᵥ .* σ¹)

    # Calculate reconstructed correlator using current field configuration
    Gᵧ = calc_correlator(SE, kernel)

    # Calculate goodness-of-fit functional χ²
    𝚾 = calc_goodness(Gᵧ, Gᵥ)
    χ², χ²min = 𝚾, 𝚾

    return StochSKContext(Gᵥ, Gᵧ, σ¹, allow, grid, mesh, kernel, Aout,
                        χ², χ²min, χ²vec, Θ, Θvec)
end

"""
    calc_fmesh(S::StochSKSolver)

Try to calculate very fine (dense) linear mesh in [wmin, wmax], which
is used internally to build the kernel function. Note that the stochastic
analytic continuation method (A. W. Sandvik's version) does not support
the self-adaptive mesh.

### Arguments
* S -> A StochSKSolver struct.

### Returns
* fmesh -> A very fine, uniform mesh in [wmin, wmax].

See also: [`LinearMesh`](@ref), [`DynamicMesh`](@ref).
"""
function calc_fmesh(S::StochSKSolver)
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nfine = get_k("nfine")

    fmesh = LinearMesh(nfine, wmin, wmax)

    return fmesh
end

"""
    calc_correlator(SE::StochSKElement, kernel::Array{F64,2})

Try to calculate correlator with the kernel function and the Monte Carlo
field configuration. This correlator will then be used to evaluate the
goodness-of-fit function χ².

### Arguments
* SE     -> A StochSKElement struct.
* kernel -> The fermionic or bosonic kernel.

### Returns
* G -> Reconstructed correlator.

See also: [`calc_goodness`](@ref).
"""
function calc_correlator(SE::StochSKElement, kernel::Array{F64,2})
    ngamm = length(SE.P)
    𝐴 = fill(SE.A, ngamm)
    𝐾 = kernel[:, SE.P]
    return 𝐾 * 𝐴
end

"""
    calc_goodness(Gₙ::Vector{F64}, Gᵥ::Vector{F64})

Try to calculate the goodness-of-fit function (i.e, χ²), which measures
the distance between input and regenerated correlators.

### Arguments
* Gₙ -> Reconstructed correlators.
* Gᵥ -> Input (original) correlators.

### Returns
* χ² -> Goodness-of-fit function.

See also: [`calc_correlator`](@ref).
"""
function calc_goodness(Gₙ::Vector{F64}, Gᵥ::Vector{F64})
    ΔG = Gₙ - Gᵥ
    return dot(ΔG, ΔG)
end

"""
    calc_theta(len::I64, SC::StochSKContext)

Try to locate the optimal Θ and χ². This function implements the `chi2min`
and `chi2kink` algorithms. Note that the `chi2min` algorithm is preferred.

### Arguments
* len -> Length of vector Θ.
* SC -> A StochSKContext struct.

### Returns
* c -> Selected index for optimal Θ.
"""
function calc_theta(len::I64, SC::StochSKContext)
    function fitfun(x, p)
        return @. p[1] + p[2] / (1.0 + exp(-p[4] * (x - p[3])))
    end

    # Which algorithm is preferred ?
    method = get_k("method")

    # Get length of Θ and χ² vectors
    c = len

    # `chi2min` algorithm, proposed by Shao and Sandvik
    if method == "chi2min"
        while c ≥ 1
            if SC.χ²vec[c] > SC.χ²min + 2.0 * sqrt(SC.χ²min)
                break
            end
            c = c - 1
        end
    end

    # `chi2kink` algorithm, inspired by the `chi2kink` algorithm
    # used in MaxEnt solver
    if method == "chi2kink"
        guess = [0.0, 5.0, 2.0, 0.0]
        fit = curve_fit(fitfun, log10.(SC.Θvec[1:c]), log10.(SC.χ²vec[1:c]), guess)
        _, _, a, b = fit.param
        #
        fit_pos = 2.5
        Θ_opt = a - fit_pos / b
        c = argmin( abs.( log10.(SC.Θvec[1:c]) .- Θ_opt ) )
    end

    return c
end

"""
    constraints(S::StochSKSolver, fmesh::AbstractMesh)

Try to implement the constrained stochastic analytic continuation
method. This function will return a collection. It contains all the
allowable indices. Be careful, `fmesh` should be a fine linear mesh.

### Arguments
* S     -> A StochSKSolver struct.
* fmesh -> Very dense mesh for the δ peaks.

### Returns
* allow -> Allowable indices.

See also: [`StochSKSolver`](@ref).
"""
function constraints(S::StochSKSolver, fmesh::AbstractMesh)
    exclude = get_b("exclude")
    nfine = get_k("nfine")
    @assert nfine == length(fmesh)

    allow = I64[]

    # Go through the fine mesh and check every mesh point.
    # Is is excluded?
    for i in eachindex(fmesh)
        is_excluded = false
        #
        if !isa(exclude, Missing)
            for j in eachindex(exclude)
                if exclude[j][1] ≤ fmesh[i] ≤ exclude[j][2]
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

### Arguments
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
N/A

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

    for _ = 1:ngamm
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
        @. ΔG = Gₙ - SC.Gᵥ
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
algorithm. In each update, only a pair of δ functions are shifted.

### Arguments
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
N/A

See also: [`try_move_s`](@ref).
"""
function try_move_p(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get parameters
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    # We have to make sure that there are at least two δ functions here.
    ngamm < 2 && return

    # Reset counters
    MC.Pacc = 0
    MC.Ptry = ngamm
    @assert 1 < SE.W ≤ nfine

    # Allocate memory for new correlator
    Gₙ = zeros(F64, size(SC.Gᵧ))
    ΔG = zeros(F64, size(SC.Gᵧ))

    for _ = 1:ngamm
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
        @. ΔG = Gₙ - SC.Gᵥ
        χ²new = dot(ΔG, ΔG)
        #
        prob = exp( 0.5 * (SC.χ² - χ²new) / SC.Θ )

        # Important sampling, if true, the δ functions are shifted and the
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

"""
    try_move_q(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)

Try to update the Monte Carlo field configurations via the Metropolis
algorithm. In each update, four different δ functions are shifted.

### Arguments
* MC -> A StochSKMC struct.
* SE -> A StochSKElement struct.
* SC -> A StochSKContext struct.

### Returns
N/A

See also: [`try_move_s`](@ref).
"""
function try_move_q(MC::StochSKMC, SE::StochSKElement, SC::StochSKContext)
    # Get parameters
    nfine = get_k("nfine")
    ngamm = get_k("ngamm")

    # We have to make sure that there are at least four δ functions here.
    ngamm < 4 && return

    # Reset counters
    MC.Qacc = 0
    MC.Qtry = ngamm
    @assert 1 < SE.W ≤ nfine

    # Allocate memory for new correlator
    Gₙ = zeros(F64, size(SC.Gᵧ))
    ΔG = zeros(F64, size(SC.Gᵧ))

    for _ = 1:ngamm
        # Choose four different δ functions
        𝑆 = nothing
        while true
            𝑆 = rand(MC.rng, 1:ngamm, 4)
            𝒮 = unique(𝑆)
            if length(𝑆) == length(𝒮)
                break
            end
        end
        s₁, s₂, s₃, s₄ = 𝑆

        # Evaluate new positions for the four δ functions
        pcurr₁ = SE.P[s₁]
        pcurr₂ = SE.P[s₂]
        pcurr₃ = SE.P[s₃]
        pcurr₄ = SE.P[s₄]
        #
        if 1 < SE.W < nfine
            δW₁ = rand(MC.rng, 1:SE.W)
            δW₂ = rand(MC.rng, 1:SE.W)
            δW₃ = rand(MC.rng, 1:SE.W)
            δW₄ = rand(MC.rng, 1:SE.W)
            #
            if rand(MC.rng) > 0.5
                pnext₁ = pcurr₁ + δW₁
                pnext₂ = pcurr₂ - δW₂
                pnext₃ = pcurr₃ + δW₃
                pnext₄ = pcurr₄ - δW₄
            else
                pnext₁ = pcurr₁ - δW₁
                pnext₂ = pcurr₂ + δW₂
                pnext₃ = pcurr₃ - δW₃
                pnext₄ = pcurr₄ + δW₄
            end
            #
            pnext₁ < 1     && (pnext₁ = pnext₁ + nfine)
            pnext₁ > nfine && (pnext₁ = pnext₁ - nfine)
            pnext₂ < 1     && (pnext₂ = pnext₂ + nfine)
            pnext₂ > nfine && (pnext₂ = pnext₂ - nfine)
            pnext₃ < 1     && (pnext₃ = pnext₃ + nfine)
            pnext₃ > nfine && (pnext₃ = pnext₃ - nfine)
            pnext₄ < 1     && (pnext₄ = pnext₄ + nfine)
            pnext₄ > nfine && (pnext₄ = pnext₄ - nfine)
        else
            pnext₁ = rand(MC.rng, 1:nfine)
            pnext₂ = rand(MC.rng, 1:nfine)
            pnext₃ = rand(MC.rng, 1:nfine)
            pnext₄ = rand(MC.rng, 1:nfine)
        end

        # Apply the constraints
        !(pnext₁ in SC.allow) && continue
        !(pnext₂ in SC.allow) && continue
        !(pnext₃ in SC.allow) && continue
        !(pnext₄ in SC.allow) && continue

        # Calculate the transition probability
        Knext₁ = view(SC.kernel, :, pnext₁)
        Kcurr₁ = view(SC.kernel, :, pcurr₁)
        Knext₂ = view(SC.kernel, :, pnext₂)
        Kcurr₂ = view(SC.kernel, :, pcurr₂)
        Knext₃ = view(SC.kernel, :, pnext₃)
        Kcurr₃ = view(SC.kernel, :, pcurr₃)
        Knext₄ = view(SC.kernel, :, pnext₄)
        Kcurr₄ = view(SC.kernel, :, pcurr₄)
        #
        @. Gₙ = SC.Gᵧ + SE.A * ( Knext₁ - Kcurr₁ +
                                 Knext₂ - Kcurr₂ +
                                 Knext₃ - Kcurr₃ +
                                 Knext₄ - Kcurr₄ )
        @. ΔG = Gₙ - SC.Gᵥ
        χ²new = dot(ΔG, ΔG)
        #
        prob = exp( 0.5 * (SC.χ² - χ²new) / SC.Θ )

        # Important sampling, if true, the δ functions are shifted and the
        # corresponding objects are updated.
        if rand(MC.rng) < min(prob, 1.0)
            SE.P[s₁] = pnext₁
            SE.P[s₂] = pnext₂
            SE.P[s₃] = pnext₃
            SE.P[s₄] = pnext₄
            @. SC.Gᵧ = Gₙ
            #
            SC.χ² = χ²new
            if χ²new < SC.χ²min
                SC.χ²min = χ²new
            end
            #
            MC.Qacc = MC.Qacc + 1
        end
    end
end
