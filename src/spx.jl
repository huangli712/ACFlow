#
# Project : Gardenia
# Source  : spx.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/09/30
#

#=
### *Customized Structs* : *StochPX Solver*
=#

"""
    StochPXElement

Mutable struct. It is used to record the field configurations, which will
be sampled by Monte Carlo sweeping procedure.

For the off-diagonal elements of the matrix-valued Green's function, the
signs of the poles (𝕊) could be negative (-1.0). However, for the other
cases, 𝕊 is always positive (+1.0).

### Members
* P  -> It means the positions of the poles.
* A  -> It means the weights / amplitudes of the poles.
* 𝕊  -> It means the signs of the poles.
"""
mutable struct StochPXElement
    P  :: Vector{I64}
    A  :: Vector{F64}
    𝕊  :: Vector{F64}
end

"""
    StochPXContext

Mutable struct. It is used within the StochPX solver only.

Note that χ² denotes the goodness-of-fit functional, and Gᵧ denotes the
reproduced correlator. They should be always compatible with P, A, and 𝕊
in the StochPXElement struct.

### Members
* Gᵥ    -> Input data for correlator.
* Gᵧ    -> Generated correlator.
* σ¹    -> Actually 1.0 / σ¹.
* allow -> Allowable indices.
* grid  -> Grid for input data.
* mesh  -> Mesh for output spectrum.
* fmesh -> Very dense mesh for the poles.
* Λ     -> Precomputed kernel matrix.
* Θ     -> Artificial inverse temperature.
* χ²    -> Goodness-of-fit functional for the current configuration.
* χ²ᵥ   -> Vector of goodness-of-fit functional.
* Pᵥ    -> Vector of poles' positions.
* Aᵥ    -> Vector of poles' amplitudes.
* 𝕊ᵥ    -> Vector of poles' signs.
"""
mutable struct StochPXContext
    Gᵥ    :: Vector{F64}
    Gᵧ    :: Vector{F64}
    σ¹    :: Vector{F64}
    allow :: Vector{I64}
    grid  :: AbstractGrid
    mesh  :: AbstractMesh
    fmesh :: AbstractMesh
    Λ     :: Array{F64,2}
    Θ     :: F64
    χ²    :: F64
    χ²ᵥ   :: Vector{F64}
    Pᵥ    :: Vector{Vector{I64}}
    Aᵥ    :: Vector{Vector{F64}}
    𝕊ᵥ    :: Vector{Vector{F64}}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochPXSolver, rd::RawData)

Solve the analytic continuation problem by the stochastic pole expansion.
It is the driver for the StochPX solver. Note that this solver is still
`experimental`. It is useful for analytic continuation of Matsubara data.
In other words, this solver should not be used to handle the imagnary-time
Green's function.

Similar to the BarRat and NevanAC solvers, this solver always returns A(ω).

### Arguments
* S -> A StochPXSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* mesh -> Real frequency mesh, ω.
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
"""
function solve(S::StochPXSolver, rd::RawData)
    ktype = get_b("ktype")
    ngrid = get_b("ngrid")
    nmesh = get_b("nmesh")

    println("[ StochPX ]")
    MC, SE, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
        #
        println("Using $(nworkers()) workers")
        #
        # Copy configuration dicts
        p1 = deepcopy(PBASE)
        p2 = deepcopy(PStochPX)
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
        Gout = zeros(C64, nmesh)
        if ktype == "bsymm"
            Gᵣ = zeros(F64, ngrid)
        else
            Gᵣ = zeros(F64, 2 * ngrid)
        end
        for i in eachindex(sol)
            a, b, c = sol[i]
            @. Aout = Aout + a / nworkers()
            @. Gout = Gout + b / nworkers()
            @. Gᵣ = Gᵣ + c / nworkers()
        end
        #
        # Postprocess the solutions
        last(SC, Aout, Gout, Gᵣ)
        #
    # Sequential version
    else
        #
        Aout, Gout, Gᵣ = run(MC, SE, SC)
        last(SC, Aout, Gout, Gᵣ)
        #
    end

    return SC.mesh.mesh, Aout, Gout
end

"""
    init(S::StochPXSolver, rd::RawData)

Initialize the StochPX solver and return the StochPXMC, StochPXElement,
and StochPXContext structs. Please don't call this function directly.

### Arguments
* S -> A StochPXSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.
"""
function init(S::StochPXSolver, rd::RawData)
    # Initialize possible constraints.
    # The array arrow contains all the possible indices for poles.
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
    Gᵧ, Λ, Θ, χ², χ²ᵥ, Pᵥ, Aᵥ, 𝕊ᵥ = init_context(S, SE, grid, fmesh, Gᵥ)
    SC = StochPXContext(Gᵥ, Gᵧ, σ¹, allow, grid, mesh, fmesh,
                        Λ, Θ, χ², χ²ᵥ, Pᵥ, Aᵥ, 𝕊ᵥ)
    println("Initialize context for the StochPX solver")

    return MC, SE, SC
end

"""
    run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Perform stochastic pole expansion simulation, sequential version.

### Arguments
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
* Gᵣ -> Reproduced Green's function, G(iωₙ).
"""
function run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Setup essential parameters
    ntry = get_x("ntry")
    nstep = get_x("nstep")

    # Warmup the Monte Carlo engine
    println("Start thermalization...")
    for _ = 1:nstep
        sample(1, MC, SE, SC)
    end

    # Sample and collect data
    println("Start stochastic sampling...")
    for t = 1:ntry
        # Reset Monte Carlo counters
        reset_mc(MC)

        # Reset Monte Carlo field configuration
        reset_element(MC.rng, SC.allow, SE)

        # Reset Gᵧ and χ² in SC (StochPXContext)
        reset_context(t, SE, SC)

        # Apply simulated annealing algorithm
        for _ = 1:nstep
            sample(t, MC, SE, SC)
        end

        # Write Monte Carlo statistics
        fwrite && write_statistics(MC)

        # Show the best χ² (the smallest) for the current attempt
        @printf("try = %6i -> [χ² = %9.4e]\n", t, SC.χ²ᵥ[t])
        flush(stdout)
    end

    # Write pole expansion coefficients
    fwrite && write_pole(SC.Pᵥ, SC.Aᵥ, SC.𝕊ᵥ, SC.χ²ᵥ, SC.fmesh)

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    prun(
        S::StochPXSolver,
        p1::Dict{String,Vector{Any}},
        p2::Dict{String,Vector{Any}},
        MC::StochPXMC,
        SE::StochPXElement,
        SC::StochPXContext
    )

Perform stochastic pole expansion simulation, parallel version.
The arguments `p1` and `p2` are copies of PBASE and PStochPX, respectively.

### Arguments
* S -> A StochPXSolver struct.
* p1 -> A copy of PBASE.
* p2 -> A copy of PStochPX.
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
* Gᵣ -> Reproduced Green's function, G(iωₙ).
"""
function prun(
    S::StochPXSolver,
    p1::Dict{String,Vector{Any}},
    p2::Dict{String,Vector{Any}},
    MC::StochPXMC,
    SE::StochPXElement,
    SC::StochPXContext
    )
    # Revise parameteric dicts
    # We have to make sure that all processes share the same parameters.
    rev_dict_b(p1)
    rev_dict_x(S, p2)

    # Initialize random number generator again
    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Setup essential parameters
    ntry = get_x("ntry")
    nstep = get_x("nstep")

    # Warmup the Monte Carlo engine
    println("Start thermalization...")
    for _ = 1:nstep
        sample(1, MC, SE, SC)
    end

    # Sample and collect data
    println("Start stochastic sampling...")
    for t = 1:ntry
        # Reset Monte Carlo counters
        reset_mc(MC)

        # Reset Monte Carlo field configuration
        reset_element(MC.rng, SC.allow, SE)

        # Reset Gᵧ and χ² in SC (StochPXContext)
        reset_context(t, SE, SC)

        # Apply simulated annealing algorithm
        for _ = 1:nstep
            sample(t, MC, SE, SC)
        end

        # Write Monte Carlo statistics
        myid() == 2 && fwrite && write_statistics(MC)

        # Show the best χ² (the smallest) for the current attempt
        @printf("try = %6i -> [χ² = %9.4e]\n", t, SC.χ²ᵥ[t])
        flush(stdout)
    end

    # Write pole expansion coefficients
    myid() == 2 && fwrite && write_pole(SC.Pᵥ, SC.Aᵥ, SC.𝕊ᵥ, SC.χ²ᵥ, SC.fmesh)

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    average(SC::StochPXContext)

Postprocess the results generated during the stochastic pole expansion
simulations. It will calculate the spectral functions, real frequency
Green's function, and imaginary frequency Green's function.

### Arguments
* SC -> A StochPXContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
* Gᵣ -> Reproduced Green's function, G(iωₙ).
"""
function average(SC::StochPXContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Setup essential parameters
    nmesh = get_b("nmesh")
    method = get_x("method")
    ntry = get_x("ntry")

    # Allocate memory
    # Gout: real frequency Green's function, G(ω).
    # Gᵣ: imaginary frequency Green's function, G(iωₙ)
    ngrid, _ = size(SC.Λ)
    Gout = zeros(C64, nmesh)
    Gᵣ = zeros(F64, ngrid)

    # Choose the best solution
    if method == "best"
        # The χ² of the best solution should be the smallest.
        p = argmin(SC.χ²ᵥ)
        @printf("Best solution: try = %6i -> [χ² = %9.4e]\n", p, SC.χ²ᵥ[p])
        #
        # Calculate G(ω)
        Gout = calc_green(p, SC, true)
        #
        # Calculate G(iωₙ)
        Gᵣ = calc_green(p, SC, false)
    #
    # Collect the `good` solutions and calculate their average.
    else
        # Calculate the median of SC.χ²ᵥ
        chi2_med = median(SC.χ²ᵥ)
        chi2_ave = mean(SC.χ²ᵥ)

        # Determine the αgood parameter, which is used to filter the
        # calculated spectra.
        αgood = 1.2
        if count(x -> x < chi2_med / αgood, SC.χ²ᵥ) ≤ ntry / 10
            αgood = 1.0
        end

        # Go through all the solutions
        c = 0.0 # A counter
        passed = I64[]
        for i = 1:ntry
            if SC.χ²ᵥ[i] < chi2_med / αgood
                # Calculate and accumulate G(ω)
                G = calc_green(i, SC, true)
                @. Gout = Gout + G
                #
                # Calculate and accumulate G(iωₙ)
                G = calc_green(i, SC, false)
                @. Gᵣ = Gᵣ + G
                #
                # Increase the counter
                c = c + 1.0
                append!(passed, i)
            end
        end
        #
        # Normalize the final results
        @. Gout = Gout / c
        @. Gᵣ = Gᵣ / c
        println("Mean value of χ²: $(chi2_ave)")
        println("Median value of χ²: $(chi2_med)")
        println("Accumulate $(round(I64,c)) solutions to get the spectral density")
        #
        # Write indices of selected solutions
        if nworkers() > 1
            myid() == 2 && fwrite && write_passed(passed, chi2_med, αgood)
        else
            fwrite && write_passed(passed, chi2_med, αgood)
        end
    #
    end

    return -imag.(Gout) / π, Gout, Gᵣ
end

"""
    last(
        SC::StochPXContext,
        Aout::Vector{F64},
        Gout::Vector{C64},
        Gᵣ::Vector{F64}
    )

It will write the calculated results by the StochPX solver, including the
final spectral function and reproduced correlator.

### Arguments
* SC   -> A StochPXContext struct.
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
* Gᵣ   -> Reconstructed Green's function, G(iωₙ).

### Returns
N/A
"""
function last(
    SC::StochPXContext,
    Aout::Vector{F64},
    Gout::Vector{C64},
    Gᵣ::Vector{F64}
    )
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Write the spectral function, A(ω).
    fwrite && write_spectrum(SC.mesh, Aout)

    # Reproduce input data and write them, G(iωₙ).
    fwrite && write_backward(SC.grid, Gᵣ)

    # Write full response function on real axis, G(ω).
    fwrite && write_complete(SC.mesh, Gout)
end

#=
### *Core Algorithms*
=#

"""
    sample(
        t::I64,
        MC::StochPXMC,
        SE::StochPXElement,
        SC::StochPXContext
    )

Try to search the configuration space to locate the minimum by using the
simulated annealing algorithm. Here, `t` means the t-th attempt.

### Arguments
* t -> Counter for attempts.
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
N/A
"""
function sample(
    t::I64,
    MC::StochPXMC,
    SE::StochPXElement,
    SC::StochPXContext
    )
    # Try to change positions of poles
    if rand(MC.rng) < 0.5
        if rand(MC.rng) < 0.9
            try_move_s(t, MC, SE, SC)
        else
            try_move_p(t, MC, SE, SC)
        end
    # Try to change amplitudes of poles
    else
        if rand(MC.rng) < 0.5
            try_move_a(t, MC, SE, SC)
        else
            try_move_x(t, MC, SE, SC)
        end
    end
end

"""
    measure(t::I64, SE::StochPXElement, SC::StochPXContext)

Store Monte Carlo field configurations (positions, amplitudes, and signs
of many poles) for the `t`-th attempt. In other words, the current field
configuration (recorded in `SE`) will be saved in `SC`.

Note that not all configurations for the `t`-th attempt will be saved.
Only the solution that exhibits the smallest χ² will be saved. For the
`t`-th attempt, the StochPX solver will do `nstep` Monte Carlo updates.
It will calculate all χ², and try to figure out the smallest one. Then
the corresponding configuration (solution) will be saved.

### Arguments
* t -> Counter for the attemps.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
N/A
"""
function measure(t::I64, SE::StochPXElement, SC::StochPXContext)
    SC.χ²ᵥ[t] = SC.χ²
    #
    @. SC.Pᵥ[t] = SE.P
    @. SC.Aᵥ[t] = SE.A
    @. SC.𝕊ᵥ[t] = SE.𝕊
end

#=
### *Service Functions*
=#

"""
    init_iodata(S::StochPXSolver, rd::RawData)

Preprocess the input data (`rd`).

### Arguments
* S -> A StochPXSolver struct.
* rd -> A RawData struct, which contains essential input data.

### Returns
* Gᵥ -> Input correlator.
* σ¹ -> 1.0 / σ¹.

See also: [`RawData`](@ref), [`GreenData`](@ref).
"""
function init_iodata(S::StochPXSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value # Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹
end

"""
    init_mc(S::StochPXSolver)

Try to create a StochPXMC struct. Some counters for Monte Carlo updates
are initialized here.

### Arguments
* S -> A StochPXSolver struct.

### Returns
* MC -> A StochPXMC struct.

See also: [`StochPXMC`](@ref).
"""
function init_mc(S::StochPXSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    #
    Sacc = 0
    Stry = 0
    #
    Pacc = 0
    Ptry = 0
    #
    Aacc = 0
    Atry = 0
    #
    Xacc = 0
    Xtry = 0
    #
    MC = StochPXMC(rng, Sacc, Stry, Pacc, Ptry, Aacc, Atry, Xacc, Xtry)

    return MC
end

"""
    init_element(
        S::StochPXSolver,
        rng::AbstractRNG,
        allow::Vector{I64}
    )

Randomize the configurations for future Monte Carlo sampling. It will
return a StochPXElement struct. Note that `allow` is generated in the
`constraints()` function.

### Arguments
* S     -> A StochPXSolver struct.
* rng   -> Random number generator.
* allow -> Allowed positions for the poles.

### Returns
* SE -> A StochPXElement struct.

See also: [`StochPXElement`](@ref).
"""
function init_element(
    S::StochPXSolver,
    rng::AbstractRNG,
    allow::Vector{I64}
    )
    offdiag = get_b("offdiag")
    npole = get_x("npole")

    if offdiag
        # We just assume that the numbers of poles for the positive and
        # negative parts are equal.
        @assert iseven(npole)
        hpole = npole ÷ 2

        # We have to make sure that both negative weights and positive
        # weights are allowed for the generated poles.
        @assert any(x -> x > 0.0, allow)
        @assert any(x -> x < 0.0, allow)

        # Initialize P, A, and 𝕊
        allow₊ = filter(x -> x > 0.0, allow)
        allow₋ = filter(x -> x < 0.0, allow)
        P₊ = rand(rng, allow₊, hpole)
        P₋ = rand(rng, allow₋, hpole)
        A₊ = rand(rng, F64, hpole)
        A₋ = rand(rng, F64, hpole)
        𝕊₊ = ones(F64, hpole)
        𝕊₋ = ones(F64, hpole) * (-1.0)

        # We have to make sure ∑ᵢ Aᵢ = 1
        s = sum(A₊)
        @. A₊ = A₊ / s
        s = sum(A₋)
        @. A₋ = A₋ / s

        # Merge positive and negative parts
        P = vcat(P₊, P₋)
        A = vcat(A₊, A₋)
        𝕊 = vcat(𝕊₊, 𝕊₋)
    else
        # Initialize P, A, and 𝕊
        P = rand(rng, allow, npole)
        A = rand(rng, F64, npole)
        𝕊 = ones(F64, npole)

        # We have to make sure ∑ᵢ Aᵢ = 1
        s = sum(A)
        @. A = A / s
    end

    SE = StochPXElement(abs.(P), A, 𝕊)

    return SE
end

"""
    init_context(
        S::StochPXSolver,
        SE::StochPXElement,
        grid::AbstractGrid,
        fmesh::AbstractMesh,
        Gᵥ::Vector{F64}
    )

Try to initialize the key members of a StochPXContext struct. It will try
to return some key variables, which should be used to construct the
StochPXContext struct.

### Arguments
* S     -> A StochPXSolver struct.
* SE    -> A StochPXElement struct.
* grid  -> Grid for input correlator.
* fmesh -> Fine mesh in [wmin, wmax], used to build the kernel matrix Λ.
* Gᵥ    -> Preprocessed input correlator.

### Returns
* Gᵧ  -> Reconstructed correlator.
* Λ   -> Precomputed kernel matrix.
* Θ   -> Artificial inverse temperature.
* χ²  -> Current goodness-of-fit functional.
* χ²ᵥ -> Vector of goodness-of-fit functional.
* Pᵥ  -> Vector of poles' positions.
* Aᵥ  -> Vector of poles' amplitudes.
* 𝕊ᵥ  -> Vector of poles' signs.

See also: [`StochPXContext`](@ref).
"""
function init_context(
    S::StochPXSolver,
    SE::StochPXElement,
    grid::AbstractGrid,
    fmesh::AbstractMesh,
    Gᵥ::Vector{F64}
    )
    # Extract some parameters
    ntry = get_x("ntry")
    npole = get_x("npole")
    Θ = get_x("theta")

    # Prepare the kernel matrix Λ. It is used to speed up the simulation.
    # Note that Λ depends on the type of kernel.
    Λ = calc_lambda(grid, fmesh, Gᵥ)

    # We have to make sure that the starting Gᵧ and χ² are consistent with
    # the current Monte Carlo configuration fields.
    Gᵧ = calc_green(SE.P, SE.A, SE.𝕊, Λ)
    χ² = calc_chi2(Gᵧ, Gᵥ)

    # χ²ᵥ is initialized by a large number. Later it will be updated by
    # the smallest χ² during the simulation.
    χ²ᵥ = zeros(F64, ntry)
    @. χ²ᵥ = 1e10

    # P, A, and 𝕊 should be always compatible with χ². They are updated
    # in the `measure()` function.
    Pᵥ = Vector{I64}[]
    Aᵥ = Vector{F64}[]
    𝕊ᵥ = Vector{F64}[]
    #
    for _ = 1:ntry
        push!(Pᵥ,  ones(I64, npole))
        push!(Aᵥ, zeros(F64, npole))
        push!(𝕊ᵥ, zeros(F64, npole))
    end

    return Gᵧ, Λ, Θ, χ², χ²ᵥ, Pᵥ, Aᵥ, 𝕊ᵥ
end

"""
    reset_mc(MC::StochPXMC)

Reset the counters in StochPXMC struct.

### Arguments
* MC -> A StochPXMC struct.

### Returns
N/A
"""
function reset_mc(MC::StochPXMC)
    MC.Sacc = 0
    MC.Stry = 0
    #
    MC.Pacc = 0
    MC.Ptry = 0
    #
    MC.Aacc = 0
    MC.Atry = 0
    #
    MC.Xacc = 0
    MC.Xtry = 0
end

"""
    reset_element(
        rng::AbstractRNG,
        allow::Vector{I64},
        SE::StochPXElement
    )

Reset the Monte Carlo field configurations (i.e. positions and amplitudes
of the poles). Note that the signs of the poles should not be changed.
Be careful, the corresponding χ² (goodness-of-fit functional) and Gᵧ
(generated correlator) will be updated in the `reset_context()` function.

### Arguments
* rng   -> Random number generator.
* allow -> Allowed positions for the poles.
* SE -> A StochPXElement struct.

### Returns
N/A
"""
function reset_element(
    rng::AbstractRNG,
    allow::Vector{I64},
    SE::StochPXElement
    )
    offdiag = get_b("offdiag")
    npole = get_x("npole")

    # For off-diagonal elements
    if offdiag
        # The number of poles must be even.
        # Here, `hpole` means half number of poles.
        hpole = npole ÷ 2

        # How many poles should be changed
        if hpole ≤ 5
            if 4 ≤ hpole ≤ 5
                hselect = 2
            else
                hselect = 1
            end
        else
            hselect = hpole ÷ 5
        end
        @assert hselect ≤ hpole

        # Which poles should be changed
        # For poles that with positive weights
        selected₊ = rand(rng, 1:hpole, hselect)
        unique!(selected₊)
        hselect₊ = length(selected₊)
        #
        # For poles that with negative weights
        selected₋ = rand(rng, hpole+1:npole, hselect)
        unique!(selected₋)
        hselect₋ = length(selected₋)

        # Change poles' positions
        if rand(rng) < 0.9
            allow₊ = filter(x -> x > 0.0, allow)
            P₊ = rand(rng, allow₊, hselect₊)
            @. SE.P[selected₊] = abs.(P₊)
            #
            allow₋ = filter(x -> x < 0.0, allow)
            P₋ = rand(rng, allow₋, hselect₋)
            @. SE.P[selected₋] = abs.(P₋)
        # Change poles' amplitudes
        else
            # For positive-weight poles
            A₁₊ = SE.A[selected₊]
            s₁₊ = sum(A₁₊)
            #
            A₂₊ = rand(rng, F64, hselect₊)
            s₂₊ = sum(A₂₊)
            @. A₂₊ = A₂₊ / s₂₊ * s₁₊
            #
            @. SE.A[selected₊] = A₂₊

            # For negative-weight poles
            A₁₋ = SE.A[selected₋]
            s₁₋ = sum(A₁₋)
            #
            A₂₋ = rand(rng, F64, hselect₋)
            s₂₋ = sum(A₂₋)
            @. A₂₋ = A₂₋ / s₂₋ * s₁₋
            #
            @. SE.A[selected₋] = A₂₋
        end
    # For diagonal elements
    else
        # How many poles should be changed
        if npole ≤ 5
            if 4 ≤ npole ≤ 5
                nselect = 2
            else
                nselect = 1
            end
        else
            nselect = npole ÷ 5
        end
        @assert nselect ≤ npole

        # Which poles should be changed
        selected = rand(rng, 1:npole, nselect)
        unique!(selected)
        nselect = length(selected)

        # Change poles' positions
        if rand(rng) < 0.5
            P = rand(rng, allow, nselect)
            @. SE.P[selected] = P
        # Change poles' amplitudes
        else
            A₁ = SE.A[selected]
            s₁ = sum(A₁)
            #
            A₂ = rand(rng, F64, nselect)
            s₂ = sum(A₂)
            @. A₂ = A₂ / s₂ * s₁
            #
            @. SE.A[selected] = A₂
        end
    end
end

"""
    reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)

Recalculate imaginary frequency Green's function Gᵧ and goodness-of-fit
functional χ² by new Monte Carlo field configurations for the `t`-th
attempts. They must be consistent with each other.

Some key variables in `SC` are also updated as well. Perhaps we should
develop a smart algorhtm to update Θ here.

### Arguments
* t -> Counter for attempts.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
N/A
"""
function reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)
    SC.Θ = get_x("theta")
    SC.Gᵧ = calc_green(SE.P, SE.A, SE.𝕊, SC.Λ)
    SC.χ² = calc_chi2(SC.Gᵧ, SC.Gᵥ)
    SC.χ²ᵥ[t] = 1e10
end

"""
    calc_fmesh(S::StochPXSolver)

Try to calculate very fine (dense) mesh in [wmin, wmax], which is used
internally to represent the possible positions of poles. Note that this
mesh could be non-uniform. If the file `fmesh.inp` exists, the code will
try to load it to initialize the mesh. Or else the code will generate
a linear mesh automatically.

### Arguments
* S -> A StochPXSolver struct.

### Returns
* fmesh -> A very fine, perhaps non-uniform mesh in [wmin, wmax].

See also: [`LinearMesh`](@ref), [`DynamicMesh`](@ref).
"""
function calc_fmesh(S::StochPXSolver)
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nfine = get_x("nfine")

    # Filename for the predefined mesh
    # This file should contain at least `nfine` lines
    fn = "fmesh.inp"

    # If the file `fmesh.inp` exists, we will use it to build the mesh.
    if isfile(fn)
        mesh = zeros(F64, nfine)
        #
        open(fn, "r") do fin
            for i = 1:nfine
                arr = line_to_array(fin)
                mesh[i] = parse(F64, arr[2])
            end
        end
        #
        fmesh = DynamicMesh(mesh)
    # Or else we will return a linear mesh directly.
    else
        fmesh = LinearMesh(nfine, wmin, wmax)
    end

    return fmesh
end

#=
*Remarks* :

Here, we would like to discuss the kernel matrix and sum-rule in the
stochastic pole expansion algorithm.

**For Fermionic Green's Functions**

The pole expansion of fermionic Green's function reads:

```math
\begin{equation}
G(i\omega_n) = \sum^{N_p}_{i = 1} \frac{A_i}{i\omega_n - P_i}.
\end{equation}
```

The kernel matrix reads:

```math
\begin{equation}
Λ(\omega_n, \omega) = \frac{1}{i\omega_n - \omega}.
\end{equation}
```

The sum-rule reads:

```math
\begin{equation}
\sum^{N_p}_{i=1} A_i = 1.
\end{equation}
```

---

**For Bosonic Green's Functions**

Supposed that ``\chi_0 \equiv -G(i\omega_n = 0)``, which must be a
positive real number, then the pole expansion of bosonic Green's
function reads:

```math
\begin{equation}
G(i\omega_n) = \sum^{N_p}_{i = 1} \frac{A_i}{i\omega_n - P_i}
= \sum^{N_p}_{i=1}
\frac{\chi_0 P_i}{i\omega_n - P_i}
\frac{A_i}{\chi_0 P_i}.
\end{equation}
```

The kernel matrix reads:

```math
\begin{equation}
Λ(\omega_n, \omega) = \frac{\chi_0 \omega}{i\omega_n - \omega}.
\end{equation}
```

Be careful, when ``\omega_n = 0``, ``\Lambda(0,\omega) = -\chi_0``.

The sum-rule reads:

```math
\begin{equation}
\sum^{N_p}_{i=1} \tilde{A}_i = 1,
\end{equation}
```

where

```math
\begin{equation}
\tilde{A}_i = \frac{A_i}{\chi_0 P_i}.
\end{equation}
```

---

**For Correlator of Hermitian Operator**

The pole expansion for correlator of Hermitian operator reads:

```math
\begin{equation}
G(i\omega_n) = \sum^{N_p}_{i = 1} A_i
\left(\frac{1}{i\omega_n - P_i} - \frac{1}{i\omega_n + P_i}\right)
= \sum^{N_p}_{i =1}
\frac{-\chi_0 P^2_i}{\omega^2_n + P^2_i} \frac{2A_i}{\chi_0 P_i}.
\end{equation}
```

Here, ``\chi_0 \equiv -G(i\omega_n = 0)``, which is a positive real
number. And ``\forall i, A_i \ge 0, P_i \ge 0``.

The kernel matrix reads:

```math
\begin{equation}
Λ(\omega_n, \omega) = \frac{-\chi_0 \omega^2}{\omega^2_n + \omega^2}.
\end{equation}
```

Be careful, when ``\omega_n = 0``, ``\Lambda(0,\omega) = -\chi_0``.

The sum-rule reads:

```math
\begin{equation}
\sum^{N_p}_{i=1} \tilde{A}_i = 1,
\end{equation}
```

where

```math
\begin{equation}
\tilde{A}_i = \frac{2A_i}{\chi_0 P_i}.
\end{equation}
```

---

=#

"""
    calc_lambda(
        grid::AbstractGrid,
        fmesh::AbstractMesh,
        Gᵥ::Vector{F64}
    )

Precompute the kernel matrix Λ (Λ ≡ 1 / (iωₙ - ϵ)). It is the essential
driver function. Note that Λ depends on the kernel's type (`ktype`).

### Arguments
* grid  -> Imaginary axis grid for input data.
* fmesh -> Very dense mesh in [wmin, wmax].
* Gᵥ -> Preprocessed input correlator.

### Returns
* Λ -> The kernel matrix, a 2D array.
"""
function calc_lambda(
    grid::AbstractGrid,
    fmesh::AbstractMesh,
    Gᵥ::Vector{F64}
    )
    ktype = get_b("ktype")
    χ₀ = -Gᵥ[1]

    @cswitch ktype begin
        #
        @case "fermi"
            Λ = calc_lambda(grid, fmesh)
            break
        #
        @case "boson"
            Λ = calc_lambda(grid, fmesh, χ₀, false)
            break
        #
        @case "bsymm"
            Λ = calc_lambda(grid, fmesh, χ₀, true)
            break
        #
    end

    return Λ
end

"""
    calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh)

Precompute the kernel matrix Λ (Λ ≡ 1 / (iωₙ - ϵ)). It is a service
function and is for the fermionic systems.

### Arguments
* grid  -> Imaginary axis grid for input data.
* fmesh -> Very dense mesh in [wmin, wmax].

### Returns
* Λ -> The kernel matrix, a 2D array.
"""
function calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh)
    ngrid = get_b("ngrid")
    nfine = get_x("nfine")

    _Λ = zeros(C64, ngrid, nfine)
    #
    for i in eachindex(grid)
        iωₙ = im * grid[i]
        for j in eachindex(fmesh)
            _Λ[i,j] = 1.0 / (iωₙ - fmesh[j])
        end
    end
    #
    Λ = vcat(real(_Λ), imag(_Λ))

    return Λ
end

"""
    calc_lambda(
        grid::AbstractGrid,
        fmesh::AbstractMesh,
        χ₀::F64,
        bsymm::Bool
    )

Precompute the kernel matrix Λ. Here, `χ₀` is actually -G(iωₙ = 0). And
the argument `bsymm` is used to distinguish two different bosonic kernels.
If `bsymm` is false, it means that the kernel is `boson`. If `bsymm` is
true, the kernel is `bsymm`. This function is for the bosonic systems.

### Arguments
* grid  -> Imaginary axis grid for input data.
* fmesh -> Very dense mesh in [wmin, wmax].
* χ₀    -> -G(0).
* bsymm -> Type of bosonic kernel.

### Returns
* Λ -> The kernel matrix, a 2D array.
"""
function calc_lambda(
    grid::AbstractGrid,
    fmesh::AbstractMesh,
    χ₀::F64,
    bsymm::Bool
    )
    ngrid = get_b("ngrid")
    nfine = get_x("nfine")

    # For standard bosonic kernel matrix
    if bsymm == false

        _Λ = zeros(C64, ngrid, nfine)
        #
        for i in eachindex(grid)
            iωₙ = im * grid[i]
            for j in eachindex(fmesh)
                _Λ[i,j] = χ₀ * fmesh[j] / (iωₙ - fmesh[j])
            end
        end
        #
        # Special treatment for iωₙ = 0
        for j in eachindex(fmesh)
            _Λ[1,j] = -χ₀
        end

        Λ = vcat(real(_Λ), imag(_Λ))

    # For symmetric bosonic kernel matrix
    else

        _Λ = zeros(F64, ngrid, nfine)
        #
        for i in eachindex(grid)
            ωₙ = grid[i]
            for j in eachindex(fmesh)
                _Λ[i,j] = -χ₀ * (fmesh[j] ^ 2.0) / (ωₙ ^ 2.0 + fmesh[j] ^ 2.0)
            end
        end
        #
        # Special treatment for ωₙ = 0
        for j in eachindex(fmesh)
            _Λ[1,j] = -χ₀
        end
        #
        Λ = copy(_Λ)

    end

    return Λ
end

"""
    calc_green(t::I64, SC::StochPXContext, real_axis::Bool)

Reconstruct Green's function at imaginary axis or real axis by using the
pole expansion. It is a driver function. If `real_axis = true`, it will
returns G(ω), or else G(iωₙ).

### Arguments
* t -> Index of the current attempt.
* SC -> A StochPXContext struct.
* real_axis -> Working at real axis (true) or imaginary axis (false)?

### Returns
* G -> Reconstructed Green's function, G(ω) or G(iωₙ).
"""
function calc_green(t::I64, SC::StochPXContext, real_axis::Bool)
    ktype = get_b("ktype")
    ntry = get_x("ntry")
    @assert t ≤ ntry

    # Calculate G(iωₙ)
    if real_axis == false
        return calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.Λ)
    end

    # Calculate G(ω). Now we don't need SC.Λ.
    χ₀ = -SC.Gᵥ[1]
    @cswitch ktype begin
        @case "fermi"
            G = calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.mesh, SC.fmesh)
            break

        @case "boson"
            G = calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.mesh, SC.fmesh, χ₀, false)
            break

        @case "bsymm"
            G = calc_green(SC.Pᵥ[t], SC.Aᵥ[t], SC.𝕊ᵥ[t], SC.mesh, SC.fmesh, χ₀, true)
            break
    end

    return G
end

"""
    calc_green(
        P::Vector{I64},
        A::Vector{F64},
        𝕊::Vector{F64},
        Λ::Array{F64,2}
    )

Reconstruct Green's function at imaginary axis by the pole expansion.

### Arguments
* P -> Positions of poles.
* A -> Amplitudes of poles.
* 𝕊 -> Signs of poles.
* Λ -> Kernel matrix Λ.

### Returns
* G -> Reconstructed Green's function, G(iωₙ).
"""
function calc_green(
    P::Vector{I64},
    A::Vector{F64},
    𝕊::Vector{F64},
    Λ::Array{F64,2}
    )
    # Note that here `ngrid` is equal to 2 × ngrid sometimes.
    ngrid, _ = size(Λ)

    G = zeros(F64, ngrid)
    for i = 1:ngrid
        G[i] = dot(A .* 𝕊, Λ[i,P])
    end

    return G
end

"""
    calc_green(
        P::Vector{I64},
        A::Vector{F64},
        𝕊::Vector{F64},
        mesh::AbstractMesh,
        fmesh::AbstractMesh
    )

Reconstruct Green's function at real axis by the pole expansion. It is
for the fermionic systems only.

### Arguments
* P     -> Positions of poles.
* A     -> Amplitudes of poles.
* 𝕊     -> Signs of poles.
* mesh  -> Real frequency mesh for spectral functions.
* fmesh -> Very dense real frequency mesh for poles.

### Returns
* G -> Retarded Green's function, G(ω).
"""
function calc_green(
    P::Vector{I64},
    A::Vector{F64},
    𝕊::Vector{F64},
    mesh::AbstractMesh,
    fmesh::AbstractMesh
    )
    η = get_x("eta")
    nmesh = length(mesh)

    iωₙ = mesh.mesh .+ im * η
    G = zeros(C64, nmesh)
    for i in eachindex(mesh)
        G[i] = sum( @. (A * 𝕊) / (iωₙ[i] - fmesh.mesh[P]) )
    end

    return G
end

"""
    calc_green(
        P::Vector{I64},
        A::Vector{F64},
        𝕊::Vector{F64},
        mesh::AbstractMesh,
        fmesh::AbstractMesh,
        χ₀::F64,
        bsymm::Bool
    )

Reconstruct Green's function at real axis by the pole expansion. Here,
`χ₀` is actually -G(iωₙ = 0). And the argument `bsymm` is used to
distinguish two different bosonic kernels. If `bsymm` is false, it means
that the kernel is `boson`. If `bsymm` is true, the kernel is `bsymm`.
It is for the bosonic systems only.

### Arguments
* P     -> Positions of poles.
* A     -> Amplitudes of poles.
* 𝕊     -> Signs of poles.
* mesh  -> Real frequency mesh for spectral functions.
* fmesh -> Very dense real frequency mesh for poles.
* χ₀    -> -G(0).
* bsymm -> Type of bosonic kernel.

### Returns
* G -> Retarded Green's function, G(ω).
"""
function calc_green(
    P::Vector{I64},
    A::Vector{F64},
    𝕊::Vector{F64},
    mesh::AbstractMesh,
    fmesh::AbstractMesh,
    χ₀::F64,
    bsymm::Bool
    )
    η = get_x("eta")
    nmesh = length(mesh)

    iωₙ = mesh.mesh .+ im * η
    G = zeros(C64, nmesh)
    if bsymm == false
        _A = A .* 𝕊 .* χ₀ .* fmesh.mesh[P]
        for i in eachindex(mesh)
            G[i] = sum( @. _A / (iωₙ[i] - fmesh.mesh[P]) )
        end
    #
    else
        _A = A .* 𝕊 .* χ₀ .* fmesh.mesh[P] .* 0.5
        for i in eachindex(mesh)
            G₊ = sum( @. _A / (iωₙ[i] - fmesh.mesh[P]) )
            G₋ = sum( @. _A / (iωₙ[i] + fmesh.mesh[P]) )
            G[i] = G₊ - G₋
        end
    #
    end

    return G
end

"""
    calc_chi2(Gₙ::Vector{F64}, Gᵥ::Vector{F64})

Try to calculate the goodness-of-fit function (i.e, χ²), which measures
the distance between input and regenerated correlators.

### Arguments
* Gₙ -> Reconstructed Green's function.
* Gᵥ -> Original Green's function.

### Returns
* chi2 -> Goodness-of-fit functional, χ².

See also: [`calc_green`](@ref).
"""
function calc_chi2(Gₙ::Vector{F64}, Gᵥ::Vector{F64})
    ΔG = Gₙ - Gᵥ
    return dot(ΔG, ΔG)
end

"""
    constraints(S::StochPXSolver, fmesh::AbstractMesh)

Try to implement the constrained stochastic pole expansion. This
function will return a collection. It contains all the allowable indices.
Be careful, the constrained stochastic pole expansion method is always
compatible with the self-adaptive mesh.

If `offdiag = true`, the regions defined by `exclude` are for A(ω) < 0,
the corresponding signs of poles are negative. While the other regions
that are not included in `exclude` are for A(ω) > 0, and the corresponding
signs of poles are positive.

### Arguments
* S     -> A StochPXSolver struct.
* fmesh -> Very dense mesh for the poles.

### Returns
* allow -> Allowable indices.

See also: [`StochPXSolver`](@ref).
"""
function constraints(S::StochPXSolver, fmesh::AbstractMesh)
    offdiag = get_b("offdiag")
    exclude = get_b("exclude")
    nfine = get_x("nfine")
    @assert nfine == length(fmesh)

    allow = I64[]
    unallow = I64[]

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
        else
            push!(unallow, -i)
        end
    end

    # If it is offdiagonal, then the spectral function can be negative.
    # Now `allow` is for A(ω) > 0 and `unallow` is for A(ω) < 0. We have
    # to distinguish them.
    if offdiag
        append!(allow, unallow)
    end

    return allow
end

"""
    try_move_s(
        t::I64,
        MC::StochPXMC,
        SE::StochPXElement,
        SC::StochPXContext
    )

Change the position of one randomly selected pole.

### Arguments
* t -> Counter of attempts.
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
N/A

See also: [`try_move_p`](@ref).
"""
function try_move_s(
    t::I64,
    MC::StochPXMC,
    SE::StochPXElement,
    SC::StochPXContext
    )
    # Get parameters
    ngrid = length(SC.Gᵥ) # get_b("ngrid")
    nfine = get_x("nfine")
    npole = get_x("npole")
    move_window = nfine ÷ 100

    # It is used to save the change of Green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

        # Select one pole randomly
        s = rand(MC.rng, 1:npole)

        # Try to change position of the s pole
        Aₛ = SE.A[s]
        𝕊ₛ = SE.𝕊[s]
        #
        δP = rand(MC.rng, 1:move_window)
        #
        P₁ = SE.P[s]
        P₂ = P₁
        if rand(MC.rng) > 0.5
            P₂ = P₁ + δP
        else
            P₂ = P₁ - δP
        end
        #
        if 𝕊ₛ > 0.0
            !(+P₂ in SC.allow) && continue
        else
            !(-P₂ in SC.allow) && continue
        end

        # Calculate change of Green's function
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        @. δG = 𝕊ₛ * Aₛ * (Λ₂ - Λ₁)

        # Calculate new Green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²

        # Simulated annealing algorithm
        MC.Stry = MC.Stry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.P[s] = P₂

            # Update reconstructed Green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ² = χ²

            # Update Monte Carlo counter
            MC.Sacc = MC.Sacc + 1

            # Save optimal solution
            if SC.χ² < SC.χ²ᵥ[t]
                measure(t, SE, SC)
            end
        end

    end
end

"""
    try_move_p(
        t::I64,
        MC::StochPXMC,
        SE::StochPXElement,
        SC::StochPXContext
    )

Change the positions of two randomly selected poles.

### Arguments
* t -> Counter of attempts.
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
N/A

See also: [`try_move_s`](@ref).
"""
function try_move_p(
    t::I64,
    MC::StochPXMC,
    SE::StochPXElement,
    SC::StochPXContext
    )
    # Get parameters
    ngrid = length(SC.Gᵥ) # get_b("ngrid")
    npole = get_x("npole")

    # Sanity check
    if npole == 1
        return
    end

    # It is used to save the change of Green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

        # Select two poles randomly
        # The two poles should not be the same.
        s₁ = 1
        s₂ = 1
        while s₁ == s₂
            s₁ = rand(MC.rng, 1:npole)
            s₂ = rand(MC.rng, 1:npole)
        end

        # Try to change position of the s₁ pole
        A₁ = SE.A[s₁]
        𝕊₁ = SE.𝕊[s₁]
        P₁ = SE.P[s₁]
        P₃ = P₁
        while P₃ == P₁ || sign(P₃) != sign(𝕊₁)
            P₃ = rand(MC.rng, SC.allow)
        end
        P₃ = abs(P₃)
        #
        # Try to change position of the s₂ pole
        A₂ = SE.A[s₂]
        𝕊₂ = SE.𝕊[s₂]
        P₂ = SE.P[s₂]
        P₄ = P₂
        while P₄ == P₂ || sign(P₄) != sign(𝕊₂)
            P₄ = rand(MC.rng, SC.allow)
        end
        P₄ = abs(P₄)

        # Calculate change of Green's function
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        Λ₃ = view(SC.Λ, :, P₃)
        Λ₄ = view(SC.Λ, :, P₄)
        @. δG = 𝕊₁ * A₁ * (Λ₃ - Λ₁) + 𝕊₂ * A₂ * (Λ₄ - Λ₂)

        # Calculate new Green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²

        # Simulated annealing algorithm
        MC.Ptry = MC.Ptry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.P[s₁] = P₃
            SE.P[s₂] = P₄

            # Update reconstructed Green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ² = χ²

            # Update Monte Carlo counter
            MC.Pacc = MC.Pacc + 1

            # Save optimal solution
            if SC.χ² < SC.χ²ᵥ[t]
                measure(t, SE, SC)
            end
        end

    end
end

"""
    try_move_a(
        t::I64,
        MC::StochPXMC,
        SE::StochPXElement,
        SC::StochPXContext
    )

Change the amplitudes of two randomly selected poles.

### Arguments
* t -> Counter of attempts.
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
N/A

See also: [`try_move_x`](@ref).
"""
function try_move_a(
    t::I64,
    MC::StochPXMC,
    SE::StochPXElement,
    SC::StochPXContext
    )
    # Get parameters
    ngrid = length(SC.Gᵥ) # get_b("ngrid")
    npole = get_x("npole")

    # Sanity check
    if npole == 1
        return
    end

    # It is used to save the change of Green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

        # Select two poles randomly
        # The two poles should not be the same.
        s₁ = 1
        s₂ = 1
        while s₁ == s₂
            s₁ = rand(MC.rng, 1:npole)
            s₂ = rand(MC.rng, 1:npole)
        end

        # Try to change amplitudes of the two poles, but their sum is kept.
        P₁ = SE.P[s₁]
        P₂ = SE.P[s₂]
        A₁ = SE.A[s₁]
        A₂ = SE.A[s₂]
        A₃ = 0.0
        A₄ = 0.0
        𝕊₁ = SE.𝕊[s₁]
        𝕊₂ = SE.𝕊[s₂]

        if 𝕊₁ == 𝕊₂
            while true
                δA = rand(MC.rng) * (A₁ + A₂) - A₁
                A₃ = A₁ + δA
                A₄ = A₂ - δA

                if 1.0 > A₃ > 0.0 && 1.0 > A₄ > 0.0
                    break
                end
            end
        else
            while true
                _δA = rand(MC.rng) * (A₁ + A₂) - A₁
                δA = rand(MC.rng) > 0.5 ? _δA * (+1.0) : _δA * (-1.0)
                A₃ = (𝕊₁ * A₁ + δA) / 𝕊₁
                A₄ = (𝕊₂ * A₂ - δA) / 𝕊₂

                if 1.0 > A₃ > 0.0 && 1.0 > A₄ > 0.0
                    break
                end
            end
        end

        # Calculate change of Green's function
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        @. δG = 𝕊₁ * (A₃ - A₁) * Λ₁ + 𝕊₂ * (A₄ - A₂) * Λ₂

        # Calculate new Green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²

        # Simulated annealing algorithm
        MC.Atry = MC.Atry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.A[s₁] = A₃
            SE.A[s₂] = A₄

            # Update reconstructed Green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ² = χ²

            # Update Monte Carlo counter
            MC.Aacc = MC.Aacc + 1

            # Save optimal solution
            if SC.χ² < SC.χ²ᵥ[t]
                measure(t, SE, SC)
            end
        end

    end
end

"""
    try_move_x(
        t::I64,
        MC::StochPXMC,
        SE::StochPXElement,
        SC::StochPXContext
    )

Exchange the amplitudes of two randomly selected poles.

### Arguments
* t -> Counter of attempts.
* MC -> A StochPXMC struct.
* SE -> A StochPXElement struct.
* SC -> A StochPXContext struct.

### Returns
N/A

See also: [`try_move_a`](@ref).
"""
function try_move_x(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Get parameters
    ngrid = length(SC.Gᵥ) # get_b("ngrid")
    offdiag = get_b("offdiag")
    npole = get_x("npole")

    # Sanity check
    if offdiag
        if npole ≤ 3
            return
        end
    else
        if npole == 1
            return
        end
    end

    # It is used to save the change of Green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

        # Select two poles randomly
        # The positions of the two poles are different,
        # but their signs should be the same.
        s₁ = 1
        s₂ = 1
        while ( s₁ == s₂ ) || ( SE.𝕊[s₁] != SE.𝕊[s₂] )
            s₁ = rand(MC.rng, 1:npole)
            s₂ = rand(MC.rng, 1:npole)
        end

        # Try to swap amplitudes of the two poles, but their sum is kept.
        P₁ = SE.P[s₁]
        P₂ = SE.P[s₂]
        A₁ = SE.A[s₁]
        A₂ = SE.A[s₂]
        A₃ = A₂
        A₄ = A₁
        𝕊₁ = SE.𝕊[s₁]
        𝕊₂ = SE.𝕊[s₂]

        # Calculate change of Green's function
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        @. δG = 𝕊₁ * (A₃ - A₁) * Λ₁ + 𝕊₂ * (A₄ - A₂) * Λ₂

        # Calculate new Green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²

        # Simulated annealing algorithm
        MC.Xtry = MC.Xtry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.A[s₁] = A₃
            SE.A[s₂] = A₄

            # Update reconstructed Green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ² = χ²

            # Update Monte Carlo counter
            MC.Xacc = MC.Xacc + 1

            # Save optimal solution
            if SC.χ² < SC.χ²ᵥ[t]
                measure(t, SE, SC)
            end
        end

    end
end
