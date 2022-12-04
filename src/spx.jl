#
# Project : Gardenia
# Source  : spx.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/12/05
#

#=
### *Customized Structs* : *StochPX Solver*
=#

"""
    StochPXElement

Mutable struct. It is used to record the field configurations, which will
be sampled by Monte Carlo sweeping procedure.

### Members

* P -> It means the positions of the poles.
* A -> It means the weights / amplitudes of the poles.
"""
mutable struct StochPXElement
    P :: Vector{I64}
    A :: Vector{F64}
end

"""
    StochPXContext

Mutable struct. It is used within the StochPX solver only.

### Members

* Gᵥ     -> Input data for correlator.
* Gᵧ     -> Generated correlator.
* σ¹     -> Actually 1.0 / σ¹.
* allow  -> Allowable indices.
* grid   -> Grid for input data.
* mesh   -> Mesh for output spectrum.
* fmesh  -> Very dense mesh for the poles.
* χ²     -> Vector of goodness function.
* Pᵥ     -> Vector of poles' positions.
* Aᵥ     -> Vector of poles' amplitudes.
"""
mutable struct StochPXContext
    Gᵥ    :: Vector{F64}
    Gᵧ    :: Vector{F64}
    σ¹    :: Vector{F64}
    allow :: Vector{I64}
    grid  :: AbstractGrid
    mesh  :: AbstractMesh
    fmesh :: AbstractMesh
    χ²    :: Vector{F64}
    Pᵥ    :: Vector{Vector{I64}}
    Aᵥ    :: Vector{Vector{F64}}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochPXSolver, rd::RawData)

Solve the analytical continuation problem by the stochastic
pole expansion. Note that this solver is still `experimental`.
"""
function solve(S::StochPXSolver, rd::RawData)
    ngrid = get_b("ngrid")
    nmesh = get_b("nmesh")

    println("[ StochPX ]")
    MC, SE, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
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
        Gᵣ = zeros(F64, 2 * ngrid)
        for i in eachindex(sol)
            a, b, c = sol[i]
            @. Aout = Aout + a / nworkers()
            @. Gout = Gout + b / nworkers()
            @. Gᵣ = Gᵣ + c / nworkers()
        end
        #
        # Postprocess the solutions
        last(SC, Aout, Gout, Gᵣ)

    # Sequential version
    else
        Aout, Gout, Gᵣ = run(MC, SE, SC)
        last(SC, Aout, Gout, Gᵣ)

    end

    return SC.mesh.mesh, Aout, Gout
end

"""
    init(S::StochPXSolver, rd::RawData)

Initialize the StochPX solver and return the StochPXMC, StochPXElement,
and StochPXContext structs.
"""
function init(S::StochPXSolver, rd::RawData)
    # Initialize possible constraints. The allow array contains all the
    # possible indices for poles.
    allow = constraints(S)

    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    SE = init_element(S, MC.rng, allow)
    println("Randomize Monte Carlo configurations")

    Gᵥ, σ¹ = init_iodata(S, rd)
    Gᵧ = deepcopy(Gᵥ)
    println("Postprocess input data: ", length(σ¹), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    fmesh = calc_fmesh(S)
    println("Build mesh for spectrum: ", length(mesh), " points")

    χ², Pᵥ, Aᵥ = init_context(S)
    SC = StochPXContext(Gᵥ, Gᵧ, σ¹, allow, grid, mesh, fmesh, χ², Pᵥ, Aᵥ)

    return MC, SE, SC
end

"""
    run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Perform stochastic pole expansion simulation, sequential version.
"""
function run(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Setup essential parameters
    ntry = get_x("ntry")
    nstep = get_x("nstep")
    threshold = 1e-6

    println("Start stochastic sampling...")
    #
    # Sample and collect data
    for t = 1:ntry
        # Reset Monte Carlo counters
        reset_mc(MC)

        # Reset Monte Carlo field configuration
        reset_element(MC.rng, SC.allow, SE)

        # Reset Gᵧ and χ²
        reset_context(t, SE, SC)

        # Apply simulated annealling algorithm
        for i = 1:nstep
            sample(t, MC, SE, SC)

            # Check convergence
            if SC.χ²[t] < threshold
                @printf("try = %6i ", t)
                @printf("[χ² = %9.4e]\n", SC.χ²[t])
                flush(stdout)
                break
            else
                if i == nstep
                    @printf("try = %6i ", t)
                    @printf("[χ² = %9.4e] FAILED\n", SC.χ²[t])
                    flush(stdout)
                end
            end
        end

        # Write Monte Carlo statistics
        write_statistics(MC)

        # Record Monte Carlo field configuration
        measure(t, SE, SC)
    end
    #
    # Print summary information
    passed = count(<(threshold), SC.χ²)
    failed = count(≥(threshold), SC.χ²)
    println("Summary: passed [$passed] failed [$failed]")

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    prun(S::StochPXSolver,
         p1::Dict{String,Vector{Any}},
         p2::Dict{String,Vector{Any}},
         MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Perform stochastic pole expansion simulation, parallel version.
The arguments `p1` and `p2` are copies of PBASE and PStochPX, respectively.
"""
function prun(S::StochPXSolver,
              p1::Dict{String,Vector{Any}},
              p2::Dict{String,Vector{Any}},
              MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Revise parameteric dicts
    rev_dict(p1)
    rev_dict(S, p2)

    # Initialize random number generator again
    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    # Setup essential parameters
    ntry = get_x("ntry")
    nstep = get_x("nstep")
    threshold = 1e-6

    println("Start stochastic sampling...")
    #
    # Sample and collect data
    for t = 1:ntry
        # Reset Monte Carlo counters
        reset_mc(MC)

        # Reset Monte Carlo field configuration
        reset_element(MC.rng, SC.allow, SE)

        # Reset Gᵧ and χ²
        reset_context(t, SE, SC)

        # Apply simulated annealling algorithm
        for i = 1:nstep
            sample(t, MC, SE, SC)

            # Check convergence
            if SC.χ²[t] < threshold
                @printf("try = %6i ", t)
                @printf("[χ² = %9.4e]\n", SC.χ²[t])
                flush(stdout)
                break
            else
                if i == nstep
                    @printf("try = %6i ", t)
                    @printf("[χ² = %9.4e] FAILED\n", SC.χ²[t])
                    flush(stdout)
                end
            end
        end

        # Write Monte Carlo statistics
        myid() == 2 && write_statistics(MC)

        # Record Monte Carlo field configuration
        measure(t, SE, SC)
    end
    #
    # Print summary information
    passed = count(<(threshold), SC.χ²)
    failed = count(≥(threshold), SC.χ²)
    println("Summary: passed [$passed] failed [$failed]")

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    average(SC::StochPXContext)

Postprocess the results generated during the stochastic pole expansion
simulations. It will generate the spectral functions, real frequency
green's function, and imaginary frequency green's function.
"""
function average(SC::StochPXContext)
    # Setup essential parameters
    nmesh = get_b("nmesh")
    ngrid = get_b("ngrid")
    ntry = get_x("ntry")

    # The threshold is used to distinguish good or bad solutions
    threshold = 1e-6

    # Allocate memory
    # Gout: real frequency green's function, G(ω).
    # Gᵣ: imaginary frequency green's function, G(iωₙ)
    Gout = zeros(C64, nmesh)
    Gᵣ = zeros(F64, 2 * ngrid)

    # Go through all the solutions
    c = 0.0
    for i = 1:ntry
        if SC.χ²[i] < threshold
            G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.mesh, SC.fmesh)
            @. Gout = Gout + G
            #
            G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.grid, SC.fmesh)
            @. Gᵣ = Gᵣ + G
            #
            # Increase the counter
            c = c + 1.0
        end
    end
    #
    # Average the final results
    @. Gout = Gout / c
    @. Gᵣ = Gᵣ / c

    return -imag.(Gout) / π, Gout, Gᵣ
end

"""
    last(SC::StochPXContext, Aout::Vector{F64}, Gout::Vector{C64}, Gᵣ::Vector{F64})

It will write the calculated results by the StochPX solver, including
final spectral function and reproduced correlator.
"""
function last(SC::StochPXContext, Aout::Vector{F64}, Gout::Vector{C64}, Gᵣ::Vector{F64})
    # Write the spectral function
    write_spectrum(SC.mesh, Aout)

    # Reproduce input data and write them
    write_backward(SC.grid, Gᵣ)

    # Write full response function on real axis
    write_complete(SC.mesh, Gout)

    # Write pole expansion coefficients
    write_pole(SC.Pᵥ, SC.Aᵥ, SC.fmesh)
end

#=
### *Core Algorithms*
=#

"""
    sample(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Try to search the configuration space to locate the minimum by using the
simulated annealling algorithm. Here, `t` means the t-th attempt.
"""
function sample(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Try to change positions of two poles
    if rand(MC.rng) < 0.5
        try_move_p(t, MC, SE, SC)
    # Try to change amplitudes of two poles
    else
        try_move_a(t, MC, SE, SC)
    end
end

"""
    measure(t::I64, SE::StochPXElement, SC::StochPXContext)

Store Monte Carlo field configurations (positions and amplitudes of many
poles) for the t-th attempt.
"""
function measure(t::I64, SE::StochPXElement, SC::StochPXContext)
    @. SC.Pᵥ[t] = SE.P
    @. SC.Aᵥ[t] = SE.A
end

#=
### *Service Functions*
=#

"""
    init_mc(S::StochPXSolver)

Try to create a StochPXMC struct.

See also: [`StochPXMC`](@ref).
"""
function init_mc(S::StochPXSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    #
    Pacc = 0
    Ptry = 0
    Aacc = 0
    Atry = 0
    Sacc = 0
    Stry = 0

    MC = StochPXMC(rng, Pacc, Ptry, Aacc, Atry, Sacc, Stry)

    return MC
end

"""
    init_element(S::StochPXSolver, rng::AbstractRNG, allow::Vector{I64})

Randomize the configurations for future Monte Carlo sampling. It will
return a StochPXElement object.

See also: [`StochPXElement`](@ref).
"""
function init_element(S::StochPXSolver, rng::AbstractRNG, allow::Vector{I64})
    npole = get_x("npole")

    P = rand(rng, allow, npole)
    A = rand(rng, F64, npole)

    # We have to make sure ∑ Aᵢ = 1
    s = sum(A)
    @. A = A / s

    SE = StochPXElement(P, A)

    return SE
end

"""
    init_iodata(S::StochPXSolver, rd::RawData)

Preprocess the input data (`rd`).

See also: [`RawData`](@ref).
"""
function init_iodata(S::StochPXSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value # Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹
end

"""
    init_context(S::StochPXSolver)

Try to initialize the key members of a StochPXContext struct.

See also: [`StochPXContext`](@ref).
"""
function init_context(S::StochPXSolver)
    ntry = get_x("ntry")
    npole = get_x("npole")

    χ² = zeros(F64, ntry)

    Pᵥ = Vector{I64}[]
    Aᵥ = Vector{F64}[]
    for i = 1:ntry
        push!(Pᵥ, zeros(I64, npole))
        push!(Aᵥ, zeros(F64, npole))
    end

    return χ², Pᵥ, Aᵥ
end

"""
    reset_mc(MC::StochPXMC)

Reset the counters in StochPXMC struct.
"""
function reset_mc(MC::StochPXMC)
    MC.Pacc = 0
    MC.Ptry = 0
    MC.Aacc = 0
    MC.Atry = 0
    MC.Sacc = 0
    MC.Stry = 0
end

"""
    reset_element(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)

Reset the Monte Carlo field configurations (i.e. positions and amplitudes
of the poles).
"""
function reset_element(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
    npole = get_x("npole")

    P = rand(rng, allow, npole)
    A = rand(rng, F64, npole)

    s = sum(A)

    @. SE.P = P
    @. SE.A = A / s
end

"""
    reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)

Recalculate imaginary frequency green's function and goodness-of-fit
function by new Monte Carlo field configurations for the t-th attempts.
"""
function reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)
    Gᵧ = calc_green(SE.P, SE.A, SC.grid, SC.fmesh)
    χ² = calc_chi2(Gᵧ, SC.Gᵥ)
    @. SC.Gᵧ = Gᵧ
    SC.χ²[t] = χ²
end

"""
    calc_fmesh(S::StochPXSolver)

Try to calculate very fine (dense) linear mesh in [wmin, wmax], which
is used internally to represent the possible positions of poles.

See also: [`LinearMesh`](@ref).
"""
function calc_fmesh(S::StochPXSolver)
    nfine = get_x("nfine")
    wmin = get_b("wmin")
    wmax = get_b("wmax")

    fmesh = LinearMesh(nfine, wmin, wmax)

    return fmesh
end

"""
    calc_green(P::Vector{I64},
               A::Vector{F64},
               grid::AbstractGrid,
               fmesh::AbstractMesh)

Reconstruct green's function at imaginary axis by the pole expansion.
"""
function calc_green(P::Vector{I64},
                    A::Vector{F64},
                    grid::AbstractGrid,
                    fmesh::AbstractMesh)
    npole = get_x("npole")
    ngrid = length(grid)

    G = zeros(C64, ngrid)
    for i in eachindex(grid)
        z = 0.0
        for j = 1:npole
            z = z + A[j] / ( im * grid[i] - fmesh[P[j]] )
        end
        G[i] = z
    end

    return vcat(real(G), imag(G))
end

"""
    calc_green(P::Vector{I64},
               A::Vector{F64},
               mesh::AbstractMesh,
               fmesh::AbstractMesh)

Reconstruct green's function at real axis by the pole expansion.
"""
function calc_green(P::Vector{I64},
                    A::Vector{F64},
                    mesh::AbstractMesh,
                    fmesh::AbstractMesh)
    npole = get_x("npole")
    η = get_x("eta")
    nmesh = length(mesh)

    G = zeros(C64, nmesh)
    for i in eachindex(mesh)
        z = 0.0
        for j = 1:npole
            z = z + A[j] / ( mesh[i] + im * η - fmesh[P[j]] )
        end
        G[i] = z
    end

    return G
end

"""
    calc_chi2(Gₙ::Vector{F64}, Gᵥ::Vector{F64})

Try to calculate the goodness function (i.e, χ²), which measures the
distance between input and regenerated correlators.

See also: [`calc_green`](@ref).
"""
function calc_chi2(Gₙ::Vector{F64}, Gᵥ::Vector{F64})
    ΔG = Gₙ - Gᵥ
    return dot(ΔG, ΔG)
end

"""
    constraints(S::StochPXSolver)

Try to implement the constrained stochastic pole expansion. This
function will return a collection. It contains all the allowable indices.

See also: [`StochPXSolver`](@ref).
"""
function constraints(S::StochPXSolver)
    exclude = get_b("exclude")
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nfine = get_x("nfine")

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
    try_move_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Change the positions of two randomly selected poles.

See also: [`try_move_a`](@ref).
"""
function try_move_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Get parameters
    ngrid = get_b("ngrid")
    npole = get_x("npole")

    # It is used to save the change of green's function
    δG = zeros(C64, ngrid)

    # Select two poles randomly
    s₁ = 1
    s₂ = 1
    while s₁ == s₂
        s₁ = rand(MC.rng, 1:npole)
        s₂ = rand(MC.rng, 1:npole)
    end

    # Try to change position of the s₁ pole
    P₁ = SE.P[s₁]
    P₃ = P₁
    while P₃ == P₁
        P₃ = rand(MC.rng, SC.allow)
    end
    A₁ = SE.A[s₁]
    #
    # Try to change position of the s₂ pole
    P₂ = SE.P[s₂]
    P₄ = P₂
    while P₄ == P₂
        P₄ = rand(MC.rng, SC.allow)
    end
    A₂ = SE.A[s₂]

    # Calculate change of green's function
    for i in eachindex(SC.grid)
        z = 0 + A₁ / ( im * SC.grid[i] - SC.fmesh[P₃] )
        z = z - A₁ / ( im * SC.grid[i] - SC.fmesh[P₁] )
        z = z + A₂ / ( im * SC.grid[i] - SC.fmesh[P₄] )
        z = z - A₂ / ( im * SC.grid[i] - SC.fmesh[P₂] )
        δG[i] = z
    end

    # Calculate new green's function and goodness-of-fit function
    Gₙ = SC.Gᵧ + vcat(real(δG), imag(δG))
    χ² = calc_chi2(Gₙ, SC.Gᵥ)

    # Simulated annealling algorithm
    MC.Ptry = MC.Ptry + 1
    if χ² < SC.χ²[t] || rand(MC.rng) < min(exp(-(χ² - SC.χ²[t])), 1.0)
        # Update Monte Carlo configuration
        SE.P[s₁] = P₃
        SE.P[s₂] = P₄

        # Update reconstructed green's function
        @. SC.Gᵧ = Gₙ

        # Update goodness-of-fit function
        SC.χ²[t] = χ²

        # Update Monte Carlo counter
        MC.Pacc = MC.Pacc + 1
    end
end

"""
    try_move_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Change the amplitudes of two randomly selected poles.

See also: [`try_move_p`](@ref).
"""
function try_move_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Get parameters
    ngrid = get_b("ngrid")
    npole = get_x("npole")

    # It is used to save the change of green's function
    δG = zeros(C64, ngrid)

    # Select two poles randomly
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
    while true
        δA = rand(MC.rng) * (A₁ + A₂) - A₁
        A₃ = A₁ + δA
        A₄ = A₂ - δA

        if A₃ > 0 && A₄ > 0
            break
        end
    end

    # Calculate change of green's function
    for i in eachindex(SC.grid)
        z = 0 + (A₃ - A₁) / ( im * SC.grid[i] - SC.fmesh[P₁] )
        z = z + (A₄ - A₂) / ( im * SC.grid[i] - SC.fmesh[P₂] )
        δG[i] = z
    end

    # Calculate new green's function and goodness-of-fit function
    Gₙ = SC.Gᵧ + vcat(real(δG), imag(δG))
    χ² = calc_chi2(Gₙ, SC.Gᵥ)

    # Simulated annealling algorithm
    MC.Atry = MC.Atry + 1
    if χ² < SC.χ²[t] || rand(MC.rng) < min(exp(-(χ² - SC.χ²[t])), 1.0)
        # Update Monte Carlo configuration
        SE.A[s₁] = A₃
        SE.A[s₂] = A₄

        # Update reconstructed green's function
        @. SC.Gᵧ = Gₙ

        # Update goodness-of-fit function
        SC.χ²[t] = χ²

        # Update Monte Carlo counter
        MC.Aacc = MC.Aacc + 1
    end
end
