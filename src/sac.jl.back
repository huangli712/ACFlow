#
# Project : Gardenia
# Source  : sac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/01/21
#

#=
### *Customized Structs* : *StochAC Solver*
=#

"""
    StochACElement

Mutable struct. It is used to record the field configurations, which will
be sampled by  Monte Carlo sweeping procedure.

### Members

* Γₚ -> It means the positions of the δ functions.
* Γₐ -> It means the weights / amplitudes of the δ functions.
"""
mutable struct StochACElement
    Γₚ :: Array{I64,2}
    Γₐ :: Array{F64,2}
end

"""
    StochACContext

Mutable struct. It is used within the StochAC solver only.

### Members

* Gᵥ     -> Input data for correlator.
* σ¹     -> Actually 1.0 / σ¹.
* allow  -> Allowable indices.
* grid   -> Grid for input data.
* mesh   -> Mesh for output spectrum.
* model  -> Default model function.
* kernel -> Default kernel function.
* Aout   -> Calculated spectral function, it is actually ⟨n(x)⟩.
* Δ      -> Precomputed δ functions.
* hτ     -> α-resolved h(τ).
* Hα     -> α-resolved Hc.
* Uα     -> α-resolved internal energy, it is actually ⟨Hα⟩.
* αₗ     -> Vector of the α parameters.
"""
mutable struct StochACContext
    Gᵥ     :: Vector{F64}
    σ¹     :: Vector{F64}
    allow  :: Vector{I64}
    grid   :: AbstractGrid
    mesh   :: AbstractMesh
    model  :: Vector{F64}
    kernel :: Array{F64,2}
    Aout   :: Array{F64,2}
    Δ      :: Array{F64,2}
    hτ     :: Array{F64,2}
    Hα     :: Vector{F64}
    Uα     :: Vector{F64}
    αₗ     :: Vector{F64}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochACSolver, rd::RawData)

Solve the analytical continuation problem by the stochastic analytical
continuation algorithm (K. S. D. Beach's version).
"""
function solve(S::StochACSolver, rd::RawData)
    nmesh = get_b("nmesh")
    nalph = get_a("nalph")

    println("[ StochAC ]")
    MC, SE, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
        println("Using $(nworkers()) workers")
        #
        # Copy configuration dicts
        p1 = deepcopy(PBASE)
        p2 = deepcopy(PStochAC)
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
        Aout = zeros(F64, nmesh, nalph)
        Uα = zeros(F64, nalph)
        for i in eachindex(sol)
            a, b = sol[i]
            @. Aout = Aout + a / nworkers()
            @. Uα = Uα + b / nworkers()
        end
        #
        # Postprocess the solutions
        Gout = last(SC, Aout, Uα)

    # Sequential version
    else
        Aout, Uα = run(MC, SE, SC)
        Gout = last(SC, Aout, Uα)

    end

    return SC.mesh.mesh, Aout, Gout
end

"""
    init(S::StochACSolver, rd::RawData)

Initialize the StochAC solver and return the StochACMC, StochACElement,
and StochACContext structs.
"""
function init(S::StochACSolver, rd::RawData)
    # Initialize possible constraints. The array allow contains all the
    # possible indices for δ functions.
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

    # Only flat model is valid for the StochAC solver.
    model = make_model(mesh)
    println("Build default model: ", get_b("mtype"))

    fmesh = calc_fmesh(S)
    kernel = make_kernel(fmesh, grid)
    println("Build default kernel: ", get_b("ktype"))

    xmesh = calc_xmesh()
    ϕ = calc_phi(mesh, model)
    Δ = calc_delta(xmesh, ϕ)
    println("Precompute δ functions")

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
    U, V, S = make_singular_space(Diagonal(σ¹) * kernel)
    Gᵥ = U' *  (Gᵥ .* σ¹)
    kernel = Diagonal(S) * V'
    hτ, Hα, Uα = calc_hamil(SE.Γₚ, SE.Γₐ, kernel, Gᵥ)
    println("Precompute hamiltonian")

    αₗ = calc_alpha()
    println("Precompute α parameters")

    SC = StochACContext(Gᵥ, σ¹, allow, grid, mesh, model,
                        kernel, Aout, Δ, hτ, Hα, Uα, αₗ)

    return MC, SE, SC
end

"""
    run(MC::StochACMC, SE::StochACElement, SC::StochACContext)

Perform stochastic analytical continuation simulation, sequential version.
"""
function run(MC::StochACMC, SE::StochACElement, SC::StochACContext)
    # Setup essential parameters
    nstep = get_a("nstep")
    output_per_steps = get_a("ndump")
    measure_per_steps = 100

    # Warmup the Monte Carlo engine
    println("Start thermalization...")
    warmup(MC, SE, SC)

    # Sample and collect data
    step = 0.0
    println("Start stochastic sampling...")
    for iter = 1:nstep
        sample(MC, SE, SC)

        if iter % measure_per_steps == 0
            step = step + 1.0
            measure(SE, SC)
        end

        if iter % output_per_steps == 0
            prog = round(I64, iter / nstep * 100)
            @printf("step = %9i ", iter)
            @printf("[progress = %3i]\n", prog)
            flush(stdout)
            write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    prun(S::StochACSolver,
         p1::Dict{String,Vector{Any}},
         p2::Dict{String,Vector{Any}},
         MC::StochACMC, SE::StochACElement, SC::StochACContext)

Perform stochastic analytical continuation simulation, parallel version.
The arguments `p1` and `p2` are copies of PBASE and PStochAC, respectively.
"""
function prun(S::StochACSolver,
              p1::Dict{String,Vector{Any}},
              p2::Dict{String,Vector{Any}},
              MC::StochACMC, SE::StochACElement, SC::StochACContext)
    # Revise parameteric dicts
    rev_dict(p1)
    rev_dict(S, p2)

    # Initialize random number generator again
    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    # Setup essential parameters
    nstep = get_a("nstep")
    output_per_steps = get_a("ndump")
    measure_per_steps = 100

    # Warmup the Monte Carlo engine
    println("Start thermalization...")
    warmup(MC, SE, SC)

    # Sample and collect data
    step = 0.0
    println("Start stochastic sampling...")
    for iter = 1:nstep
        sample(MC, SE, SC)

        if iter % measure_per_steps == 0
            step = step + 1.0
            measure(SE, SC)
        end

        if iter % output_per_steps == 0
            prog = round(I64, iter / nstep * 100)
            @printf("step = %9i ", iter)
            @printf("[progress = %3i]\n", prog)
            flush(stdout)
            myid() == 2 && write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    average(step::F64, SC::StochACContext)

Postprocess the results generated during the stochastic analytical
continuation simulations. It will calculate the spectral functions, and
internal energies.
"""
function average(step::F64, SC::StochACContext)
    # Get key parameters
    nmesh = length(SC.mesh)
    nalph = length(SC.αₗ)

    # Renormalize the spectral functions
    Aout = zeros(F64, nmesh, nalph)
    for i = 1:nalph
        for j = 1:nmesh
            Aout[j,i] = SC.Aout[j,i] * SC.model[j] / π / step
        end
    end

    # Renormalize the internal energies
    Uα = SC.Uα / step

    return Aout, Uα
end

"""
    last(SC::StochACContext, Aout::Array{F64,2}, Uα::Vector{F64})

It will process and write the calculated results by the StochAC solver,
including effective hamiltonian, final spectral function, reproduced
correlator.
"""
function last(SC::StochACContext, Aout::Array{F64,2}, Uα::Vector{F64})
    function fitfun(x, p)
        return @. p[1] * x + p[2]
    end

    # Get dimensional parameters
    nmesh, nalph = size(Aout)

    # Try to fit the internal energies to find out optimal α
    guess = [1.0, 1.0]
    fit_l = curve_fit(fitfun, SC.αₗ[1:5], log10.(Uα[1:5]), guess)
    fit_r = curve_fit(fitfun, SC.αₗ[end-4:end], log10.(Uα[end-4:end]), guess)
    a, b = fit_l.param
    c, d = fit_r.param
    aopt = (d - b) / (a - c)
    close = argmin( abs.( SC.αₗ .- aopt ) )
    println("Fitting parameters [a,b] are: [ $a, $b ]")
    println("Fitting parameters [c,d] are: [ $c, $d ]")
    println("Perhaps the optimal α is: ", aopt)
    write_hamiltonian(SC.αₗ, Uα)

    # Calculate final spectral function and write them
    Asum = zeros(F64, nmesh)
    for i = close : nalph - 1
        @. Asum = Asum + (Uα[i] - Uα[i+1]) * Aout[:,i]
    end
    @. Asum = Asum / (Uα[close] - Uα[end])
    write_spectrum(SC.mesh, Asum)
    write_spectrum(SC.mesh, SC.αₗ, Aout)
    write_model(SC.mesh, SC.model)

    # Reproduce input data and write them
    kernel = make_kernel(SC.mesh, SC.grid)
    G = reprod(SC.mesh, kernel, Asum)
    write_backward(SC.grid, G)

    # Calculate full response function on real axis and write them
    if get_b("ktype") == "fermi"
        _G = kramers(SC.mesh, Asum)
    else
        _G = kramers(SC.mesh, Asum .* SC.mesh)
    end
    write_complete(SC.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

"""
    warmup(MC::StochACMC, SE::StochACElement, SC::StochACContext)

Warmup the Monte Carlo engine to acheieve thermalized equilibrium.
"""
function warmup(MC::StochACMC, SE::StochACElement, SC::StochACContext)
    # Set essential parameter
    nwarm = get_a("nwarm")
    output_per_steps = 100

    # Shuffle the Monte Carlo field configuration
    for iter = 1:nwarm
        sample(MC, SE, SC)

        if iter % output_per_steps == 0
            prog = round(I64, iter / nwarm * 100)
            @printf("step = %9i ", iter)
            @printf("swap( 1 ) -> %8.4f ", MC.Macc[1] / MC.Mtry[1])
            @printf("swap(end) -> %8.4f ", MC.Macc[end] / MC.Mtry[end])
            @printf("[progress = %3i]\n", prog)
            flush(stdout)
        end
    end

    # Reset the counters
    fill!(MC.Macc, 0.0)
    fill!(MC.Mtry, 0.0)

    fill!(MC.Sacc, 0.0)
    fill!(MC.Stry, 0.0)
end

"""
    sample(MC::StochACMC, SE::StochACElement, SC::StochACContext)

Perform Monte Carlo sweeps and sample the field configurations.
"""
function sample(MC::StochACMC, SE::StochACElement, SC::StochACContext)
    nalph = get_a("nalph")

    if rand(MC.rng) < 0.9
        if rand(MC.rng) > 0.5
            for i = 1:nalph
                try_move_a(i, MC, SE, SC)
            end
        else
            for i = 1:nalph
                try_move_p(i, MC, SE, SC)
            end
        end
    else
        if nalph > 1
            try_move_x(MC, SE, SC)
        end
    end
end

"""
    measure(SE::StochACElement, SC::StochACContext)

Measure the spectral functions and internal energies.
"""
function measure(SE::StochACElement, SC::StochACContext)
    nalph = get_a("nalph")

    for ia = 1:nalph
        da = view(SE.Γₐ, :, ia)
        dp = view(SE.Γₚ, :, ia)
        SC.Aout[:,ia] = SC.Aout[:,ia] .+ SC.Δ[:,dp] * da
        SC.Uα[ia] = SC.Uα[ia] + SC.Hα[ia]
    end
end

#=
### *Service Functions*
=#

"""
    init_mc(S::StochACSolver)

Try to create a StochACMC struct.

See also: [`StochACMC`](@ref).
"""
function init_mc(S::StochACSolver)
    nalph = get_a("nalph")

    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    Macc = zeros(F64, nalph)
    Mtry = zeros(F64, nalph)
    Sacc = zeros(F64, nalph)
    Stry = zeros(F64, nalph)

    MC = StochACMC(rng, Macc, Mtry, Sacc, Stry)

    return MC
end

"""
    init_element(S::StochACSolver, rng::AbstractRNG, allow::Vector{I64})

Randomize the configurations for future Monte Carlo sampling. It will
return a StochACElement object.

See also: [`StochACElement`](@ref).
"""
function init_element(S::StochACSolver, rng::AbstractRNG, allow::Vector{I64})
    nalph = get_a("nalph")
    ngamm = get_a("ngamm")

    Γₚ = rand(rng, allow, (ngamm, nalph))
    Γₐ = rand(rng, F64, (ngamm, nalph))

    for j = 1:nalph
        Γⱼ = view(Γₐ, :, j)
        s = sum(Γⱼ)
        @. Γⱼ = Γⱼ / s
    end

    SE = StochACElement(Γₚ, Γₐ)

    return SE
end

"""
    init_iodata(S::StochACSolver, rd::RawData)

Preprocess the input data (`rd`), then allocate memory for the α-resolved
spectral functions.

See also: [`RawData`](@ref).
"""
function init_iodata(S::StochACSolver, rd::RawData)
    nalph = get_a("nalph")
    nmesh = get_b("nmesh")

    Aout = zeros(F64, nmesh, nalph)

    G = make_data(rd)
    Gᵥ = G.value # Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹, Aout
end

"""
    calc_fmesh(S::StochACSolver)

Try to calculate very fine (dense) linear mesh in [wmin, wmax], which
is used internally to build the kernel function.

See also: [`LinearMesh`](@ref).
"""
function calc_fmesh(S::StochACSolver)
    nfine = get_a("nfine")
    wmin = get_b("wmin")
    wmax = get_b("wmax")

    fmesh = LinearMesh(nfine, wmin, wmax)

    return fmesh
end

"""
    calc_xmesh()

Try to calculate very fine (dense) linear mesh in [0, 1], which is used
internally to build the δ functions.

See also: [`calc_delta`](@ref).
"""
function calc_xmesh()
    nfine = get_a("nfine")

    _mesh = fill(1.0/nfine, nfine)
    xmesh = cumsum(_mesh)

    return xmesh
end

"""
    calc_phi(am::AbstractMesh, model::Vector{F64})

Try to calculate ϕ(ω) function. `am` is the mesh for calculated spectrum,
and `model` means the default model function.

See also: [`AbstractMesh`](@ref), [`calc_delta`](@ref).
"""
function calc_phi(am::AbstractMesh, model::Vector{F64})
    ϕ = cumsum(model .* am.weight)
    return ϕ
end

"""
    calc_delta(xmesh::Vector{F64}, ϕ::Vector{F64})

Precompute the δ functions. `xmesh` is a very dense linear mesh in [0, 1]
and `ϕ` is the ϕ function.

See also: [`calc_xmesh`](@ref), [`calc_phi`](@ref).
"""
function calc_delta(xmesh::Vector{F64}, ϕ::Vector{F64})
    nmesh = length(ϕ)
    nfine = length(xmesh)

    η₁ = 0.001
    η₂ = 0.001 ^ 2.0

    Δ = zeros(F64, nmesh, nfine)
    s = similar(ϕ)
    for i = 1:nfine
        @. s = (ϕ - xmesh[i]) ^ 2.0 + η₂
        @. Δ[:,i] = η₁ / s
    end

    return Δ
end

"""
    calc_hamil(Γₚ, Γₐ, kernel, Gᵥ, σ¹)

Initialize h(τ) and H(α) using Eq.(35) and Eq.(36), respectively. `Γₚ`
and `Γₐ` represent n(x), `kernel` means the kernel function, `Gᵥ` is the
correlator. Note that `kernel` and `Gᵥ` have been rotated into singular
space. Please see comments in `init()` for more details.

See also: [`calc_htau`](@ref).
"""
function calc_hamil(Γₚ::Array{I64,2}, Γₐ::Array{F64,2},
                    kernel::Matrix{F64},
                    Gᵥ::Vector{F64})
    nalph = get_a("nalph")
    ngrid = length(Gᵥ) # It is not equal to get_b("ngrid") any more!

    hτ = zeros(F64, ngrid, nalph)
    Hα = zeros(F64, nalph)
    Uα = zeros(F64, nalph)

    for i = 1:nalph
        hτ[:,i] = calc_htau(Γₚ[:,i], Γₐ[:,i], kernel, Gᵥ)
        Hα[i] = dot(hτ[:,i], hτ[:,i])
    end

    return hτ, Hα, Uα
end

"""
    calc_htau(Γₚ, Γₐ, kernel, Gᵥ)

Try to calculate α-dependent h(τ) via Eq.(36). `Γₚ` and `Γₐ` represent
n(x), `kernel` means the kernel function, `Gᵥ` is the correlator. Note
that `kernel` and `Gᵥ` have been rotated into singular space. Please
see comments in `init()` for more details.

See also: [`calc_hamil`](@ref).
"""
function calc_htau(Γₚ::Vector{I64}, Γₐ::Vector{F64},
                   kernel::Matrix{F64},
                   Gᵥ::Vector{F64})
    hτ = similar(Gᵥ)
    #
    for i in eachindex(Gᵥ)
        hτ[i] = dot(Γₐ, view(kernel, i, Γₚ)) - Gᵥ[i]
    end
    #
    return hτ
end

"""
    calc_alpha()

Generate a list for the α parameters
"""
function calc_alpha()
    nalph = get_a("nalph")
    alpha = get_a("alpha")
    ratio = get_a("ratio")

    αₗ = collect(alpha * (ratio ^ (x - 1)) for x in 1:nalph)

    return αₗ
end

"""
    constraints(S::StochACSolver)

Try to implement the constrained stochastic analytical continuation
method. This function will return a collection. It contains all the
allowable indices.

See also: [`StochACSolver`](@ref).
"""
function constraints(S::StochACSolver)
    exclude = get_b("exclude")
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nfine = get_a("nfine")

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
    try_move_a(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)

Select two δ functions randomly and then change their weights. Here `i`
means the index for α parameters.

See also: [`try_move_p`](@ref).
"""
function try_move_a(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)
    # Get current number of δ functions
    ngamm = get_a("ngamm")

    # Choose two δ functions, they are labelled as γ₁ and γ₂, respectively.
    γ₁ = 1
    γ₂ = 1
    while γ₁ == γ₂
        γ₁ = rand(MC.rng, 1:ngamm)
        γ₂ = rand(MC.rng, 1:ngamm)
    end

    # Extract weights for the two δ functions (a₃ and a₄), then try to
    # calculate new weights for them (a₁ and a₂).
    a₁ = 0.0
    a₂ = 0.0
    a₃ = SE.Γₐ[γ₁,i]
    a₄ = SE.Γₐ[γ₂,i]
    δa = 0.0
    while true
        δa = rand(MC.rng) * (a₃ + a₄) - a₃
        a₁ = a₃ + δa
        a₂ = a₄ - δa
        if a₁ > 0 && a₂ > 0
            break
        end
    end

    # Try to calculate the change of Hc using Eq.~(42).
    hc = view(SC.hτ, :, i)
    K₁ = view(SC.kernel, :, SE.Γₚ[γ₁,i])
    K₂ = view(SC.kernel, :, SE.Γₚ[γ₂,i])
    #
    δhc = δa * (K₁ - K₂)
    δH = dot(δhc, 2.0 * hc + δhc)

    # Apply Metropolis algorithm
    MC.Mtry[i] = MC.Mtry[i] + 1
    if δH ≤ 0.0 || exp(-SC.αₗ[i] * δH) > rand(MC.rng)
        # Update Monte Carlo configurations
        SE.Γₐ[γ₁,i] = a₁
        SE.Γₐ[γ₂,i] = a₂

        # Update h(τ)
        @. hc = hc + δhc

        # Update Hc
        SC.Hα[i] = SC.Hα[i] + δH

        # Update Monte Carlo counter
        MC.Macc[i] = MC.Macc[i] + 1
    end
end

"""
    try_move_p(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)

Select two δ functions randomly and then change their positions. Here `i`
means the index for α parameters.

See also: [`try_move_a`](@ref).
"""
function try_move_p(i::I64, MC::StochACMC, SE::StochACElement, SC::StochACContext)
    # Get current number of δ functions
    ngamm = get_a("ngamm")

    # Choose two δ functions, they are labelled as γ₁ and γ₂, respectively.
    γ₁ = 1
    γ₂ = 1
    while γ₁ == γ₂
        γ₁ = rand(MC.rng, 1:ngamm)
        γ₂ = rand(MC.rng, 1:ngamm)
    end

    # Extract weights for the two δ functions (a₁ and a₂)
    a₁ = SE.Γₐ[γ₁,i]
    a₂ = SE.Γₐ[γ₂,i]

    # Choose new positions for the two δ functions (p₁ and p₂).
    # Note that their old positions are SE.Γₚ[γ₁,i] and SE.Γₚ[γ₂,i].
    p₁ = rand(MC.rng, SC.allow)
    p₂ = rand(MC.rng, SC.allow)

    # Try to calculate the change of Hc using Eq.~(42).
    hc = view(SC.hτ, :, i)
    K₁ = view(SC.kernel, :, p₁)
    K₂ = view(SC.kernel, :, p₂)
    K₃ = view(SC.kernel, :, SE.Γₚ[γ₁,i])
    K₄ = view(SC.kernel, :, SE.Γₚ[γ₂,i])
    #
    δhc = a₁ * (K₁ - K₃) + a₂ * (K₂ - K₄)
    δH = dot(δhc, 2.0 * hc + δhc)

    # Apply Metropolis algorithm
    MC.Mtry[i] = MC.Mtry[i] + 1
    if δH ≤ 0.0 || exp(-SC.αₗ[i] * δH) > rand(MC.rng)
        # Update Monte Carlo configurations
        SE.Γₚ[γ₁,i] = p₁
        SE.Γₚ[γ₂,i] = p₂

        # Update h(τ)
        @. hc = hc + δhc

        # Update Hc
        SC.Hα[i] = SC.Hα[i] + δH

        # Update Monte Carlo counter
        MC.Macc[i] = MC.Macc[i] + 1
    end
end

"""
    try_move_x(MC::StochACMC, SE::StochACElement, SC::StochACContext)

Try to exchange field configurations between two adjacent layers.
"""
function try_move_x(MC::StochACMC, SE::StochACElement, SC::StochACContext)
    # Get number of α parameters
    nalph = get_a("nalph")

    # Select two adjacent layers (two adjacent α parameters)
    i = rand(MC.rng, 1:nalph)
    j = rand(MC.rng) > 0.5 ? i + 1 : i - 1
    i == 1 && (j = i + 1)
    i == nalph && (j = i - 1)

    # Calculate change of Hc
    δα = SC.αₗ[i] - SC.αₗ[j]
    δH = SC.Hα[i] - SC.Hα[j]

    # Apply Metropolis algorithm
    MC.Stry[i] = MC.Stry[i] + 1
    MC.Stry[j] = MC.Stry[j] + 1
    if exp(δα * δH) > rand(MC.rng)
        # Update Monte Carlo configurations
        SE.Γₚ[:,i], SE.Γₚ[:,j] = SE.Γₚ[:,j], SE.Γₚ[:,i]
        SE.Γₐ[:,i], SE.Γₐ[:,j] = SE.Γₐ[:,j], SE.Γₐ[:,i]

        # Update h(τ) and Hc
        SC.hτ[:,i], SC.hτ[:,j] = SC.hτ[:,j], SC.hτ[:,i]
        SC.Hα[i], SC.Hα[j] = SC.Hα[j], SC.Hα[i]

        # Update Monte Carlo counters
        MC.Sacc[i] = MC.Sacc[i] + 1
        MC.Sacc[j] = MC.Sacc[j] + 1
    end
end
