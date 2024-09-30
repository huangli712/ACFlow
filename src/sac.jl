#
# Project : Gardenia
# Source  : sac.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/09/30
#

#=
### *Customized Structs* : *StochAC Solver*
=#

"""
    StochACElement

Mutable struct. It is used to record the field configurations, which will
be sampled by Monte Carlo sweeping procedure.

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
* grid   -> Imaginary axis grid for input data.
* mesh   -> Real frequency mesh for output spectrum.
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

Solve the analytic continuation problem by the stochastic analytic
continuation algorithm (K. S. D. Beach's version). This is the driver for
the StochAC solver.

If the input correlators are bosonic, this solver will return A(ω) / ω
via `Asum`, instead of A(ω). At this time, `Asum` is not compatible with
`Gout`. If the input correlators are fermionic, this solver will return
A(ω) in `Asum`. Now it is compatible with `Gout`. These behaviors are just
similar to the MaxEnt, StochSK, and StochOM solvers.

Now the StochAC solver supports both continuous and δ-like spectra.

### Arguments
* S -> A StochACSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* mesh -> Real frequency mesh, ω.
* Asum -> Final spectral function, A(ω). Note that it is α-averaged.
* Gout -> Retarded Green's function, G(ω).
"""
function solve(S::StochACSolver, rd::RawData)
    nmesh = get_b("nmesh")
    nalph = get_a("nalph")

    println("[ StochAC ]")
    MC, SE, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
        #
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
        Asum, Gout = last(SC, Aout, Uα)
        #
    # Sequential version
    else
        #
        Aout, Uα = run(MC, SE, SC)
        Asum, Gout = last(SC, Aout, Uα)
        #
    end

    return SC.mesh.mesh, Asum, Gout
end

"""
    init(S::StochACSolver, rd::RawData)

Initialize the StochAC solver and return the StochACMC, StochACElement,
and StochACContext structs. Please don't call this function directly.

### Arguments
* S -> A StochACSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.
"""
function init(S::StochACSolver, rd::RawData)
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
    println("Initialize context for the StochAC solver")

    return MC, SE, SC
end

"""
    run(MC::StochACMC, SE::StochACElement, SC::StochACContext)

Perform stochastic analytic continuation simulation, sequential version.

### Arguments
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* Uα -> α-resolved internal energy.
"""
function run(MC::StochACMC, SE::StochACElement, SC::StochACContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

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
            fwrite && write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    prun(
        S::StochACSolver,
        p1::Dict{String,Vector{Any}},
        p2::Dict{String,Vector{Any}},
        MC::StochACMC,
        SE::StochACElement,
        SC::StochACContext
    )

Perform stochastic analytic continuation simulation, parallel version.
The arguments `p1` and `p2` are copies of PBASE and PStochAC, respectively.

### Arguments
* S -> A StochACSolver struct.
* p1 -> A copy of PBASE.
* p2 -> A copy of PStochAC.
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
* Aout -> Spectral function, A(ω).
* Uα -> α-resolved internal energy.
"""
function prun(
    S::StochACSolver,
    p1::Dict{String,Vector{Any}},
    p2::Dict{String,Vector{Any}},
    MC::StochACMC,
    SE::StochACElement,
    SC::StochACContext
    )
    # Revise parameteric dicts
    rev_dict_b(p1)
    rev_dict_a(S, p2)

    # Initialize random number generator again
    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

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
            myid() == 2 && fwrite && write_statistics(MC)
        end
    end

    return average(step, SC)
end

"""
    average(step::F64, SC::StochACContext)

Postprocess the results generated during the stochastic analytic
continuation simulations. It will calculate the spectral functions, and
α-resolved internal energies.

### Arguments
* step -> How many steps are there in the Monte Carlo samplings.
* SC   -> A StochACContext struct.

### Returns
* Aout -> Spectral function, A(ω,α).
* Uα -> α-resolved internal energy.
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

### Arguments
* SC   -> A StochACContext struct.
* Aout -> α-dependent spectral functions.
* Uα   -> α-dependent internal energies.

### Returns
* Asum -> Final spectral function (α-averaged), A(ω).
* G -> Retarded Green's function, G(ω).
"""
function last(SC::StochACContext, Aout::Array{F64,2}, Uα::Vector{F64})
    function fitfun(x, p)
        return @. p[1] * x + p[2]
    end

    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

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
    fwrite && write_hamiltonian(SC.αₗ, Uα)

    # Calculate final spectral functions and write them
    Asum = zeros(F64, nmesh)
    for i = close : nalph - 1
        @. Asum = Asum + (Uα[i] - Uα[i+1]) * Aout[:,i]
    end
    @. Asum = Asum / (Uα[close] - Uα[end])
    fwrite && write_spectrum(SC.mesh, Asum)
    fwrite && write_spectrum(SC.mesh, SC.αₗ, Aout)
    fwrite && write_model(SC.mesh, SC.model)

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

    return Asum, _G
end

#=
### *Core Algorithms*
=#

"""
    warmup(MC::StochACMC, SE::StochACElement, SC::StochACContext)

Warmup the Monte Carlo engine to acheieve thermalized equilibrium. After
that, the Monte Carlo counters will be reset.

### Arguments
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
N/A
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

### Arguments
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
N/A
"""
function sample(MC::StochACMC, SE::StochACElement, SC::StochACContext)
    nalph = get_a("nalph")

    if rand(MC.rng) < 0.9
        if rand(MC.rng) > 0.5
            for i = 1:nalph
                try_move_a(i, MC, SE, SC)
            end
        else
            if rand(MC.rng) > 0.2
                for i = 1:nalph
                    try_move_s(i, MC, SE, SC)
                end
            else
                for i = 1:nalph
                    try_move_p(i, MC, SE, SC)
                end
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

Accumulate the α-resolved spectral functions and internal energies.

### Arguments
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
N/A
"""
function measure(SE::StochACElement, SC::StochACContext)
    nalph = get_a("nalph")

    # Loop over each α parameter
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
    init_iodata(S::StochACSolver, rd::RawData)

Preprocess the input data (`rd`).

### Arguments
* S -> A StochACSolver struct.
* rd -> A RawData struct, which contains essential input data.

### Returns
* Gᵥ -> Input correlator.
* σ¹ -> 1.0 / σ¹.

See also: [`RawData`](@ref).
"""
function init_iodata(S::StochACSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value # Gᵥ = abs.(G.value)
    σ¹ = 1.0 ./ sqrt.(G.covar)

    return Gᵥ, σ¹
end

"""
    init_mc(S::StochACSolver)

Try to create a StochACMC struct. Some counters for Monte Carlo updates
are initialized here.

### Arguments
* S -> A StochACSolver struct.

### Returns
* MC -> A StochACMC struct.

See also: [`StochACMC`](@ref).
"""
function init_mc(S::StochACSolver)
    nalph = get_a("nalph")
    #
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    #
    Macc = zeros(F64, nalph)
    Mtry = zeros(F64, nalph)
    Sacc = zeros(F64, nalph)
    Stry = zeros(F64, nalph)
    #
    MC = StochACMC(rng, Macc, Mtry, Sacc, Stry)

    return MC
end

"""
    init_element(
        S::StochACSolver,
        rng::AbstractRNG,
        allow::Vector{I64}
    )

Randomize the configurations for future Monte Carlo sampling. It will
return a StochACElement struct.

### Arguments
* S     -> A StochACSolver struct.
* rng   -> Random number generator.
* allow -> Allowed positions for the δ peaks.

### Returns
* SE -> A StochACElement struct.

See also: [`StochACElement`](@ref).
"""
function init_element(
    S::StochACSolver,
    rng::AbstractRNG,
    allow::Vector{I64}
    )
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
    init_context(
        SE::StochACElement,
        Gᵥ::Vector{F64},
        σ¹::Vector{F64},
        allow::Vector{I64},
        grid::AbstractGrid,
        mesh::AbstractMesh,
        fmesh::AbstractMesh
    )

Try to create a StochACContext struct, which contains some key variables,
including grid, mesh, input correlator and the corresponding standard
deviation, kernel matrix, spectral function, and α-resolved Hamiltonian.

### Arguments
* SE -> A StochACElement struct.
* Gᵥ -> Input correlator. It will be changed in this function.
* σ¹ -> Standard deviation for input correlator.
* allow -> Allowable indices for δ-like peaks.
* grid -> Imaginary axis grid for input data.
* mesh -> Real frequency mesh for output spectrum.
* fmesh -> Very fine mesh in [wmin, wmax].

### Returns
* SC -> A StochACContext struct.
"""
function init_context(
    SE::StochACElement,
    Gᵥ::Vector{F64},
    σ¹::Vector{F64},
    allow::Vector{I64},
    grid::AbstractGrid,
    mesh::AbstractMesh,
    fmesh::AbstractMesh
    )
    # Get parameters
    nmesh = get_b("nmesh")
    nalph = get_a("nalph")

    # Allocate memory for spectral function, A(ω,α)
    Aout = zeros(F64, nmesh, nalph)

    # Prepare some key variables
    # Only flat model is valid for the StochAC solver.
    model = make_model(mesh)

    # Precompute δ functions
    ϕ = calc_phi(mesh, model)
    Δ = calc_delta(fmesh, ϕ)

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

    # Precompute hamiltonian
    hτ, Hα, Uα = calc_hamil(SE.Γₚ, SE.Γₐ, kernel, Gᵥ)

    # Precompute α parameters
    αₗ = calc_alpha()

    return StochACContext(Gᵥ, σ¹, allow, grid, mesh, model,
                        kernel, Aout, Δ, hτ, Hα, Uα, αₗ)
end

"""
    calc_fmesh(S::StochACSolver)

Try to calculate very fine (dense) mesh in [wmin, wmax], which is used
internally to build the kernel function. Note that this mesh could be
non-uniform. If the file `fmesh.inp` exists, the code will try to load
it to initialize the mesh. Or else the code will try to create a linear
mesh automatically.

### Arguments
* S -> A StochACSolver struct.

### Returns
* fmesh -> A very fine, perhaps non-uniform mesh in [wmin, wmax].

See also: [`LinearMesh`](@ref), [`DynamicMesh`](@ref).
"""
function calc_fmesh(S::StochACSolver)
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nfine = get_a("nfine")

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

"""
    calc_phi(am::AbstractMesh, model::Vector{F64})

Try to calculate ϕ(ω) function. `am` is the mesh for calculated spectrum,
and `model` means the default model function.

For the definition of the ϕ(ω) function, see Eq.(16) in arXiv:0403055. It
creates a smooth mapping from ℝ to [0,1].

### Arguments
See above explanations.

### Returns
* ϕ -> The ϕ(ω) function.

See also: [`calc_delta`](@ref).
"""
function calc_phi(am::AbstractMesh, model::Vector{F64})
    ϕ = cumsum(model .* am.weight)
    return ϕ
end

"""
    calc_delta(fmesh::AbstractMesh, ϕ::Vector{F64})

Precompute the Δ functions. `fmesh` is a very dense mesh in [wmin, wmax]
and `ϕ` is the ϕ function.

Here we just use f(x) = η / (x² + η²) to approximate the δ function, where
η is a small parameter.

### Arguments
See above explanations.

### Returns
* Δ -> The Δ(ω) function.

See also: [`calc_phi`](@ref).
"""
function calc_delta(fmesh::AbstractMesh, ϕ::Vector{F64})
    nmesh = length(ϕ)
    #
    nfine = length(fmesh)
    wmax = fmesh.wmax
    wmin = fmesh.wmin
    #
    η₁ = 0.001
    η₂ = 0.001 ^ 2.0

    Δ = zeros(F64, nmesh, nfine)
    s = similar(ϕ)
    for i = 1:nfine
        # We should convert the mesh `fmesh` from [wmin,wmax] to [0,1].
        𝑥 = (fmesh[i] - wmin) / (wmax - wmin)
        @. s = (ϕ - 𝑥) ^ 2.0 + η₂
        @. Δ[:,i] = η₁ / s
    end

    return Δ
end

"""
    calc_hamil(
        Γₚ::Array{I64,2},
        Γₐ::Array{I64,2},
        kernel::Matrix{F64},
        Gᵥ::Vector{F64}
    )

Initialize h(τ) and H(α) using Eq.(35) and Eq.(36), respectively. `Γₚ`
and `Γₐ` represent n(x), `kernel` means the kernel function, `Gᵥ` is the
correlator. Note that `kernel` and `Gᵥ` have been rotated into singular
space. Please see comments in `init()` for more details.

### Arguments
See above explanations.

### Returns
* hτ -> α-resolved h(τ).
* Hα -> α-resolved Hc.
* Uα -> α-resolved internal energy, it is actually ⟨Hα⟩.

See also: [`calc_htau`](@ref).
"""
function calc_hamil(
    Γₚ::Array{I64,2},
    Γₐ::Array{F64,2},
    kernel::Matrix{F64},
    Gᵥ::Vector{F64}
    )
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
    calc_htau(
        Γₚ::Vector{I64},
        Γₐ::Vector{F64},
        kernel::Matrix{F64},
        Gᵥ:Vector{F64}
    )

Try to calculate α-dependent h(τ) via Eq.(36). `Γₚ` and `Γₐ` represent
n(x), `kernel` means the kernel function, `Gᵥ` is the correlator. Note
that `kernel` and `Gᵥ` have been rotated into singular space. Please
see comments in `init_context()` for more details.

### Arguments
See above explanations.

### Returns
* hτ -> α-resolved h(τ).

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

Generate a list for the α parameters.

### Arguments
N/A

### Returns
* αₗ -> List of the α parameters.
"""
function calc_alpha()
    nalph = get_a("nalph")
    alpha = get_a("alpha")
    ratio = get_a("ratio")

    αₗ = collect(alpha * (ratio ^ (x - 1)) for x in 1:nalph)

    return αₗ
end

"""
    constraints(S::StochACSolver, fmesh::AbstractMesh)

Try to implement the constrained stochastic analytic continuation
method. This function will return a collection. It contains all the
allowable indices. Be careful, the constrained stochastic analytic
continuation method is compatible with the self-adaptive mesh.

### Arguments
* S     -> A StochACSolver struct.
* fmesh -> Very dense mesh for the δ peaks.

### Returns
* allow -> Allowable indices.

See also: [`StochACSolver`](@ref).
"""
function constraints(S::StochACSolver, fmesh::AbstractMesh)
    exclude = get_b("exclude")
    nfine = get_a("nfine")
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
    try_move_s(
        i::I64,
        MC::StochACMC,
        SE::StochACElement,
        SC::StochACContext
    )

Select one δ function randomly and then change its position.

### Arguments
* i -> Index for α parameters.
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
N/A

See also: [`try_move_p`](@ref).
"""
function try_move_s(
    i::I64,
    MC::StochACMC,
    SE::StochACElement,
    SC::StochACContext
    )
    # Get current number of δ functions
    ngamm = get_a("ngamm")

    # Choose one δ function
    γ = rand(MC.rng, 1:ngamm)

    # Extract weight for the δ function
    a = SE.Γₐ[γ,i]

    # Choose new position for the δ function
    p = rand(MC.rng, SC.allow)

    # Try to calculate the change of Hc using Eq.~(42).
    hc = view(SC.hτ, :, i)
    Kₚ = view(SC.kernel, :, p)
    Kᵧ = view(SC.kernel, :, SE.Γₚ[γ,i])
    #
    δhc = a * (Kₚ - Kᵧ)
    δH = dot(δhc, 2.0 * hc + δhc)

    # Apply Metropolis algorithm
    MC.Mtry[i] = MC.Mtry[i] + 1
    if δH ≤ 0.0 || exp(-SC.αₗ[i] * δH) > rand(MC.rng)
        # Update Monte Carlo configurations
        SE.Γₚ[γ,i] = p

        # Update h(τ)
        @. hc = hc + δhc

        # Update Hc
        SC.Hα[i] = SC.Hα[i] + δH

        # Update Monte Carlo counter
        MC.Macc[i] = MC.Macc[i] + 1
    end
end

"""
    try_move_p(
        i::I64,
        MC::StochACMC,
        SE::StochACElement,
        SC::StochACContext
    )

Select two δ functions randomly and then change their positions.

### Arguments
* i -> Index for α parameters.
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
N/A

See also: [`try_move_s`](@ref).
"""
function try_move_p(
    i::I64,
    MC::StochACMC,
    SE::StochACElement,
    SC::StochACContext
    )
    # Get current number of δ functions
    ngamm = get_a("ngamm")
    #
    if ngamm < 2
        return
    end

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
    try_move_a(
        i::I64,
        MC::StochACMC,
        SE::StochACElement,
        SC::StochACContext
    )

Select two δ functions randomly and then change their weights.

### Arguments
* i -> Index for α parameters.
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
N/A

See also: [`try_move_x`](@ref).
"""
function try_move_a(
    i::I64,
    MC::StochACMC,
    SE::StochACElement,
    SC::StochACContext
    )
    # Get current number of δ functions
    ngamm = get_a("ngamm")
    #
    if ngamm < 2
        return
    end

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
    try_move_x(
        MC::StochACMC,
        SE::StochACElement,
        SC::StochACContext
    )

Try to exchange field configurations between two adjacent layers. Because
this function involves two layers, so it doesn't need the argument `i`.

### Arguments
* MC -> A StochACMC struct.
* SE -> A StochACElement struct.
* SC -> A StochACContext struct.

### Returns
N/A

See also: [`try_move_a`](@ref).
"""
function try_move_x(
    MC::StochACMC,
    SE::StochACElement,
    SC::StochACContext
    )
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
