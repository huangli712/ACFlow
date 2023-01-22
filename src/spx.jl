#
# Project : Gardenia
# Source  : spx.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2023/01/22
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
* Λ      -> Precomputed kernel matrix.
* Θ      -> Artificial inverse temperature.
* χ²min  -> Minimum of χ²min.
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
    Λ     :: Array{F64,2}
    Θ     :: F64
    χ²min :: F64
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
    # Initialize possible constraints. The array arrow contains all the
    # possible indices for poles.
    allow = constraints(S)

    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    SE = init_element(S, MC.rng, allow)
    println("Randomize Monte Carlo configurations")

    Gᵥ, σ¹ = init_iodata(S, rd)
    println("Postprocess input data: ", length(σ¹), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    fmesh = calc_fmesh(S)
    println("Build mesh for spectrum: ", length(mesh), " points")

    # Prepare the kernel matrix Λ. It is used to speed up the simulation.
    # Note that Λ depends on the type of kernel.
    ktype = get_b("ktype")
    χ₀ = -Gᵥ[1]
    #
    if     ktype == "fermi"
        Λ = calc_lambda(grid, fmesh)
    #
    elseif ktype == "boson"
        Λ = calc_lambda(grid, fmesh, χ₀, false)
    #
    elseif ktype == "bsymm"
        Λ = calc_lambda(grid, fmesh, χ₀, true)
    #
    end

    # Prepare some key variables
    Θ, χ²min, χ², Pᵥ, Aᵥ = init_context(S)

    # We have to make sure that the starting Gᵧ and χ² (i.e. χ²[1]) are
    # consistent with the current Monte Carlo configuration fields.
    Gᵧ = calc_green(SE.P, SE.A, Λ)
    χ²[1] = calc_chi2(Gᵧ, Gᵥ)

    SC = StochPXContext(Gᵥ, Gᵧ, σ¹, allow, grid, mesh, fmesh,
                        Λ, Θ, χ²min, χ², Pᵥ, Aᵥ)

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

        # Reset Gᵧ and χ²
        reset_context(t, SE, SC)

        # Apply simulated annealing algorithm
        for _ = 1:nstep
            sample(t, MC, SE, SC)
        end

        # Write Monte Carlo statistics
        write_statistics(MC)

        # Update χ²[t] to be consistent with SC.Pᵥ[t] and SC.Aᵥ[t]
        SC.χ²[t] = SC.χ²min
        @printf("try = %6i -> [χ² = %9.4e]\n", t, SC.χ²min)
        flush(stdout)
    end

    # Write pole expansion coefficients
    write_pole(SC.Pᵥ, SC.Aᵥ, SC.χ², SC.fmesh)

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

        # Reset Gᵧ and χ²
        reset_context(t, SE, SC)

        # Apply simulated annealing algorithm
        for _ = 1:nstep
            sample(t, MC, SE, SC)
        end

        # Write Monte Carlo statistics
        myid() == 2 && write_statistics(MC)

        # Update χ²[t] to be consistent with SC.Pᵥ[t] and SC.Aᵥ[t]
        SC.χ²[t] = SC.χ²min
        @printf("try = %6i -> [χ² = %9.4e]\n", t, SC.χ²min)
        flush(stdout)
    end

    # Write pole expansion coefficients
    myid() == 2 && write_pole(SC.Pᵥ, SC.Aᵥ, SC.χ², SC.fmesh)

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
    ktype = get_b("ktype")
    nmesh = get_b("nmesh")
    method = get_x("method")
    ntry = get_x("ntry")

    # Allocate memory
    # Gout: real frequency green's function, G(ω).
    # Gᵣ: imaginary frequency green's function, G(iωₙ)
    ngrid, _ = size(SC.Λ)
    Gout = zeros(C64, nmesh)
    Gᵣ = zeros(F64, ngrid)

    # Choose the best solution
    if method == "best"
        p = argmin(SC.χ²)
        χ₀ = -SC.Gᵥ[1]

        if     ktype == "fermi"
            Gout = calc_green(SC.Pᵥ[p], SC.Aᵥ[p], SC.mesh, SC.fmesh)
        #
        elseif ktype == "boson"
            Gout = calc_green(SC.Pᵥ[p], SC.Aᵥ[p], SC.mesh, SC.fmesh, χ₀, false)
        #
        elseif ktype == "bsymm"
            Gout = calc_green(SC.Pᵥ[p], SC.Aᵥ[p], SC.mesh, SC.fmesh, χ₀, true)
        #
        end

        Gᵣ = calc_green(SC.Pᵥ[p], SC.Aᵥ[p], SC.Λ)
        @printf("Best solution: try = %6i -> [χ² = %9.4e]\n", p, SC.χ²[p])
    #
    # Collect the `good` solutions and calculate their average.
    else
        # Calculate the median of SC.χ²
        chi2_med = median(SC.χ²)
        chi2_ave = mean(SC.χ²)

        # Determine the αgood parameter, which is used to filter the
        # calculated spectra.
        αgood = 1.2
        if count(x -> x < chi2_med / αgood, SC.χ²) ≤ ntry / 10
            αgood = 1.0
        end

        # Go through all the solutions
        c = 0.0
        χ₀ = -SC.Gᵥ[1]
        for i = 1:ntry
            if SC.χ²[i] < chi2_med / αgood
                if     ktype == "fermi"
                    G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.mesh, SC.fmesh)
                #
                elseif ktype == "boson"
                    G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.mesh, SC.fmesh, χ₀, false)
                #
                elseif ktype == "bsymm"
                    G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.mesh, SC.fmesh, χ₀, true)
                #
                end
                @. Gout = Gout + G
                #
                G = calc_green(SC.Pᵥ[i], SC.Aᵥ[i], SC.Λ)
                @. Gᵣ = Gᵣ + G
                #
                # Increase the counter
                c = c + 1.0
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
    end

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
end

#=
### *Core Algorithms*
=#

"""
    sample(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Try to search the configuration space to locate the minimum by using the
simulated annealing algorithm. Here, `t` means the t-th attempt.
"""
function sample(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
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
    Sacc = 0
    Stry = 0
    Pacc = 0
    Ptry = 0
    Aacc = 0
    Atry = 0
    Xacc = 0
    Xtry = 0

    MC = StochPXMC(rng, Sacc, Stry, Pacc, Ptry, Aacc, Atry, Xacc, Xtry)

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

    # We have to make sure ∑ᵢ Aᵢ = 1
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
    Θ = get_x("theta")

    χ²min = 1e10
    χ² = zeros(F64, ntry)

    Pᵥ = Vector{I64}[]
    Aᵥ = Vector{F64}[]
    for _ = 1:ntry
        push!(Pᵥ,  ones(I64, npole))
        push!(Aᵥ, zeros(F64, npole))
    end

    return Θ, χ²min, χ², Pᵥ, Aᵥ
end

"""
    reset_mc(MC::StochPXMC)

Reset the counters in StochPXMC struct.
"""
function reset_mc(MC::StochPXMC)
    MC.Sacc = 0
    MC.Stry = 0
    MC.Pacc = 0
    MC.Ptry = 0
    MC.Aacc = 0
    MC.Atry = 0
    MC.Xacc = 0
    MC.Xtry = 0
end

"""
    reset_element(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)

Reset the Monte Carlo field configurations (i.e. positions and amplitudes
of the poles).
"""
function reset_element(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
    npole = get_x("npole")
    if npole ≤ 5
        if 4 ≤ npole ≤ 5
            nselect = 2
        else
            nselect = 1
        end
    else
        nselect = ceil(I64, npole / 5)
    end
    @assert nselect ≤ npole

    selected = rand(rng, 1:npole, nselect)
    unique!(selected)
    nselect = length(selected)

    if rand(rng) < 0.5
        P = rand(rng, allow, nselect)
        @. SE.P[selected] = P
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

"""
    reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)

Recalculate imaginary frequency green's function and goodness-of-fit
function by new Monte Carlo field configurations for the t-th attempts.
"""
function reset_context(t::I64, SE::StochPXElement, SC::StochPXContext)
    Gᵧ = calc_green(SE.P, SE.A, SC.Λ)
    χ² = calc_chi2(Gᵧ, SC.Gᵥ)

    @. SC.Gᵧ = Gᵧ
    SC.χ²[t] = χ²
    SC.χ²min = 1e10
    SC.Θ = get_x("theta")
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
    calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh)

Precompute the kernel matrix Λ (Λ ≡ 1 / (iωₙ - ϵ)).
It is for the fermionic systems.
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
    calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)

Precompute the kernel matrix Λ. Here, `χ₀` is actually -G(iωₙ = 0). And
the argument `bsymm` is used to distinguish two different bosonic kernels.
If `bsymm` is false, it means that the kernel is `boson`. If `bsymm` is
true, the kernel is `bsymm`. This function is for the bosonic systems.
"""
function calc_lambda(grid::AbstractGrid, fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)
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
    calc_green(P::Vector{I64},
               A::Vector{F64},
               Λ::Array{F64,2})

Reconstruct green's function at imaginary axis by the pole expansion.
"""
function calc_green(P::Vector{I64}, A::Vector{F64}, Λ::Array{F64,2})
    # Note that here `ngrid` is equal to 2 × ngrid sometimes.
    ngrid, _ = size(Λ)

    G = zeros(F64, ngrid)
    for i = 1:ngrid
        G[i] = dot(A, Λ[i,P])
    end

    return G
end

"""
    calc_green(P::Vector{I64},
               A::Vector{F64},
               mesh::AbstractMesh,
               fmesh::AbstractMesh)

Reconstruct green's function at real axis by the pole expansion. It is
for the fermionic systems only.
"""
function calc_green(P::Vector{I64},
                    A::Vector{F64},
                    mesh::AbstractMesh,
                    fmesh::AbstractMesh)
    η = get_x("eta")
    nmesh = length(mesh)

    iωₙ = mesh.mesh .+ im * η
    G = zeros(C64, nmesh)
    for i in eachindex(mesh)
        G[i] = sum( @. A / (iωₙ[i] - fmesh.mesh[P]) )
    end

    return G
end

"""
    calc_green(P::Vector{I64},
               A::Vector{F64},
               mesh::AbstractMesh,
               fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)

Reconstruct green's function at real axis by the pole expansion. Here,
`χ₀` is actually -G(iωₙ = 0). And the argument `bsymm` is used to
distinguish two different bosonic kernels. If `bsymm` is false, it means
that the kernel is `boson`. If `bsymm` is true, the kernel is `bsymm`.
It is for the bosonic systems only.
"""
function calc_green(P::Vector{I64},
                    A::Vector{F64},
                    mesh::AbstractMesh,
                    fmesh::AbstractMesh, χ₀::F64, bsymm::Bool)
    η = get_x("eta")
    nmesh = length(mesh)

    iωₙ = mesh.mesh .+ im * η
    G = zeros(C64, nmesh)
    if bsymm == false
        _A = A .* χ₀ .* fmesh.mesh[P]
        for i in eachindex(mesh)
            G[i] = sum( @. _A / (iωₙ[i] - fmesh.mesh[P]) )
        end
    #
    else
        _A = A .* χ₀ .* fmesh.mesh[P] .* 0.5
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
    try_move_s(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Change the position of one randomly selected pole.

See also: [`try_move_p`](@ref).
"""
function try_move_s(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Get parameters
    ngrid = length(SC.Gᵧ) # get_b("ngrid")
    nfine = get_x("nfine")
    npole = get_x("npole")
    move_window = ceil(I64, nfine / 100)

    # It is used to save the change of green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

        # Select one pole randomly
        s = rand(MC.rng, 1:npole)

        # Try to change position of the s pole
        Aₛ = SE.A[s]
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
        !(P₂ in SC.allow) && continue

        # Calculate change of green's function
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        @. δG = Aₛ * (Λ₂ - Λ₁)

        # Calculate new green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²[t]

        # Simulated annealing algorithm
        MC.Stry = MC.Stry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.P[s] = P₂

            # Update reconstructed green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ²[t] = χ²

            # Update Monte Carlo counter
            MC.Sacc = MC.Sacc + 1

            # Save optimal solution
            if χ² < SC.χ²min
                SC.χ²min = χ²
                measure(t, SE, SC)
            end
        end

    end
end

"""
    try_move_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Change the positions of two randomly selected poles.

See also: [`try_move_s`](@ref).
"""
function try_move_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Get parameters
    ngrid = length(SC.Gᵧ) # get_b("ngrid")
    npole = get_x("npole")
    #
    if npole == 1
        return
    end

    # It is used to save the change of green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

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
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        Λ₃ = view(SC.Λ, :, P₃)
        Λ₄ = view(SC.Λ, :, P₄)
        @. δG = A₁ * (Λ₃ - Λ₁) + A₂ * (Λ₄ - Λ₂)

        # Calculate new green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²[t]

        # Simulated annealing algorithm
        MC.Ptry = MC.Ptry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.P[s₁] = P₃
            SE.P[s₂] = P₄

            # Update reconstructed green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ²[t] = χ²

            # Update Monte Carlo counter
            MC.Pacc = MC.Pacc + 1

            # Save optimal solution
            if χ² < SC.χ²min
                SC.χ²min = χ²
                measure(t, SE, SC)
            end
        end

    end
end

"""
    try_move_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Change the amplitudes of two randomly selected poles.

See also: [`try_move_x`](@ref).
"""
function try_move_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Get parameters
    ngrid = length(SC.Gᵧ) # get_b("ngrid")
    npole = get_x("npole")
    #
    if npole == 1
        return
    end

    # It is used to save the change of green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

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
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        @. δG = (A₃ - A₁) * Λ₁ + (A₄ - A₂) * Λ₂

        # Calculate new green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²[t]

        # Simulated annealing algorithm
        MC.Atry = MC.Atry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.A[s₁] = A₃
            SE.A[s₂] = A₄

            # Update reconstructed green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ²[t] = χ²

            # Update Monte Carlo counter
            MC.Aacc = MC.Aacc + 1

            # Save optimal solution
            if χ² < SC.χ²min
                SC.χ²min = χ²
                measure(t, SE, SC)
            end
        end

    end
end

"""
    try_move_x(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Exchange the amplitudes of two randomly selected poles.

See also: [`try_move_a`](@ref).
"""
function try_move_x(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # Get parameters
    ngrid = length(SC.Gᵧ) # get_b("ngrid")
    npole = get_x("npole")
    #
    if npole == 1
        return
    end

    # It is used to save the change of green's function
    δG = zeros(F64, ngrid)
    Gₙ = zeros(F64, ngrid)

    # Try to go through each pole
    for _ = 1:npole

        # Select two poles randomly
        s₁ = 1
        s₂ = 1
        while s₁ == s₂
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

        # Calculate change of green's function
        Λ₁ = view(SC.Λ, :, P₁)
        Λ₂ = view(SC.Λ, :, P₂)
        @. δG = (A₃ - A₁) * Λ₁ + (A₄ - A₂) * Λ₂

        # Calculate new green's function and goodness-of-fit function
        @. Gₙ = δG + SC.Gᵧ
        χ² = calc_chi2(Gₙ, SC.Gᵥ)
        δχ² = χ² - SC.χ²[t]

        # Simulated annealing algorithm
        MC.Xtry = MC.Xtry + 1
        if δχ² < 0 || min(1.0, exp(-δχ² * SC.Θ)) > rand(MC.rng)
            # Update Monte Carlo configuration
            SE.A[s₁] = A₃
            SE.A[s₂] = A₄

            # Update reconstructed green's function
            @. SC.Gᵧ = Gₙ

            # Update goodness-of-fit function
            SC.χ²[t] = χ²

            # Update Monte Carlo counter
            MC.Xacc = MC.Xacc + 1

            # Save optimal solution
            if χ² < SC.χ²min
                SC.χ²min = χ²
                measure(t, SE, SC)
            end
        end

    end
end
