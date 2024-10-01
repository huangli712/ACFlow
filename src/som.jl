#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/01
#

#=
### *Customized Structs* : *StochOM Solver*
=#

"""
    Box

Rectangle. The field configuration consists of many boxes. They exhibit
various areas (width × height). We use the Metropolis important sampling
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
be sampled by Monte Carlo sweeping procedure.

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
* 𝕊ᵥ    -> It is used to interpolate the Λ functions.
"""
mutable struct StochOMContext
    Gᵥ   :: Vector{F64}
    σ¹   :: Vector{F64}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    Cᵥ   :: Vector{Vector{Box}}
    Δᵥ   :: Vector{F64}
    𝕊ᵥ   :: Vector{CubicSplineInterpolation}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochOMSolver, rd::RawData)

Solve the analytic continuation problem by the stochastic optimization
method. This solver requires a lot of computational resources to get
reasonable results. It is suitable for both Matsubara and imaginary
time correlators. It is the driver for the StochOM solver.

If the input correlators are bosonic, this solver will return A(ω) / ω
via `Aout`, instead of A(ω). At this time, `Aout` is not compatible with
`Gout`. If the input correlators are fermionic, this solver will return
A(ω) in `Aout`. Now it is compatible with `Gout`. These behaviors are just
similar to the MaxEnt, StochAC, and StochSK solvers.

Now the StochOM solver supports both continuous and δ-like spectra.

### Arguments
* S -> A StochOMSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* mesh -> Real frequency mesh, ω.
* Aout -> Spectral function, A(ω).
* Gout -> Retarded Green's function, G(ω).
"""
function solve(S::StochOMSolver, rd::RawData)
    println("[ StochOM ]")
    MC, SC = init(S, rd)

    # Parallel version
    if nworkers() > 1
        #
        println("Using $(nworkers()) workers")
        #
        # Copy configuration dicts
        p1 = deepcopy(PBASE)
        p2 = deepcopy(PStochOM)
        #
        #  Launch the tasks one by one
        𝐹 = Future[]
        for i = 1:nworkers()
            𝑓 = @spawnat i + 1 prun(S, p1, p2, MC, SC)
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
        Aout = similar(sol[end])
        fill!(Aout, 0.0)
        for i in eachindex(sol)
            @. Aout = Aout + sol[i] / nworkers()
        end
        #
        # Postprocess the solutions
        Gout = last(SC, Aout)
        #
    # Sequential version
    else
        #
        Aout = run(MC, SC)
        Gout = last(SC, Aout)
        #
    end

    return SC.mesh.mesh, Aout, Gout
end

"""
    init(S::StochOMSolver, rd::RawData)

Initialize the StochOM solver and return the StochOMMC and StochOMContext
structs. Please don't call this function directly.

### Arguments
* S -> A StochOMSolver struct.
* rd -> A RawData struct, containing raw data for input correlator.

### Returns
* MC -> A StochOMMC struct.
* SC -> A StochOMContext struct.
"""
function init(S::StochOMSolver, rd::RawData)
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

    # Prepare some key variables
    Cᵥ, Δᵥ, 𝕊ᵥ = init_context(S, grid)
    SC = StochOMContext(Gᵥ, σ¹, grid, mesh, Cᵥ, Δᵥ, 𝕊ᵥ)
    println("Initialize context for the StochOM solver")

    return MC, SC
end

"""
    run(MC::StochOMMC, SC::StochOMContext)

Perform stochastic optimization simulation, sequential version.

### Arguments
* MC -> A StochOMMC struct.
* SC -> A StochOMContext struct.

### Returns
* Aout -> Spectral function, A(ω).
"""
function run(MC::StochOMMC, SC::StochOMContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Setup essential parameters
    ntry = get_s("ntry")
    nstep = get_s("nstep")

    # Sample and collect data
    println("Start stochastic sampling...")
    for l = 1:ntry
        # Re-initialize the simulation
        SE = init_element(MC, SC)

        # For each attempt, we should perform `nstep × N` Monte Carlo
        # updates, where `N` means length of the Markov chain.
        for _ = 1:nstep
            update(MC, SE, SC)
        end

        # Write Monte Carlo statistics
        l % 10 == 0 && fwrite && write_statistics(MC)

        # Accumulate the data
        SC.Δᵥ[l] = SE.Δ
        SC.Cᵥ[l] = deepcopy(SE.C)

        # Show error function for the current attempt
        @printf("try -> %6i (%6i) Δ -> %8.4e \n", l, ntry, SE.Δ)
        flush(stdout)
    end

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    prun(
        S::StochOMSolver,
        p1::Dict{String,Vector{Any}},
        p2::Dict{String,Vector{Any}},
        MC::StochOMMC,
        SC::StochOMContext
    )

Perform stochastic optimization simulation, parallel version.
The arguments `p1` and `p2` are copies of PBASE and PStochOM, respectively.

### Arguments
* S -> A StochOMSolver struct.
* p1 -> A copy of PBASE.
* p2 -> A copy of PStochOM.
* MC -> A StochOMMC struct.
* SC -> A StochOMContext struct.

### Returns
* Aout -> Spectral function, A(ω).
"""
function prun(
    S::StochOMSolver,
    p1::Dict{String,Vector{Any}},
    p2::Dict{String,Vector{Any}},
    MC::StochOMMC,
    SC::StochOMContext
    )
    # Revise parameteric dicts
    rev_dict_b(p1)
    rev_dict_s(S, p2)

    # Initialize random number generator again
    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Setup essential parameters
    ntry = get_s("ntry")
    nstep = get_s("nstep")

    # Sample and collect data
    println("Start stochastic sampling...")
    for l = 1:ntry
        # Re-initialize the simulation
        SE = init_element(MC, SC)

        # For each attempt, we should perform `nstep × N` Monte Carlo
        # updates, where `N` means length of the Markov chain.
        for _ = 1:nstep
            update(MC, SE, SC)
        end

        # Write Monte Carlo statistics
        myid() == 2 && l % 10 == 0 && fwrite && write_statistics(MC)

        # Accumulate the data
        SC.Δᵥ[l] = SE.Δ
        SC.Cᵥ[l] = deepcopy(SE.C)

        # Show error function for the current attempt
        @printf("try -> %6i (%6i) Δ -> %8.4e \n", l, ntry, SE.Δ)
        flush(stdout)
    end

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    average(SC::StochOMContext)

Postprocess the collected results after the stochastic optimization
simulations. It will generate the spectral functions.

### Arguments
* SC -> A StochOMContext struct.

### Returns
* Aom -> Spectral function, A(ω).
"""
function average(SC::StochOMContext)
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    nmesh = get_b("nmesh")
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
    passed = I64[]
    for l = 1:ntry
        # Filter the reasonable spectra
        if SC.Δᵥ[l] < dev_ave / αgood
            # Generate the spectrum, and add it to Aom.
            for w = 1:nmesh
                _omega = SC.mesh[w]
                # Scan all boxes
                for r = 1:length(SC.Cᵥ[l])
                    R = SC.Cᵥ[l][r]
                    # Yes, this box contributes. The point, _omega, is
                    # covered by this box.
                    if R.c - 0.5 * R.w ≤ _omega ≤ R.c + 0.5 * R.w
                        Aom[w] = Aom[w] + R.h
                    end
                end
            end
            #
            # Record which spectrum is used
            append!(passed, l)
        end
    end

    # Normalize the spectrum
    Lgood = count(x -> x < dev_ave / αgood, SC.Δᵥ)
    @assert Lgood == length(passed)
    @. Aom = Aom / Lgood
    @printf("Median χ² : %16.12e\n", dev_ave)
    @printf("Accepted configurations : %5i\n", Lgood)

    # Write indices of selected solutions
    if nworkers() > 1
        myid() == 2 && fwrite && write_passed(passed, dev_ave, αgood)
    else
        fwrite && write_passed(passed, dev_ave, αgood)
    end

    return Aom
end

"""
    last(SC::StochOMContext, Aout::Vector{F64})

It will process and write the calculated results by the StochOM solver,
including final spectral function and reproduced correlator.

### Arguments
* SC   -> A StochOMContext struct.
* Aout -> Spectral function, A(ω).

### Returns
* G -> Retarded Green's function, G(ω).
"""
function last(SC::StochOMContext, Aout::Vector{F64})
    # By default, we should write the analytic continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Write the spectral function
    fwrite && write_spectrum(SC.mesh, Aout)

    # Reproduce input data and write them
    kernel = make_kernel(SC.mesh, SC.grid)
    G = reprod(SC.mesh, kernel, Aout)
    fwrite && write_backward(SC.grid, G)

    # Calculate full response function on real axis and write them
    if get_b("ktype") == "fermi"
        _G = kramers(SC.mesh, Aout)
    else
        _G = kramers(SC.mesh, Aout .* SC.mesh)
    end
    fwrite && write_complete(SC.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

"""
    update(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext)

Using the Metropolis algorithm to update the field configuration, i.e, a
collection of hundreds of boxes. Be careful, this function only updates
the Monte Carlo configurations (in other words, `SE`). It doesn't record
them. Measurements are done in `run()` and `prun()`. This is the reason
why this function is named as `update()`, instead of `sample()`.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.

### Returns
N/A
"""
function update(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext)
    Tmax = 100 # Length of the Markov chain
    nbox = get_s("nbox")

    ST = deepcopy(SE)

    # The Markov chain is divided into two stages
    T1 = rand(MC.rng, 1:Tmax)
    d1 = rand(MC.rng, F64)
    d2 = 1.0 + rand(MC.rng, F64)
    #
    # The first stage
    for _ = 1:T1
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if 1 ≤ length(ST.C) ≤ nbox - 1
                    try_insert(MC, ST, SC, d1)
                end
                break

            @case 2
                if length(ST.C) ≥ 2
                    try_remove(MC, ST, SC, d1)
                end
                break

            @case 3
                if length(ST.C) ≥ 1
                    try_shift(MC, ST, SC, d1)
                end
                break

            @case 4
                if length(ST.C) ≥ 1
                    try_width(MC, ST, SC, d1)
                end
                break

            @case 5
                if length(ST.C) ≥ 2
                    try_height(MC, ST, SC, d1)
                end
                break

            @case 6
                if 1 ≤ length(ST.C) ≤ nbox - 1
                    try_split(MC, ST, SC, d1)
                end
                break

            @case 7
                if length(ST.C) ≥ 2
                    try_merge(MC, ST, SC, d1)
                end
                break
        end

    end
    #
    # The second stage
    for _ = T1+1:Tmax
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if 1 ≤ length(ST.C) ≤ nbox - 1
                    try_insert(MC, ST, SC, d2)
                end
                break

            @case 2
                if length(ST.C) ≥ 2
                    try_remove(MC, ST, SC, d2)
                end
                break

            @case 3
                if length(ST.C) ≥ 1
                    try_shift(MC, ST, SC, d2)
                end
                break

            @case 4
                if length(ST.C) ≥ 1
                    try_width(MC, ST, SC, d2)
                end
                break

            @case 5
                if length(ST.C) ≥ 2
                    try_height(MC, ST, SC, d2)
                end
                break

            @case 6
                if 1 ≤ length(ST.C) ≤ nbox - 1
                    try_split(MC, ST, SC, d2)
                end
                break

            @case 7
                if length(ST.C) ≥ 2
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
    init_iodata(S::StochOMSolver, rd::RawData)

Preprocess the input data (`rd`).

### Arguments
* S -> A StochOMSolver struct.
* rd -> A RawData struct, which contains essential input data.

### Returns
* Gᵥ -> Input correlator.
* σ¹ -> 1.0 / σ¹.

See also: [`RawData`](@ref), [`GreenData`](@ref).
"""
function init_iodata(S::StochOMSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value
    σ¹ = 1.0 ./ G.error

    return Gᵥ, σ¹
end

"""
    init_mc(S::StochOMSolver)

Try to create a StochOMMC struct. Some counters for Monte Carlo updates
are initialized here.

### Arguments
* S -> A StochOMSolver struct.

### Returns
* MC -> A StochOMMC struct.

See also: [`StochOMMC`](@ref).
"""
function init_mc(S::StochOMSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    #
    Macc = zeros(I64, 7)
    Mtry = zeros(I64, 7)
    #
    MC = StochOMMC(rng, Macc, Mtry)

    return MC
end

"""
    init_element(MC::StochOMMC, SC::StochOMContext)

Try to initialize a StochOMElement struct. In other words, we should
randomize the configurations for future Monte Carlo sampling here.

### Arguments
* MC -> A StochOMMC struct.
* SC -> A StochOMContext struct.

### Returns
* SE -> A StochOMElement struct.

See also: [`StochOMElement`](@ref).
"""
function init_element(MC::StochOMMC, SC::StochOMContext)
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    nbox = get_s("nbox")
    sbox = get_s("sbox")
    wbox = get_s("wbox")

    # Generate weights randomly
    _Know = rand(MC.rng, 2:nbox)
    _weight = zeros(F64, _Know)
    for i = 1:_Know
        _weight[i] = rand(MC.rng, F64)
    end
    _weight[end] = 1.0

    # Sort weights, make sure the sum of weights is always 1.0.
    sort!(_weight)
    weight = diff(_weight)
    insert!(weight, 1, _weight[1])
    sort!(weight)

    # Make sure that each weight is larger than sbox.
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

    # Create some boxes with random c, w, and h.
    C = Box[]
    Λ = zeros(F64, length(SC.Gᵥ), nbox)
    Δ = 0.0
    #
    for k = 1:_Know
        c = wmin + wbox / 2.0 + (wmax - wmin - wbox) * rand(MC.rng, F64)
        w = wbox + (min(2.0 * (c - wmin), 2.0 * (wmax - c)) - wbox) * rand(MC.rng, F64)
        while !constraints(c - w/2, c + w/2)
            c = wmin + wbox / 2.0 + (wmax - wmin - wbox) * rand(MC.rng, F64)
            w = wbox + (min(2.0 * (c - wmin), 2.0 * (wmax - c)) - wbox) * rand(MC.rng, F64)
        end
        h = weight[k] / w
        R = Box(h, w, c)
        push!(C, R)
        Λ[:,k] .= eval_lambda(R, SC.grid, SC.𝕊ᵥ)
    end
    #
    # Calculate Green's function and relative error using boxes
    G = calc_green(Λ, _Know)
    Δ = calc_error(G, SC.Gᵥ, SC.σ¹)

    return StochOMElement(C, Λ, G, Δ)
end

"""
    init_context(S::StochOMSolver, grid::AbstractGrid)

Try to initialize the key members of a StochOMContext struct.

### Arguments
* S -> A StochOMSolver struct.
* grid -> Grid for input data.

### Returns
* Cᵥ -> Field configurations for all attempts.
* Δᵥ -> Errors for all attempts.
* 𝕊ᵥ -> Interpolators for the Λ functions.

See also: [`StochOMContext`](@ref).
"""
function init_context(S::StochOMSolver, grid::AbstractGrid)
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    ntry = get_s("ntry")
    nbox = get_s("nbox")

    # If we increase nmesh gradually, perhaps we could get more precise
    # interpolants 𝕊ᵥ.
    nmesh = get_b("nmesh") # nmesh = 101, 201, 301, ...
    ngrid = get_b("ngrid")
    @assert ngrid == length(grid)

    # Initialize errors
    Δᵥ = zeros(F64, ntry)

    # Initialize field configurations (boxes)
    Cᵥ = []
    for _ = 1:ntry
        C = Box[]
        for _ = 1:nbox
            push!(C, Box(0.0, 0.0, 0.0))
        end
        push!(Cᵥ, C)
    end

    # Initialize interpolants 𝕊ᵥ
    # It is useful only when the input data is in imaginary time axis.
    𝕊ᵥ = Vector{CubicSplineInterpolation}(undef, ngrid)
    #
    if get_b("grid") in ("ftime", "fpart", "btime", "bpart")
        # Create linear mesh for the interpolants
        am = LinearMesh(nmesh, wmin, wmax)

        # Calculate the interpolants at the nodes
        #
        # Initialize memory
        # See below remarks for the Λ function
        Λ = zeros(F64, ngrid, nmesh)
        #
        # Just evaluate the interpolants
        for m in eachindex(am)
            if m > 1
                # Create linear mesh for the integrand
                cm = LinearMesh(nmesh, wmin, am[m])
                #
                # Calculate the integrand, i.e., the kernel.
                K = make_kernel(cm, grid)
                #
                # Calculate the integral using trapz rule. Perhaps more
                # precise algorithms should be used.
                for i = 1:ngrid
                    Λ[i,m] = trapz(cm, K[i,:])
                end
            end
        end

        # Create CubicSplineInterpolation structs in time grid τ
        for i = 1:ngrid
            𝕊ᵥ[i] = CubicSplineInterpolation(Λ[i,:], am.mesh)
        end
    end

    return Cᵥ, Δᵥ, 𝕊ᵥ
end

#=
*Remarks* :

Here, we would like to calculate the ``\Lambda`` function:

```math
\begin{equation}
\Lambda(\omega_n) = \int^{\infty}_{-\infty}
    d\omega~K(\omega_n,\omega) A(\omega).
\end{equation}
```

For a given box ``R``, its contribution to the ``\Lambda`` function is

```math
\begin{equation}
Λ_{R}(\omega_n) = h \int^{c+w/2}_{c-w/2}
    d\omega~K(\omega_n,\omega),
\end{equation}
```

where ``h``, ``w``, and ``c`` denote the height, width, and center of the
box ``R``. Next, we will show you how to calculate ``\Lambda_R(\omega_n)``.
We can use `sympy` to derive the integral formula.

---

**A** For fermionic Matsubara frequency system.

```math
\begin{equation}
Λ_{R}(\omega_n) = h \int^{c+w/2}_{c-w/2}
    d\omega~\frac{1}{i\omega_n - \omega}.
\end{equation}
```

The python codes are as follows:

```python
>>> from sympy import *
>>> w, x, a, b, h = symbols("w x a b h")
>>> integrate(h/(w - x), (x, a, b))
h*log(a - w) - h*log(b - w)
```

So, we have

```math
\begin{equation}
Λ_{R}(\omega_n) = h \log
    \left(\frac{i\omega_n - c + w/2}{i\omega_n - c - w/2}\right).
\end{equation}
```

---

**B** For bosonic Matsubara frequency system.

```math
\begin{equation}
Λ_{R}(\omega_n) = h \int^{c+w/2}_{c-w/2}
    d\omega~\frac{\omega}{i\omega_n - \omega}.
\end{equation}
```

The python codes are as follows:

```python
>>> from sympy import *
>>> w, x, a, b, h = symbols("w x a b h")
>>> f = integrate(h*x/(w - x), (x, a, b))
>>> simplify(f)
h*(a - b + w*log(a - w) - w*log(b - w))
```

So, we have

```math
\begin{equation}
Λ_{R}(\omega_n) = h * \left[-w + i\omega_n * \log
    \left(\frac{i\omega_n - c + w/2}{i\omega_n - c - w/2}\right)
\right].
\end{equation}
```

---

**C** For bosonic Matsubara frequency system (symmetric version).

```math
\begin{equation}
Λ_{R}(\omega_n) = h \int^{c+w/2}_{c-w/2}
    d\omega~\frac{-2\omega^2}{\omega_n^2 + \omega^2}.
\end{equation}
```

We can utilize the following formula:

```math
\begin{equation}
\int \frac{x^2}{a^2 + x^2} dx =
    x - a~\tan^{-1}\left(\frac{x}{a}\right)
\end{equation}
```

Thus,

```math
\begin{equation}
Λ_{R}(\omega_n) = -2h
    \left[
        \omega - \omega_n \tan^{-1} \left(\frac{\omega}{\omega_n}\right)
    \right]
    \bigg|^{c+w/2}_{c-w/2}
\end{equation}
```

Finally, we have

```math
\begin{equation}
Λ_{R}(\omega_n) = -2h
    \left[
        w +
        \omega_n \tan^{-1} \left(\frac{c - w/2}{\omega_n}\right)
        -
        \omega_n \tan^{-1} \left(\frac{c + w/2}{\omega_n}\right)
    \right]
\end{equation}
```

---

**D** For imaginary time axis.

If the kernel function ``K`` is defined in imaginary time axis, then the
``\Lambda`` and ``\Lambda_R`` functions become:

```math
\begin{equation}
\Lambda(\tau) = \int^{\infty}_{-\infty}
    d\omega~K(\tau,\omega) A(\omega),
\end{equation}
```

```math
\begin{equation}
Λ_{R}(\tau) = h \int^{c+w/2}_{c-w/2}
    d\omega~K(\tau,\omega).
\end{equation}
```

We just evaluate the following integral numerically

```math
\begin{equation}
𝕊_R(\tau,Ω) = \int^{Ω}_{\omega_{\text{min}}}
    d\omega~K(\tau,\omega).
\end{equation}
```

We notice that `` \Omega \in (\omega_{\text{min}},\omega_{\text{max}}]``,
and ``\Omega`` is defined in a dense linear mesh. Then we reach

```math
\begin{equation}
Λ_{R}(\tau) = h \left[ 𝕊_R(\tau,c + w/2) - 𝕊_R(\tau,c - w/2) \right].
\end{equation}
```

In the present implementation, Eq.(13) is evaluated by trapz algorithm,
and Eq.(14) is evaluate using cubic spline interpolation.

We have implemented the above formulas in `eval_lambda()`.
=#

"""
    eval_lambda(
        r::Box,
        grid::FermionicMatsubaraGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for FermionicMatsubaraGrid only. Because there is an
analytic expression for this case, 𝕊 is useless.

Actually, 𝕊 is undefined here. See init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(iωₙ) function, 1D function.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::FermionicMatsubaraGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Evaluate Λ
    iw = im * grid.ω
    Λ = @. r.h * log((iw - e₁) / (iw - e₂))

    return vcat(real(Λ), imag(Λ))
end

"""
    eval_lambda(
        r::Box,
        grid::FermionicFragmentMatsubaraGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for FermionicFragmentMatsubaraGrid only. Because there
is an analytic expression for this case, 𝕊 is useless.

Actually, 𝕊 is undefined here. See init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(iωₙ) function, 1D function.

See also: [`FermionicFragmentMatsubaraGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::FermionicFragmentMatsubaraGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Evaluate Λ
    iw = im * grid.ω
    Λ = @. r.h * log((iw - e₁) / (iw - e₂))

    return vcat(real(Λ), imag(Λ))
end

"""
    eval_lambda(
        r::Box,
        grid::FermionicImaginaryTimeGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for FermionicImaginaryTimeGrid only. Since there is
no analytic expressions for this case, the cubic spline interpolation
algorithm is adopted. Here, 𝕊 is initialized in init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(τ) function, 1D function.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::FermionicImaginaryTimeGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Initialize Λ function
    ntime = grid.ntime
    Λ = zeros(F64, ntime)

    # 𝕊ᵢ(e₂): integral boundary is from wmin to e₂
    # 𝕊ᵢ(e₁): integral boundary is from wmin to e₁
    for i = 1:ntime
        Λ[i] = ( 𝕊[i](e₂) - 𝕊[i](e₁) ) * r.h
    end

    return Λ
end

"""
    eval_lambda(
        r::Box,
        grid::FermionicFragmentTimeGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for FermionicFragmentTimeGrid only. Since there is
no analytic expressions for this case, the cubic spline interpolation
algorithm is adopted. Here, 𝕊 is initialized in init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(τ) function, 1D function.

See also: [`FermionicFragmentTimeGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::FermionicFragmentTimeGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Initialize Λ function
    ntime = grid.ntime
    Λ = zeros(F64, ntime)

    # 𝕊ᵢ(e₂): integral boundary is from wmin to e₂
    # 𝕊ᵢ(e₁): integral boundary is from wmin to e₁
    for i = 1:ntime
        Λ[i] = ( 𝕊[i](e₂) - 𝕊[i](e₁) ) * r.h
    end

    return Λ
end

"""
    eval_lambda(
        r::Box,
        grid::BosonicMatsubaraGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for BosonicMatsubaraGrid only. Because there is an
analytic expression for this case, 𝕊 is useless.

Actually, 𝕊 is undefined here. See init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(iωₙ) function, 1D function.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::BosonicMatsubaraGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get type of bosonic kernel
    ktype = get_b("ktype")

    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Evaluate Λ
    if ktype == "bsymm"
        Λ = @. atan( e₁ / grid.ω ) - atan( e₂ / grid.ω )
        Λ = -2.0 * r.h * (r.w .+ grid.ω .* Λ)
        return Λ
    else
        iw = im * grid.ω
        Λ = @. r.h * (-r.w + iw * log((iw - e₁) / (iw - e₂)))
        return vcat(real(Λ), imag(Λ))
    end
end

"""
    eval_lambda(
        r::Box,
        grid::BosonicFragmentMatsubaraGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for BosonicFragmentMatsubaraGrid only. Because there
is an analytic expression for this case, 𝕊 is useless.

Actually, 𝕊 is undefined here. See init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(iωₙ) function, 1D function.

See also: [`BosonicFragmentMatsubaraGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::BosonicFragmentMatsubaraGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get type of bosonic kernel
    ktype = get_b("ktype")

    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Evaluate Λ
    if ktype == "bsymm"
        Λ = @. atan( e₁ / grid.ω ) - atan( e₂ / grid.ω )
        Λ = -2.0 * r.h * (r.w .+ grid.ω .* Λ)
        return Λ
    else
        iw = im * grid.ω
        Λ = @. r.h * (-r.w + iw * log((iw - e₁) / (iw - e₂)))
        return vcat(real(Λ), imag(Λ))
    end
end

"""
    eval_lambda(
        r::Box,
        grid::BosonicImaginaryTimeGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for BosonicImaginaryTimeGrid only. Since there is
no analytic expressions for this case, the cubic spline interpolation
algorithm is adopted. Here, 𝕊 is initialized in init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(τ) function, 1D function.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::BosonicImaginaryTimeGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Initialize Λ function
    ntime = grid.ntime
    Λ = zeros(F64, ntime)

    # 𝕊ᵢ(e₂): integral boundary is from wmin to e₂
    # 𝕊ᵢ(e₁): integral boundary is from wmin to e₁
    for i = 1:ntime
        Λ[i] = ( 𝕊[i](e₂) - 𝕊[i](e₁) ) * r.h
    end

    return Λ
end

"""
    eval_lambda(
        r::Box,
        grid::BosonicFragmentTimeGrid,
        𝕊::Vector{<:AbstractInterpolation}
    )

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for BosonicFragmentTimeGrid only. Since there is
no analytic expressions for this case, the cubic spline interpolation
algorithm is adopted. Here, 𝕊 is initialized in init_context().

### Arguments
* r    -> A box or rectangle.
* grid -> Imaginary axis grid for input data.
* 𝕊    -> An interpolant.

### Returns
* Λ -> Λ(τ) function, 1D function.

See also: [`BosonicFragmentTimeGrid`](@ref).
"""
function eval_lambda(
    r::Box,
    grid::BosonicFragmentTimeGrid,
    𝕊::Vector{<:AbstractInterpolation}
    )
    # Get left and right boundaries of the given box
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    # Initialize Λ function
    ntime = grid.ntime
    Λ = zeros(F64, ntime)

    # 𝕊ᵢ(e₂): integral boundary is from wmin to e₂
    # 𝕊ᵢ(e₁): integral boundary is from wmin to e₁
    for i = 1:ntime
        Λ[i] = ( 𝕊[i](e₂) - 𝕊[i](e₁) ) * r.h
    end

    return Λ
end

"""
    calc_error(G::Vector{F64}, Gᵥ::Vector{F64}, σ¹::Vector{F64})

Try to calculate χ². Here `Gᵥ` and `σ¹` denote the raw correlator and
related standard deviation. `G` means the reproduced correlator.

### Arguments
See above explanations.

### Returns
* Δ -> χ², distance between reconstructed and raw correlators.

See also: [`calc_green`](@ref).
"""
function calc_error(G::Vector{F64}, Gᵥ::Vector{F64}, σ¹::Vector{F64})
    return sum( ( (G .- Gᵥ) .* σ¹ ) .^ 2.0 )
end

"""
    calc_green(Λ::Array{F64,2}, nk::I64)

Try to reconstruct the correlator via the field configuration. Now this
function is called by init_element(). But perhaps we can use it in last().

### Arguments
* Λ -> The Λ function. See above remarks.
* nk -> Current number of boxes.

### Returns
* G -> Reconstructed Green's function.

See also: [`calc_error`](@ref).
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

Calculate the total area of all boxes. Now this function is not used.

### Arguments
* C -> The current Monte Carlo field configuration.

### Returns
* norm -> Area of all boxes.
"""
function calc_norm(C::Vector{Box})
    norm = sum(map(x -> x.h * x.w, C))
    return norm
end

"""
    constraints(e₁::F64, e₂::F64)

This function is used to judge whether a given box overlapes with the
forbidden zone. Here `e₁` and `e₂` denote the left and right boundaries
of the box.

### Arguments
See above explanations.

### Returns
* ex -> Boolean, whether a given box is valid.

"""
function constraints(e₁::F64, e₂::F64)
    exclude = get_b("exclude")
    @assert e₁ ≤ e₂

    if !isa(exclude, Missing)
        for i in eachindex(exclude)
            if e₁ ≤ exclude[i][1] ≤ e₂ ≤ exclude[i][2]
                return false
            end
            if exclude[i][1] ≤ e₁ ≤ exclude[i][2] ≤ e₂
                return false
            end
            if exclude[i][1] ≤ e₁ ≤ e₂ ≤ exclude[i][2]
                return false
            end
            if e₁ ≤ exclude[i][1] ≤ exclude[i][2] ≤ e₂
                return false
            end
        end
    end

    return true
end

"""
    try_insert(
        MC::StochOMMC,
        SE::StochOMElement,
        SC::StochOMContext,
        dacc::F64
    )

Insert a new box into the field configuration.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.
* dacc -> A predefined parameter used to calculate transition probability.

### Returns
N/A
"""
function try_insert(
    MC::StochOMMC,
    SE::StochOMElement,
    SC::StochOMContext,
    dacc::F64
    )
    sbox = get_s("sbox")
    wbox = get_s("wbox")
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    csize = length(SE.C)

    # Choose a box randomly
    t = rand(MC.rng, 1:csize)

    # Check area of this box
    R = SE.C[t]
    if R.h * R.w ≤ 2.0 * sbox
        return
    end

    # Determine minimum and maximum areas of the new box
    dx_min = sbox
    dx_max = R.h * R.w - sbox
    if dx_max ≤ dx_min
        return
    end

    # Determine parameters for the new box
    r₁ = rand(MC.rng, F64)
    r₂ = rand(MC.rng, F64)
    #
    c = (wmin + wbox / 2.0) + (wmax - wmin - wbox) * r₁
    #
    w_new_max = 2.0 * min(wmax - c, c - wmin)
    dx = Pdx(dx_min, dx_max, MC.rng)
    #
    h = dx / w_new_max + (dx / wbox - dx / w_new_max) * r₂
    w = dx / h

    # Rnew will be used to update Box t, while Radd is the new box.
    if !constraints(c - w/2, c + w/2)
        return
    end
    Rnew = Box(R.h - dx / R.w, R.w, R.c)
    Radd = Box(h, w, c)

    # Calculate update for Λ
    G₁ = SE.Λ[:,t]
    G₂ = eval_lambda(Rnew, SC.grid, SC.𝕊ᵥ)
    G₃ = eval_lambda(Radd, SC.grid, SC.𝕊ᵥ)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G₁ + G₂ + G₃, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t
        SE.C[t] = Rnew

        # Add new Box
        push!(SE.C, Radd)

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G₁ + G₂ + G₃
        @. SE.Λ[:,t] = G₂
        @. SE.Λ[:,csize+1] = G₃

        # Update the counter
        MC.Macc[1] = MC.Macc[1] + 1
    end

    # Update the counter
    MC.Mtry[1] = MC.Mtry[1] + 1
end

"""
    try_remove(
        MC::StochOMMC,
        SE::StochOMElement,
        SC::StochOMContext,
        dacc::F64
    )

Remove an old box from the field configuration.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.
* dacc -> A predefined parameter used to calculate transition probability.

### Returns
N/A
"""
function try_remove(
    MC::StochOMMC,
    SE::StochOMElement,
    SC::StochOMContext,
    dacc::F64
    )
    csize = length(SE.C)

    # Choose two boxes randomly
    # Box t₁ will be removed, while box t₂ will be modified.
    t₁ = rand(MC.rng, 1:csize)
    t₂ = rand(MC.rng, 1:csize)
    #
    while t₁ == t₂
        t₂ = rand(MC.rng, 1:csize)
    end
    #
    if t₁ < t₂
        t₁, t₂ = t₂, t₁
    end

    # Get box t₁ and box t₂
    R₁ = SE.C[t₁]
    R₂ = SE.C[t₂]
    Rₑ = SE.C[end]

    # Generate new box t₂
    dx = R₁.h * R₁.w
    R₂ₙ = Box(R₂.h + dx / R₂.w, R₂.w, R₂.c)

    # Calculate update for Λ
    G₁ = SE.Λ[:,t₁]
    G₂ = SE.Λ[:,t₂]
    Gₑ = SE.Λ[:,csize]
    G₂ₙ = eval_lambda(R₂ₙ, SC.grid, SC.𝕊ᵥ)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G₁ - G₂ + G₂ₙ, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t₂
        SE.C[t₂] = R₂ₙ

        # Backup the last box in box t₁
        if t₁ < csize
            SE.C[t₁] = Rₑ
        end

        # Delete the last box, since its value has been stored in t₁.
        pop!(SE.C)

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G₁ - G₂ + G₂ₙ
        @. SE.Λ[:,t₂] = G₂ₙ
        if t₁ < csize
            @. SE.Λ[:,t₁] = Gₑ
        end

        # Update the counter
        MC.Macc[2] = MC.Macc[2] + 1
    end

    # Update the counter
    MC.Mtry[2] = MC.Mtry[2] + 1
end

"""
    try_shift(
        MC::StochOMMC,
        SE::StochOMElement,
        SC::StochOMContext,
        dacc::F64
    )

Change the position of given box in the field configuration.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.
* dacc -> A predefined parameter used to calculate transition probability.

### Returns
N/A
"""
function try_shift(
    MC::StochOMMC,
    SE::StochOMElement,
    SC::StochOMContext,
    dacc::F64
    )
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    csize = length(SE.C)

    # Choose a box randomly
    t = rand(MC.rng, 1:csize)

    # Retreive the box t
    R = SE.C[t]

    # Determine left and right boundaries for the center of the box
    dx_min = wmin + R.w / 2.0 - R.c
    dx_max = wmax - R.w / 2.0 - R.c
    if dx_max ≤ dx_min
        return
    end

    # Calculate δc and generate shifted box
    δc = Pdx(dx_min, dx_max, MC.rng)
    if !constraints(R.c + δc - R.w/2, R.c + δc + R.w/2)
        return
    end
    Rₙ = Box(R.h, R.w, R.c + δc)

    # Calculate update for Λ
    G₁ = SE.Λ[:,t]
    G₂ = eval_lambda(Rₙ, SC.grid, SC.𝕊ᵥ)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G₁ + G₂, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t
        SE.C[t] = Rₙ

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G₁ + G₂
        @. SE.Λ[:,t] = G₂

        # Update the counter
        MC.Macc[3] = MC.Macc[3] + 1
    end

    # Update the counter
    MC.Mtry[3] = MC.Mtry[3] + 1
end

"""
    try_width(
        MC::StochOMMC,
        SE::StochOMElement,
        SC::StochOMContext,
        dacc::F64
    )

Change the width and height of given box in the field configuration. Note
that the box's area is kept.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.
* dacc -> A predefined parameter used to calculate transition probability.

### Returns
N/A
"""
function try_width(
    MC::StochOMMC,
    SE::StochOMElement,
    SC::StochOMContext,
    dacc::F64
    )
    wbox = get_s("wbox")
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    csize = length(SE.C)

    # Choose a box randomly
    t = rand(MC.rng, 1:csize)

    # Retreive the box t
    R = SE.C[t]

    # Determine left and right boundaries for the width of the box
    weight = R.h * R.w
    dx_min = wbox - R.w
    dx_max = min(2.0 * (R.c - wmin), 2.0 * (wmax - R.c)) - R.w
    if dx_max ≤ dx_min
        return
    end

    # Calculate δw and generate new box
    dw = Pdx(dx_min, dx_max, MC.rng)
    w = R.w + dw
    h = weight / w
    c = R.c
    if !constraints(c - w/2, c + w/2)
        return
    end
    Rₙ = Box(h, w, c)

    # Calculate update for Λ
    G₁ = SE.Λ[:,t]
    G₂ = eval_lambda(Rₙ, SC.grid, SC.𝕊ᵥ)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G₁ + G₂, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t
        SE.C[t] = Rₙ

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G₁ + G₂
        @. SE.Λ[:,t] = G₂

        # Update the counter
        MC.Macc[4] = MC.Macc[4] + 1
    end

    # Update the counter
    MC.Mtry[4] = MC.Mtry[4] + 1
end

"""
    try_height(
        MC::StochOMMC,
        SE::StochOMElement,
        SC::StochOMContext,
        dacc::F64
    )

Change the heights of two given boxes in the field configuration.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.
* dacc -> A predefined parameter used to calculate transition probability.

### Returns
N/A
"""
function try_height(
    MC::StochOMMC,
    SE::StochOMElement,
    SC::StochOMContext,
    dacc::F64
    )
    sbox  = get_s("sbox")
    csize = length(SE.C)

    # Choose two boxes randomly
    t₁ = rand(MC.rng, 1:csize)
    t₂ = rand(MC.rng, 1:csize)
    #
    while t₁ == t₂
        t₂ = rand(MC.rng, 1:csize)
    end

    # Get box t₁ and box t₂
    R₁ = SE.C[t₁]
    R₂ = SE.C[t₂]

    # Determine left and right boundaries for the height of the box t₁
    w₁ = R₁.w
    w₂ = R₂.w
    h₁ = R₁.h
    h₂ = R₂.h
    dx_min = sbox / w₁ - h₁
    dx_max = (h₂ - sbox / w₂) * w₂ / w₁
    if dx_max ≤ dx_min
        return
    end

    # Calculate δh and generate new box t₁ and box t₂
    dh = Pdx(dx_min, dx_max, MC.rng)
    R₁ₙ = Box(R₁.h + dh, R₁.w, R₁.c)
    R₂ₙ = Box(R₂.h - dh * w₁ / w₂, R₂.w, R₂.c)

    # Calculate update for Λ
    G₁A = SE.Λ[:,t₁]
    G₁B = eval_lambda(R₁ₙ, SC.grid, SC.𝕊ᵥ)
    G₂A = SE.Λ[:,t₂]
    G₂B = eval_lambda(R₂ₙ, SC.grid, SC.𝕊ᵥ)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G₁A + G₁B - G₂A + G₂B, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t₁ and box t₂
        SE.C[t₁] = R₁ₙ
        SE.C[t₂] = R₂ₙ

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G₁A + G₁B - G₂A + G₂B
        @. SE.Λ[:,t₁] = G₁B
        @. SE.Λ[:,t₂] = G₂B

        # Update the counter
        MC.Macc[5] = MC.Macc[5] + 1
    end

    # Update the counter
    MC.Mtry[5] = MC.Mtry[5] + 1
end

"""
    try_split(
        MC::StochOMMC,
        SE::StochOMElement,
        SC::StochOMContext,
        dacc::F64
    )

Split a given box into two boxes in the field configuration.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.
* dacc -> A predefined parameter used to calculate transition probability.

### Returns
N/A
"""
function try_split(
    MC::StochOMMC,
    SE::StochOMElement,
    SC::StochOMContext,
    dacc::F64
    )
    wbox = get_s("wbox")
    sbox = get_s("sbox")
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    csize = length(SE.C)

    # Choose a box randomly
    t = rand(MC.rng, 1:csize)

    # Retreive the box t
    R₁ = SE.C[t]
    if R₁.w ≤ 2 * wbox || R₁.w * R₁.h ≤ 2.0 * sbox
        return
    end

    # Determine height for new boxes (h and h)
    h = R₁.h

    # Determine width for new boxes (w₁ and w₂)
    w₁ = wbox + (R₁.w - 2.0 * wbox) * rand(MC.rng, F64)
    w₂ = R₁.w - w₁
    if w₁ > w₂
        w₁, w₂ = w₂, w₁
    end

    # Determine center for new boxes (c₁ + δc₁ and c₂ + δc₂)
    c₁ = R₁.c - R₁.w / 2.0 + w₁ / 2.0
    c₂ = R₁.c + R₁.w / 2.0 - w₂ / 2.0
    dx_min = wmin + w₁ / 2.0 - c₁
    dx_max = wmax - w₁ / 2.0 - c₁
    if dx_max ≤ dx_min
        return
    end
    δc₁ = Pdx(dx_min, dx_max, MC.rng)
    δc₂ = -1.0 * w₁ * δc₁ / w₂
    if !constraints(c₁ + δc₁ - w₁/2, c₁ + δc₁ + w₁/2) ||
       !constraints(c₂ + δc₂ - w₂/2, c₂ + δc₂ + w₂/2)
        return
    end

    if (c₁ + δc₁ ≥ wmin + w₁ / 2.0) &&
       (c₁ + δc₁ ≤ wmax - w₁ / 2.0) &&
       (c₂ + δc₂ ≥ wmin + w₂ / 2.0) &&
       (c₂ + δc₂ ≤ wmax - w₂ / 2.0)

        # Generate two new boxes
        R₂ = Box(h, w₁, c₁ + δc₁)
        R₃ = Box(h, w₂, c₂ + δc₂)

        # Calculate update for Λ
        G₁ = SE.Λ[:,t]
        G₂ = eval_lambda(R₂, SC.grid, SC.𝕊ᵥ)
        G₃ = eval_lambda(R₃, SC.grid, SC.𝕊ᵥ)

        # Calculate new Δ function, it is actually the error function.
        Δ = calc_error(SE.G - G₁ + G₂ + G₃, SC.Gᵥ, SC.σ¹)

        # Apply the Metropolis algorithm
        if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
            # Remove old box t and insert two new boxes
            SE.C[t] = R₂
            push!(SE.C, R₃)

            # Update Δ, G, and Λ.
            SE.Δ = Δ
            @. SE.G = SE.G - G₁ + G₂ + G₃
            @. SE.Λ[:,t] = G₂
            @. SE.Λ[:,csize+1] = G₃

            # Update the counter
            MC.Macc[6] = MC.Macc[6] + 1
        end
    end

    # Update the counter
    MC.Mtry[6] = MC.Mtry[6] + 1
end

"""
    try_merge(
        MC::StochOMMC,
        SE::StochOMElement,
        SC::StochOMContext,
        dacc::F64
    )

Merge two given boxes into one box in the field configuration.

### Arguments
* MC -> A StochOMMC struct. It containts some counters.
* SE -> A StochOMElement struct. It contains Monte Carlo configurations.
* SC -> A StochOMContext struct. It contains grid, mesh, and Gᵥ.
* dacc -> A predefined parameter used to calculate transition probability.

### Returns
N/A
"""
function try_merge(
    MC::StochOMMC,
    SE::StochOMElement,
    SC::StochOMContext,
    dacc::F64
    )
    wmin = get_b("wmin")
    wmax = get_b("wmax")
    csize = length(SE.C)

    # Choose two boxes randomly
    # Box t₂ will be removed, while box t₁ will be modified.
    t₁ = rand(MC.rng, 1:csize)
    t₂ = rand(MC.rng, 1:csize)
    #
    while t₁ == t₂
        t₂ = rand(MC.rng, 1:csize)
    end
    #
    if t₁ > t₂
        t₁, t₂ = t₂, t₁
    end

    # Get box t₁ and box t₂
    R₁ = SE.C[t₁]
    R₂ = SE.C[t₂]

    # Determine h, w, and c for new box
    weight = R₁.h * R₁.w + R₂.h * R₂.w
    w_new = 0.5 * (R₁.w + R₂.w)
    h_new = weight / w_new
    c_new = R₁.c + (R₂.c - R₁.c) * R₂.h * R₂.w / weight

    # Determine left and right boundaries for the center of the new box
    dx_min = wmin + w_new / 2.0 - c_new
    dx_max = wmax - w_new / 2.0 - c_new
    if dx_max ≤ dx_min
        return
    end

    # Calculate δc and generate new box
    δc = Pdx(dx_min, dx_max, MC.rng)
    if !constraints(c_new + δc - w_new/2, c_new + δc + w_new/2)
        return
    end
    Rₙ = Box(h_new, w_new, c_new + δc)

    # Calculate update for Λ
    G₁ = SE.Λ[:,t₁]
    G₂ = SE.Λ[:,t₂]
    Gₑ = SE.Λ[:,csize]
    Gₙ = eval_lambda(Rₙ, SC.grid, SC.𝕊ᵥ)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G₁ - G₂ + Gₙ, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t₁ with new box
        SE.C[t₁] = Rₙ

        # Delete box t₂
        if t₂ < csize
            SE.C[t₂] = SE.C[end]
        end
        pop!(SE.C)

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G₁ - G₂ + Gₙ
        @. SE.Λ[:,t₁] = Gₙ
        if t₂ < csize
            @. SE.Λ[:,t₂] = Gₑ
        end

        # Update the counter
        MC.Macc[7] = MC.Macc[7] + 1
    end

    # Update the counter
    MC.Mtry[7] = MC.Mtry[7] + 1
end

#=
*Remarks* : *Probability Density Function*

Every proposed elementary update is parametrized by a real number
``\delta\xi \in [\delta\xi_{min}:\delta\xi_{max}]``. A concrete meaning
of ``\delta\xi`` depends on the elementary update in question. For
instance, ``\delta\xi`` can be a shift of the centre of an existing
rectangle, or the weight of a rectangle to be added. In general, ``\delta
\xi`` are defined so that larger ``|\delta\xi|`` correspond to more
prominent changes in the configuration. `StochOM` randomly generates
values of ``\delta\xi`` according to the following probability density
function:

```math
\begin{equation}
\mathcal{P}(\delta\xi) = N
\exp\left(-\gamma \frac{|\delta\xi|}{X}\right),
\end{equation}
```

```math
\begin{equation}
X \equiv \max(|\delta\xi_{min}|, |\delta\xi_{max}|),
\end{equation}
```

```math
\begin{equation}
N = \frac{\gamma}{X}
\left[
\text{sign}(\delta\xi_{min})\left(e^{-\gamma|\delta\xi_{min}|/X}-1\right) +
\text{sign}(\delta\xi_{max})\left(1-e^{-\gamma|\delta\xi_{max}|/X}\right)
\right]^{-1}
\end{equation}
```

User can change parameter ``\gamma`` to control non-uniformity of the PDF.
=#

"""
    Pdx(xmin::F64, xmax::F64, rng::AbstractRNG)

Try to calculate the value of δξ for every elementary update according to
the probability density function. The actual meaning of δξ depends on the
elementary update.

### Arguments
* xmin -> Minimum value of δξ.
* xmax -> Maximum value of δξ
* rng -> Random number generator.

### Returns
* N -> Value of δξ.
"""
function Pdx(xmin::F64, xmax::F64, rng::AbstractRNG)
    xmin_abs = abs(xmin)
    xmax_abs = abs(xmax)
    X = max(xmin_abs, xmax_abs)

    γ = 2.0
    γ_X = γ / X

    η = rand(rng, F64)
    𝑁  = (1 - η) * copysign(expm1(-γ_X * xmin_abs), xmin)
    𝑁 +=      η  * copysign(expm1(-γ_X * xmax_abs), xmax)

    return copysign( log1p(-abs(𝑁)) / γ_X, 𝑁)
end
