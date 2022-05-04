#
# Project : Gardenia
# Source  : som.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/05/05
#

#=
### *Customized Structs* : *StochOM Solver*
=#

"""
    Box

Rectangle. The field configuration consists of many boxes. They exhibit
various areas (width × height). We used the Metropolis important sampling
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
be sampled by monte carlo procedure.

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
"""
mutable struct StochOMContext
    Gᵥ   :: Vector{F64}
    σ¹   :: Vector{F64}
    grid :: AbstractGrid
    mesh :: AbstractMesh
    Cᵥ   :: Vector{Vector{Box}}
    Δᵥ   :: Vector{F64}
end

#=
### *Global Drivers*
=#

"""
    solve(S::StochOMSolver, rd::RawData)

Solve the analytical continuation problem by the stochastic optimization
method.
"""
function solve(S::StochOMSolver, rd::RawData)
    println("[ StochOM ]")
    MC, SC = init(S, rd)

    if nworkers() > 1
        println("Using $(nworkers()) workers")
        #
        p1 = deepcopy(PCOMM)
        p2 = deepcopy(PStochOM)
        #
        sol = pmap((x) -> prun(S, p1, p2, MC, SC), 1:nworkers())
        @assert length(sol) == nworkers()
        #
        Aout = similar(sol[end])
        fill!(Aout, 0.0)
        for i in eachindex(sol)
            @. Aout = Aout + sol[i] / nworkers()
        end
        #
        Gout = last(SC, Aout)
    else
        Aout = run(S, MC, SC)
        Gout = last(SC, Aout)
    end

    return SC.mesh.mesh, Aout, Gout
end

"""
    init(S::StochOMSolver, rd::RawData)

Initialize the StochOM solver and return the StochOMMC and StochOMContext
structs.
"""
function init(S::StochOMSolver, rd::RawData)
    MC = init_mc(S)
    println("Create infrastructure for Monte Carlo sampling")

    Gᵥ, σ¹ = init_iodata(S, rd)
    println("Postprocess input data: ", length(σ¹), " points")

    grid = make_grid(rd)
    println("Build grid for input data: ", length(grid), " points")

    mesh = make_mesh()
    println("Build mesh for spectrum: ", length(mesh), " points")

    Cᵥ, Δᵥ = init_context(S)
    SC = StochOMContext(Gᵥ, σ¹, grid, mesh, Cᵥ, Δᵥ)

    return MC, SC
end

"""
    run(S::StochOMSolver, MC::StochOMMC, SC::StochOMContext)

Perform stochastic optimization simulation, sequential version.
"""
function run(S::StochOMSolver, MC::StochOMMC, SC::StochOMContext)
    ntry = get_s("ntry")
    nstep = get_s("nstep")

    for l = 1:ntry
        SE = init_element(MC, SC)

        for _ = 1:nstep
            update(MC, SE, SC)
        end

        SC.Δᵥ[l] = SE.Δ
        SC.Cᵥ[l] = deepcopy(SE.C)
        @printf("try -> %5i (%5i) Δ -> %8.4e \n", l, ntry, SE.Δ)
        flush(stdout)
        (l % 10 == 0) && write_statistics(MC)
    end

    return average(SC)
end

"""
    prun(S::StochOMSolver,
         p1::Dict{String,Vector{Any}},
         p2::Dict{String,Vector{Any}},
         MC::StochOMMC, SC::StochOMContext)

Perform stochastic optimization simulation, parallel version.
The arguments `p1` and `p2` are copies of PCOMM and PStochOM, respectively.
"""
function prun(S::StochOMSolver,
              p1::Dict{String,Vector{Any}},
              p2::Dict{String,Vector{Any}},
              MC::StochOMMC, SC::StochOMContext)
    rev_dict(p1)
    rev_dict(S, p2)

    MC.rng = MersenneTwister(rand(1:10000) * myid() + 1981)

    ntry = get_s("ntry")
    nstep = get_s("nstep")

    for l = 1:ntry
        SE = init_element(MC, SC)

        for _ = 1:nstep
            update(MC, SE, SC)
        end

        SC.Δᵥ[l] = SE.Δ
        SC.Cᵥ[l] = deepcopy(SE.C)
        @printf("try -> %5i (%5i) Δ -> %8.4e \n", l, ntry, SE.Δ)
        flush(stdout)
        (myid() == 2) && (l % 10 == 0) && write_statistics(MC)
    end

    return average(SC)
end

"""
    average(SC::StochOMContext)

Postprocess the collected results after the stochastic optimization
simulations. It will calculate real spectral functions.
"""
function average(SC::StochOMContext)
    nmesh = get_c("nmesh")
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
    for l = 1:ntry
        if SC.Δᵥ[l] < dev_ave / αgood
            for w = 1:nmesh
                _omega = SC.mesh[w]
                for r = 1:length(SC.Cᵥ[l])
                    R = SC.Cᵥ[l][r]
                    if R.c - 0.5 * R.w ≤ _omega ≤ R.c + 0.5 * R.w
                        Aom[w] = Aom[w] + R.h
                    end
                end
            end
        end
    end

    # Normalize the spectrum
    Lgood = count(x -> x < dev_ave / αgood, SC.Δᵥ)
    @. Aom = Aom / Lgood

    @printf("Median χ² : %16.12e Accepted configurations : %5i \n", dev_ave, Lgood)

    return Aom
end

"""
    last(SC::StochOMContext, Aout::Vector{F64})

It will process and write the calculated results by the StochOM solver,
including final spectral function and reproduced correlator.
"""
function last(SC::StochOMContext, Aout::Vector{F64})
    # Write the spectral function
    write_spectrum(SC.mesh, Aout)

    # Reproduce input data and write them
    kernel = make_kernel(SC.mesh, SC.grid)
    G = reprod(SC.mesh, kernel, Aout)
    write_backward(SC.grid, G)

    # Calculate full response function on real axis and write them
    _G = kramers(SC.mesh, Aout)
    write_complete(SC.mesh, _G)

    return _G
end

#=
### *Core Algorithms*
=#

"""
    update(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext)

Using the Metropolis algorithm to update the field configuration, i.e, a
collection of hundreds of boxes.
"""
function update(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext)
    Tmax = 100 # Length of the Markov chain
    nbox = get_s("nbox")

    T1 = rand(MC.rng, 1:Tmax)
    d1 = rand(MC.rng, F64)
    d2 = 1.0 + rand(MC.rng, F64)

    ST = deepcopy(SE)

    for _ = 1:T1
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < nbox
                    try_insert(MC, ST, SC, d1)
                end
                break

            @case 2
                if length(ST.C) > 1
                    try_remove(MC, ST, SC, d1)
                end
                break

            @case 3
                try_shift(MC, ST, SC, d1)
                break

            @case 4
                try_width(MC, ST, SC, d1)
                break

            @case 5
                if length(ST.C) > 1
                    try_height(MC, ST, SC, d1)
                end
                break

            @case 6
                if length(ST.C) < nbox
                    try_split(MC, ST, SC, d1)
                end
                break

            @case 7
                if length(ST.C) > 1
                    try_merge(MC, ST, SC, d1)
                end
                break
        end

    end

    for _ = T1+1:Tmax
        update_type = rand(MC.rng, 1:7)

        @cswitch update_type begin
            @case 1
                if length(ST.C) < nbox
                    try_insert(MC, ST, SC, d2)
                end
                break

            @case 2
                if length(ST.C) > 1
                    try_remove(MC, ST, SC, d2)
                end
                break

            @case 3
                try_shift(MC, ST, SC, d2)
                break

            @case 4
                try_width(MC, ST, SC, d2)
                break

            @case 5
                if length(ST.C) > 1
                    try_height(MC, ST, SC, d2)
                end
                break

            @case 6
                if length(ST.C) < nbox
                    try_split(MC, ST, SC, d2)
                end
                break

            @case 7
                if length(ST.C) > 1
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
    init_mc(S::StochOMSolver)

Try to create a StochOMMC struct.

See also: [`StochOMMC`](@ref).
"""
function init_mc(S::StochOMSolver)
    seed = rand(1:100000000)
    rng = MersenneTwister(seed)
    Macc = zeros(I64, 7)
    Mtry = zeros(I64, 7)

    MC = StochOMMC(rng, Macc, Mtry)

    return MC
end

"""
    init_element(MC::StochOMMC, SC::StochOMContext)

Try to initialize a StochOMElement struct.

See also: [`StochOMElement`](@ref).
"""
function init_element(MC::StochOMMC, SC::StochOMContext)
    wmin = get_c("wmin")
    wmax = get_c("wmax")
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
        Λ[:,k] .= calc_lambda(R, SC.grid)
    end
    #
    # Calculate green's function and relative error using boxes
    G = calc_green(Λ, _Know)
    Δ = calc_error(G, SC.Gᵥ, SC.σ¹)

    return StochOMElement(C, Λ, G, Δ)
end

"""
    init_context(S::StochOMSolver)

Try to initialize the key members of a StochOMContext struct.

See also: [`StochOMContext`](@ref).
"""
function init_context(S::StochOMSolver)
    ntry = get_s("ntry")
    nbox = get_s("nbox")

    Δv = zeros(F64, ntry)

    Cv = []
    for _ = 1:ntry
        C = Box[]
        for _ = 1:nbox
            push!(C, Box(0.0, 0.0, 0.0))
        end
        push!(Cv, C)
    end

    return Cv, Δv
end

"""
    init_iodata(S::StochACSolver, rd::RawData)

Preprocess the input data (`rd`).

See also: [`RawData`](@ref), [`GreenData`](@ref).
"""
function init_iodata(S::StochOMSolver, rd::RawData)
    G = make_data(rd)
    Gᵥ = G.value
    σ¹ = 1.0 ./ G.error

    return Gᵥ, σ¹
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

**C** For bosonic Matsubara frequency system (symmetric version).

```math
\begin{equation}
Λ_{R}(\omega_n) = h \int^{c+w/2}_{c-w/2}
    d\omega~\frac{\omega^2}{\omega_n^2 + \omega^2}.
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
Λ_{R}(\omega_n) = h
    \left[
        \omega - \omega_n \tan^{-1} \left(\frac{\omega}{\omega_n}\right)
    \right]
    \bigg|^{c+w/2}_{c-w/2}
\end{equation}
```

Finally, we have

```math
\begin{equation}
Λ_{R}(\omega_n) = h
    \left[
        w +
        \omega_n \tan^{-1} \left(\frac{c - w/2}{\omega_n}\right)
        -
        \omega_n \tan^{-1} \left(\frac{c + w/2}{\omega_n}\right)
    \right]
\end{equation}
```

We have implemented the above formulas in `calc_lambda()`.
=#

"""
    calc_lambda(r::Box, grid::FermionicMatsubaraGrid)

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for FermionicMatsubaraGrid only.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
function calc_lambda(r::Box, grid::FermionicMatsubaraGrid)
    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w
    iw = im * grid.ω
    Λ = @. r.h * log((iw - e₁) / (iw - e₂))
    return vcat(real(Λ), imag(Λ))
end

"""
    calc_lambda(r::Box, grid::FermionicImaginaryTimeGrid)

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for FermionicImaginaryTimeGrid only.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
function calc_lambda(r::Box, grid::FermionicImaginaryTimeGrid)
    sorry()
end

"""
    calc_lambda(r::Box, grid::BosonicMatsubaraGrid)

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for BosonicMatsubaraGrid only.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
function calc_lambda(r::Box, grid::BosonicMatsubaraGrid)
    ktype = get_c("ktype")

    e₁ = r.c - 0.5 * r.w
    e₂ = r.c + 0.5 * r.w

    if ktype == "bsymm"
        Λ = @. atan( e₁ / grid.ω ) - atan( e₂ / grid.ω )
        Λ = r.h * (r.w .+ grid.ω .* Λ)
        return Λ
    else
        iw = im * grid.ω
        Λ = @. r.h * (-r.w + iw * log((iw - e₁) / (iw - e₂)))
        return vcat(real(Λ), imag(Λ))
    end
end

"""
    calc_lambda(r::Box, grid::BosonicImaginaryTimeGrid)

Try to calculate the contribution of a given box `r` to the Λ function.
This function works for BosonicImaginaryTimeGrid only.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
function calc_lambda(r::Box, grid::BosonicImaginaryTimeGrid)
    ktype = get_c("ktype")

    if ktype == "bsymm"
        Λ = @. r.h * exp(-1 * grid.τ * r.c) * sinh(0.5 * grid.τ * r.w)
        Λ = Λ .* π ./ grid.τ
        return Λ
    else
        sorry()
    end
end

"""
    calc_error(G::Vector{F64}, Gᵥ::Vector{F64}, σ¹::Vector{F64})

Try to calculate χ². Here `Gᵥ` and `σ¹` denote the raw correlator and
related standard deviation. `G` means the reproduced correlator.
"""
function calc_error(G::Vector{F64}, Gᵥ::Vector{F64}, σ¹::Vector{F64})
    return sum( abs.((G .- Gᵥ) .* σ¹) )
end

"""
    calc_green(Λ::Array{F64,2}, nk::I64)

Try to reconstruct the correlator via the field configuration.
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

Calculate the total area of all boxes.
"""
function calc_norm(C::Vector{Box})
    norm = sum(map(x -> x.h * x.w, C))
    return norm
end

"""
    constraints(e₁::F64, e₂::F64)
"""
function constraints(e₁::F64, e₂::F64)
    exclude = get_c("exclude")
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
    try_insert(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Insert a new box into the field configuration.
"""
function try_insert(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    sbox = get_s("sbox")
    wbox = get_s("wbox")
    wmin = get_c("wmin")
    wmax = get_c("wmax")
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
    r1 = rand(MC.rng, F64)
    r2 = rand(MC.rng, F64)
    c = (wmin + wbox / 2.0) + (wmax - wmin - wbox) * r1
    w_new_max = 2.0 * min(wmax - c, c - wmin)
    dx = Pdx(dx_min, dx_max, MC.rng)
    h = dx / w_new_max + (dx / wbox - dx / w_new_max) * r2
    w = dx / h

    # Rnew will be used to update Box t, while Radd is the new box.
    if !constraints(c - w/2, c + w/2)
        return
    end
    Rnew = Box(R.h - dx / R.w, R.w, R.c)
    Radd = Box(h, w, c)

    # Calculate update for Λ
    G1 = SE.Λ[:,t]
    G2 = calc_lambda(Rnew, SC.grid)
    G3 = calc_lambda(Radd, SC.grid)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G1 + G2 + G3, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t
        SE.C[t] = Rnew

        # Add new Box
        push!(SE.C, Radd)

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G1 + G2 + G3
        @. SE.Λ[:,t] = G2
        @. SE.Λ[:,csize+1] = G3

        # Update the counter
        MC.Macc[1] = MC.Macc[1] + 1
    end

    # Update the counter
    MC.Mtry[1] = MC.Mtry[1] + 1
end

"""
    try_remove(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Remove an old box from the field configuration.
"""
function try_remove(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    csize = length(SE.C)

    # Choose two boxes randomly
    # Box t1 will be removed, while box t2 will be modified.
    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    #
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end
    #
    if t1 < t2
        t1, t2 = t2, t1
    end

    # Get box t1 and box t2
    R1 = SE.C[t1]
    R2 = SE.C[t2]
    Re = SE.C[end]

    # Generate new box t2
    dx = R1.h * R1.w
    R2n = Box(R2.h + dx / R2.w, R2.w, R2.c)

    # Calculate update for Λ
    G1 = SE.Λ[:,t1]
    G2 = SE.Λ[:,t2]
    Ge = SE.Λ[:,csize]
    G2n = calc_lambda(R2n, SC.grid)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G1 - G2 + G2n, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t2
        SE.C[t2] = R2n

        # Backup the last box in box t1
        if t1 < csize
            SE.C[t1] = Re
        end

        # Delete the last box, since its value has been stored in t1.
        pop!(SE.C)

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G1 - G2 + G2n
        @. SE.Λ[:,t2] = G2n
        if t1 < csize
            @. SE.Λ[:,t1] = Ge
        end

        # Update the counter
        MC.Macc[2] = MC.Macc[2] + 1
    end

    # Update the counter
    MC.Mtry[2] = MC.Mtry[2] + 1
end

"""
    try_shift(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Change the position of given box in the field configuration.
"""
function try_shift(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wmin = get_c("wmin")
    wmax = get_c("wmax")
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
    dc = Pdx(dx_min, dx_max, MC.rng)
    if !constraints(R.c + dc - R.w/2, R.c + dc + R.w/2)
        return
    end
    Rn = Box(R.h, R.w, R.c + dc)

    # Calculate update for Λ
    G1 = SE.Λ[:,t]
    G2 = calc_lambda(Rn, SC.grid)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G1 + G2, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t
        SE.C[t] = Rn

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G1 + G2
        @. SE.Λ[:,t] = G2

        # Update the counter
        MC.Macc[3] = MC.Macc[3] + 1
    end

    # Update the counter
    MC.Mtry[3] = MC.Mtry[3] + 1
end

"""
    try_width(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Change the width and height of given box in the field configuration. Note
that the box's area is kept.
"""
function try_width(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wbox = get_s("wbox")
    wmin = get_c("wmin")
    wmax = get_c("wmax")
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
    Rn = Box(h, w, c)

    # Calculate update for Λ
    G1 = SE.Λ[:,t]
    G2 = calc_lambda(Rn, SC.grid)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G1 + G2, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t
        SE.C[t] = Rn

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G1 + G2
        @. SE.Λ[:,t] = G2

        # Update the counter
        MC.Macc[4] = MC.Macc[4] + 1
    end

    # Update the counter
    MC.Mtry[4] = MC.Mtry[4] + 1
end

"""
    try_height(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Change the heights of two given boxes in the field configuration.
"""
function try_height(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    sbox  = get_s("sbox")
    csize = length(SE.C)

    # Choose two boxes randomly
    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end

    # Get box t1 and box t2
    R1 = SE.C[t1]
    R2 = SE.C[t2]

    # Determine left and right boundaries for the height of the box t1
    w1 = R1.w
    w2 = R2.w
    h1 = R1.h
    h2 = R2.h
    dx_min = sbox / w1 - h1
    dx_max = (h2 - sbox / w2) * w2 / w1
    if dx_max ≤ dx_min
        return
    end

    # Calculate δh and generate new box t1 and box t2
    dh = Pdx(dx_min, dx_max, MC.rng)
    R1n = Box(R1.h + dh, R1.w, R1.c)
    R2n = Box(R2.h - dh * w1 / w2, R2.w, R2.c)

    # Calculate update for Λ
    G1A = SE.Λ[:,t1]
    G1B = calc_lambda(R1n, SC.grid)
    G2A = SE.Λ[:,t2]
    G2B = calc_lambda(R2n, SC.grid)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G1A + G1B - G2A + G2B, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t1 and box t2
        SE.C[t1] = R1n
        SE.C[t2] = R2n

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G1A + G1B - G2A + G2B
        @. SE.Λ[:,t1] = G1B
        @. SE.Λ[:,t2] = G2B

        # Update the counter
        MC.Macc[5] = MC.Macc[5] + 1
    end

    # Update the counter
    MC.Mtry[5] = MC.Mtry[5] + 1
end

"""
    try_split(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Split a given box into two boxes in the field configuration.
"""
function try_split(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wbox = get_s("wbox")
    sbox = get_s("sbox")
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    csize = length(SE.C)

    # Choose a box randomly
    t = rand(MC.rng, 1:csize)

    # Retreive the box t
    R1 = SE.C[t]
    if R1.w ≤ 2 * wbox || R1.w * R1.h ≤ 2.0 * sbox
        return
    end

    # Determine height for new boxes (h and h)
    h = R1.h

    # Determine width for new boxes (w1 and w2)
    w1 = wbox + (R1.w - 2.0 * wbox) * rand(MC.rng, F64)
    w2 = R1.w - w1
    if w1 > w2
        w1, w2 = w2, w1
    end

    # Determine center for new boxes (c1 + dc1 and c2 + dc2)
    c1 = R1.c - R1.w / 2.0 + w1 / 2.0
    c2 = R1.c + R1.w / 2.0 - w2 / 2.0
    dx_min = wmin + w1 / 2.0 - c1
    dx_max = wmax - w1 / 2.0 - c1
    if dx_max ≤ dx_min
        return
    end
    dc1 = Pdx(dx_min, dx_max, MC.rng)
    dc2 = -1.0 * w1 * dc1 / w2
    if !constraints(c1 + dc1 - w1/2, c1 + dc1 + w1/2) ||
       !constraints(c2 + dc2 - w2/2, c2 + dc2 + w2/2)
        return
    end

    if (c1 + dc1 ≥ wmin + w1 / 2.0) &&
       (c1 + dc1 ≤ wmax - w1 / 2.0) &&
       (c2 + dc2 ≥ wmin + w2 / 2.0) &&
       (c2 + dc2 ≤ wmax - w2 / 2.0)

        # Generate two new boxes
        R2 = Box(h, w1, c1 + dc1)
        R3 = Box(h, w2, c2 + dc2)

        # Calculate update for Λ
        G1 = SE.Λ[:,t]
        Ge = SE.Λ[:,csize]
        G2 = calc_lambda(R2, SC.grid)
        G3 = calc_lambda(R3, SC.grid)

        # Calculate new Δ function, it is actually the error function.
        Δ = calc_error(SE.G - G1 + G2 + G3, SC.Gᵥ, SC.σ¹)

        # Apply the Metropolis algorithm
        if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
            # Remove old box t and insert two new boxes
            SE.C[t] = SE.C[end]
            pop!(SE.C)
            push!(SE.C, R2)
            push!(SE.C, R3)

            # Update Δ, G, and Λ.
            SE.Δ = Δ
            @. SE.G = SE.G - G1 + G2 + G3
            if t < csize
                @. SE.Λ[:,t] = Ge
            end
            @. SE.Λ[:,csize] = G2
            @. SE.Λ[:,csize+1] = G3

            # Update the counter
            MC.Macc[6] = MC.Macc[6] + 1
        end
    end

    # Update the counter
    MC.Mtry[6] = MC.Mtry[6] + 1
end

"""
    try_merge(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)

Merge two given boxes into one box in the field configuration.
"""
function try_merge(MC::StochOMMC, SE::StochOMElement, SC::StochOMContext, dacc::F64)
    wmin = get_c("wmin")
    wmax = get_c("wmax")
    csize = length(SE.C)

    # Choose two boxes randomly
    t1 = rand(MC.rng, 1:csize)
    t2 = rand(MC.rng, 1:csize)
    while t1 == t2
        t2 = rand(MC.rng, 1:csize)
    end
    if t1 > t2
        t1, t2 = t2, t1
    end

    # Get box t1 and box t2
    R1 = SE.C[t1]
    R2 = SE.C[t2]

    # Determine h, w, and c for new box
    weight = R1.h * R1.w + R2.h * R2.w
    w_new = 0.5 * (R1.w + R2.w)
    h_new = weight / w_new
    c_new = R1.c + (R2.c - R1.c) * R2.h * R2.w / weight

    # Determine left and right boundaries for the center of the new box
    dx_min = wmin + w_new / 2.0 - c_new
    dx_max = wmax - w_new / 2.0 - c_new
    if dx_max ≤ dx_min
        return
    end

    # Calculate δc and generate new box
    dc = Pdx(dx_min, dx_max, MC.rng)
    if !constraints(c_new + dc - w_new/2, c_new + dc + w_new/2)
        return
    end
    Rn = Box(h_new, w_new, c_new + dc)

    # Calculate update for Λ
    G1 = SE.Λ[:,t1]
    G2 = SE.Λ[:,t2]
    Ge = SE.Λ[:,csize]
    Gn = calc_lambda(Rn, SC.grid)

    # Calculate new Δ function, it is actually the error function.
    Δ = calc_error(SE.G - G1 - G2 + Gn, SC.Gᵥ, SC.σ¹)

    # Apply the Metropolis algorithm
    if rand(MC.rng, F64) < ((SE.Δ/Δ) ^ (1.0 + dacc))
        # Update box t1 with new box
        SE.C[t1] = Rn

        # Delete box t2
        if t2 < csize
            SE.C[t2] = SE.C[end]
        end
        pop!(SE.C)

        # Update Δ, G, and Λ.
        SE.Δ = Δ
        @. SE.G = SE.G - G1 - G2 + Gn
        @. SE.Λ[:,t1] = Gn
        if t2 < csize
            @. SE.Λ[:,t2] = Ge
        end

        # Update the counter
        MC.Macc[7] = MC.Macc[7] + 1
    end

    # Update the counter
    MC.Mtry[7] = MC.Mtry[7] + 1
end

#=
*Remarks* : *Probability Density Function*

Every Every proposed elementary update is parametrized by a real number
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
=#

"""
    Pdx(xmin::F64, xmax::F64, rng::AbstractRNG)

Try to calculate the probability density function.
"""
function Pdx(xmin::F64, xmax::F64, rng::AbstractRNG)
    γ = 2.0
    y = rand(rng, F64)

    _X = max(abs(xmin), abs(xmax))
    _λ = γ / _X
    _elx = exp(-1.0 * _λ * abs(xmin))
    _N = _λ / ( (xmin / abs(xmin)) * (exp(-1.0 * _λ * abs(xmin)) - 1.0)
              + (xmax / abs(xmax)) * (1.0 - exp(-1.0 * _λ * abs(xmax))) )
    _lysn = _λ * y / _N

    if xmin ≥ 0
        return -1.0 * log(_elx - _lysn) / _λ
    elseif xmax ≤ 0
        return log(_lysn + _elx) / _λ
    else
        _C1 = _N * (1.0 - _elx) / _λ
        if y ≤ _C1
            return log(_lysn + _elx) / _λ
        else
            return -1.0 * log(1.0 - _lysn + _λ * _C1 / _N) / _λ
        end
    end
end
