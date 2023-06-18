#!/usr/bin/env julia

#
# This script is used to examine the constrained stochastic pole
# expansion method. It will try to fix the positions or amplitudes
# of the poles, and then optimize the rest ones. It will launch
# only 1 process. Note that this script does not support matrix-valued
# green's function.
#
# Usage:
#
#     $ c_spx.jl ac.toml
#

using Random
using Printf
using ACFlow

"""
    sample_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Try to sample the positions of poles only.
"""
function sample_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    if rand(MC.rng) < 0.9
        try_move_s(t, MC, SE, SC)
    else
        try_move_p(t, MC, SE, SC)
    end
end

"""
    sample_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Try to sample the amplitudes of poles only.
"""
function sample_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    if rand(MC.rng) < 0.5
        try_move_a(t, MC, SE, SC)
    else
        try_move_x(t, MC, SE, SC)
    end
end

"""
    reset_element_p(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)

Reset the positions of poles.
"""
function reset_element_p(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
    npole = get_x("npole")

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
    P = rand(rng, allow, nselect)
    @. SE.P[selected] = P
end

"""
    reset_element_a(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)

Reset the amplitudes of poles.
"""
function reset_element_a(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
    npole = get_x("npole")

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

    # Change poles' amplitudes
    A₁ = SE.A[selected]
    s₁ = sum(A₁)
    #
    A₂ = rand(rng, F64, nselect)
    s₂ = sum(A₂)
    @. A₂ = A₂ / s₂ * s₁
    #
    @. SE.A[selected] = A₂
end

"""
    run_p(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Perform stochastic pole expansion simulation, sequential version.
Only the positions of poles are updated.
"""
function run_p(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # By default, we should write the analytical continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Setup essential parameters
    ntry = get_x("ntry")
    nstep = get_x("nstep")

    # Warmup the Monte Carlo engine
    println("Start thermalization...")
    for _ = 1:nstep
        sample_p(1, MC, SE, SC)
    end

    # Sample and collect data
    println("Start stochastic sampling...")
    for t = 1:ntry
        # Reset Monte Carlo counters
        reset_mc(MC)

        # Reset Monte Carlo field configuration
        reset_element_p(MC.rng, SC.allow, SE)

        # Reset Gᵧ and χ²
        reset_context(t, SE, SC)

        # Apply simulated annealing algorithm
        for _ = 1:nstep
            sample_p(t, MC, SE, SC)
        end

        # Write Monte Carlo statistics
        fwrite && write_statistics(MC)

        # Update χ²[t] to be consistent with SC.Pᵥ[t], SC.Aᵥ[t], and SC.𝕊ᵥ[t].
        SC.χ²[t] = SC.χ²min
        @printf("try = %6i -> [χ² = %9.4e]\n", t, SC.χ²min)
        flush(stdout)
    end

    # Write pole expansion coefficients
    fwrite && write_pole(SC.Pᵥ, SC.Aᵥ, SC.𝕊ᵥ, SC.χ², SC.fmesh)

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    run_a(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)

Perform stochastic pole expansion simulation, sequential version.
Only the amplitudes of poles are updated.
"""
function run_a(MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    # By default, we should write the analytical continuation results
    # into the external files.
    _fwrite = get_b("fwrite")
    fwrite = isa(_fwrite, Missing) || _fwrite ? true : false

    # Setup essential parameters
    ntry = get_x("ntry")
    nstep = get_x("nstep")

    # Warmup the Monte Carlo engine
    println("Start thermalization...")
    for _ = 1:nstep
        sample_a(1, MC, SE, SC)
    end

    # Sample and collect data
    println("Start stochastic sampling...")
    for t = 1:ntry
        # Reset Monte Carlo counters
        reset_mc(MC)

        # Reset Monte Carlo field configuration
        reset_element_a(MC.rng, SC.allow, SE)

        # Reset Gᵧ and χ²
        reset_context(t, SE, SC)

        # Apply simulated annealing algorithm
        for _ = 1:nstep
            sample_a(t, MC, SE, SC)
        end

        # Write Monte Carlo statistics
        fwrite && write_statistics(MC)

        # Update χ²[t] to be consistent with SC.Pᵥ[t], SC.Aᵥ[t], and SC.𝕊ᵥ[t].
        SC.χ²[t] = SC.χ²min
        @printf("try = %6i -> [χ² = %9.4e]\n", t, SC.χ²min)
        flush(stdout)
    end

    # Write pole expansion coefficients
    fwrite && write_pole(SC.Pᵥ, SC.Aᵥ, SC.𝕊ᵥ, SC.χ², SC.fmesh)

    # Generate spectral density from Monte Carlo field configuration
    return average(SC)
end

"""
    solve_p()

    Solve the analytical continuation problem by the stochastic
    pole expansion. Note that the amplitudes of poles are fixed.
    The positions of poles are optimized.
"""
function solve_p()
    npole = get_x("npole")

    S = StochPXSolver()
    rd = read_data()

    println("[ StochPX ]")
    MC, SE, SC = init(S, rd)

    # Setup the amplitudes of poles
    @. SE.A = 1.0 / npole
    reset_context(1, SE, SC)

    Aout, Gout, Gᵣ = run_p(MC, SE, SC)
    ACFlow.last(SC, Aout, Gout, Gᵣ)
end

"""
    solve_p()

    Solve the analytical continuation problem by the stochastic
    pole expansion. Note that the positions of poles are fixed.
    The amplitudes of poles are optimized.
"""
function solve_a()
    # to be done
end

welcome()
overview()
read_param()
solve_p()
