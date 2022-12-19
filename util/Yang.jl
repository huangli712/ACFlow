#!/usr/bin/env julia
push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using ACFlow

function sample_p(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    if rand(MC.rng) < 0.9
        try_move_s(t, MC, SE, SC)
    else
        try_move_p(t, MC, SE, SC)
    end    
end

function sample_a(t::I64, MC::StochPXMC, SE::StochPXElement, SC::StochPXContext)
    if rand(MC.rng) < 0.5
        try_move_a(t, MC, SE, SC)
    else
        try_move_x(t, MC, SE, SC)
    end
end

function reset_element_p(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
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
    #
    selected = rand(rng, 1:npole, nselect)
    unique!(selected)
    nselect = length(selected)

    P = rand(rng, allow, nselect)
    @. SE.P[selected] = P
end

function reset_element_a(rng::AbstractRNG, allow::Vector{I64}, SE::StochPXElement)
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
    #
    selected = rand(rng, 1:npole, nselect)
    unique!(selected)
    nselect = length(selected)

    A₁ = SE.A[selected]
    s₁ = sum(A₁)
    #
    A₂ = rand(rng, F64, nselect)
    s₂ = sum(A₂)
    @. A₂ = A₂ / s₂ * s₁
    #
    @. SE.A[selected] = A₂
end

welcome()
overview()
read_param()