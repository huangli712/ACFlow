#
# Project : Gardenia
# Source  : sac.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/05
#

const P_SAC = Dict{String,Any}(
    "ommax" => 10.0,
    "ommin" => -10.0,
    "grid_interval" => 1.0e-5,
    "spec_interval" => 1.0e-2,
    "ndelta" => 1000,
    "beta" => 4.0,
)

mutable struct SACElement
    C :: Vector{I64}
    A :: F64
    W :: I64
end

mutable struct SACContext
    Gr :: Vector{F64}
    G1 :: Vector{F64}
    G2 :: Vector{F64}
    χ2 :: F64
    χ2min :: F64
    acc :: F64
    freq :: Vector{F64}
    spectrum :: Vector{F64}
end

struct SACGrid
    ommax :: F64
    ommin :: F64
    grid_interval :: F64
    spec_interval :: F64
    num_grid_index :: I64
    num_spec_index :: I64
end

function calc_grid()
    ommax = P_SAC["ommax"]
    ommin = P_SAC["ommin"]
    grid_interval = P_SAC["grid_interval"]
    spec_interval = P_SAC["spec_interval"]
    num_grid_index = ceil(I64, (ommax - ommin) / grid_interval)
    num_spec_index = ceil(I64, (ommax - ommin) / spec_interval)
    #@show num_grid_index, num_spec_index
    return SACGrid(ommax, ommin, grid_interval, spec_interval, num_grid_index, num_spec_index)
end

function GridIndex2Freq(grid_index::I64, SG::SACGrid)
    @assert 1 ≤ grid_index ≤ SG.num_grid_index
    return SG.ommin + (grid_index - 1) * SG.grid_interval
end

function Freq2GridIndex(freq::F64, SG::SACGrid)
    @assert SG.ommin ≤ freq ≤ SG.ommax
    grid = ceil(I64, (freq - SG.ommin) / SG.grid_interval) + 1
    @assert 1 ≤ grid ≤ SG.num_grid_index
    return grid
end

function SpecIndex2Freq(spec_index::I64, SG::SACGrid)
    @assert 1 ≤ spec_index ≤ SG.num_spec_index
    return SG.ommin + (spec_index - 1) * SG.spec_interval
end

function Grid2Spec(grid_index::I64, SG::SACGrid)
    @assert 1 ≤ grid_index ≤ SG.num_grid_index
    #@show (grid_index - 1) * SG.grid_interval / SG.spec_interval
    return floor(I64, (grid_index - 1) * SG.grid_interval / SG.spec_interval)
end

function init_sac(scale_factor::F64, 𝐺::GreenData, τ::ImaginaryTimeGrid, Mrot::AbstractMatrix)
    SG = calc_grid()

    ntau = length(τ.grid)
    Gr = Mrot * 𝐺.value
    G1 = zeros(F64, ntau)
    G2 = zeros(F64, ntau)
    χ2 = 0.0
    χ2min = 0.0
    acc = 0.0
    freq = zeros(F64, SG.num_spec_index)
    spectrum = zeros(F64, SG.num_spec_index)
    SC = SACContext(Gr, G1, G2, χ2, χ2min, acc, freq, spectrum)

    SE = init_spectrum(scale_factor, SG, 𝐺, τ)

    kernel = init_kernel(τ, SG, Mrot)

    compute_corr_from_spec(kernel, SE, SC)

    χ = compute_goodness(SC.G1, SC.Gr, 𝐺.covar)
    SC.χ2 = χ
    SC.χ2min = χ
    @show SG.num_spec_index
end

function init_spectrum(scale_factor::F64, SG::SACGrid, 𝐺::GreenData, τ::ImaginaryTimeGrid)
    ndelta = P_SAC["ndelta"]

    left_edge = -2.0
    right_edge = 2.0
    left_edge_index = Freq2GridIndex(left_edge, SG)
    right_edge_index = Freq2GridIndex(right_edge, SG)
    rectangle_len = right_edge_index - left_edge_index + 1
    #@show left_edge_index, right_edge_index, rectangle_len

    position = I64[]
    for i = 1:ndelta
        push!(position, left_edge_index + (i - 1) % rectangle_len)
    end
    #@show position

    amplitude = 1.0 / (scale_factor * ndelta)
    #@show amplitude

    average_freq = abs(log(1.0/𝐺.value[end]) / τ.grid[end])
    #@show 𝐺.value[end]
    #@show τ.grid[end]
    #@show average_freq

    window_width = ceil(I64, 0.1 * average_freq / SG.grid_interval)
    #@show window_width

    return SACElement(position, amplitude, window_width)
end

function init_kernel(τ::ImaginaryTimeGrid, SG::SACGrid, Mrot::AbstractMatrix)
    @show size(τ.grid)
    @show SG.num_grid_index
    beta = P_SAC["beta"]

    ntau = length(τ.grid)
    nfreq = SG.num_grid_index
    kernel = zeros(F64, ntau, nfreq)

    for f = 1:nfreq
        ω = GridIndex2Freq(f, SG)
        de = 1.0 + exp(-beta * ω)
        #for t = 1:ntau
        #    kernel[t,f] = exp(-ω * τ.grid[t]) / de
        #end
        kernel[:,f] = exp.(-ω * τ.grid) / de
    end

    kernel = Mrot * kernel

    #for t = 1:ntau
    #    @show t, kernel[t,3]
    #end

    return kernel
end

function compute_corr_from_spec(kernel::AbstractMatrix, SE::SACElement, SC::SACContext)
    ndelta = P_SAC["ndelta"]
    @show size(kernel)

    tmp_kernel = kernel[:, SE.C]
    @show size(tmp_kernel), typeof(tmp_kernel)

    amplitude = fill(SE.A, ndelta)

    SC.G1 = tmp_kernel * amplitude
    #@show amplitude
    @show SC.G1
end

function compute_goodness(G::Vector{F64,}, Gr::Vector{F64}, Sigma::Vector{N64})
    #@show size(G), size(Gr), size(Sigma)

    χ = sum(((G .- Gr) .* Sigma) .^ 2.0)
    @show χ
end