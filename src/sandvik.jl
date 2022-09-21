#
# Project : Gardenia
# Source  : sandvik.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/09/22
#



struct GreenData
    value :: Vector{N64}
    error :: Vector{N64}
    covar :: Vector{N64}
end

struct ImaginaryTimeGrid
    grid :: Vector{F64}
end

struct FermionicMatsubaraGrid
    grid :: Vector{F64}
end

function read_data!(::Type{FermionicMatsubaraGrid})
    grid  = F64[]
    value = C64[]
    error = C64[]
    covar = F64[]

    niw = 64
    #
    open("giw.data", "r") do fin
        for i = 1:niw
            arr = parse.(F64, line_to_array(fin))
            push!(grid, arr[1])
            push!(value, arr[2] + arr[3] * im)
        end
    end
    #
    open("err.data", "r") do fin
        for i = 1:niw
            arr = parse.(F64, line_to_array(fin))
            @assert grid[i] == arr[1]
            push!(error, arr[2] + arr[3] * im)
            push!(covar, arr[2]^2)
            push!(covar, arr[3]^2)
        end
    end

    return FermionicMatsubaraGrid(grid), GreenData(value, error, covar)
end

mutable struct SACMonteCarlo
    rng :: AbstractRNG
    acc :: F64
    sample_acc  :: Vector{F64}
    sample_chi2 :: Vector{F64}
    bin_acc :: Vector{F64}
    bin_chi2 :: Vector{F64}
end

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
    Θ :: F64
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

struct SACAnnealing
    Conf  :: Vector{SACElement}
    Theta :: Vector{F64}
    chi2  :: Vector{F64}
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

function SpecIndex2Freq(spec_index::I64, SG::SACGrid)
    @assert 1 ≤ spec_index ≤ SG.num_spec_index
    return SG.ommin + (spec_index - 1) * SG.spec_interval
end

function Freq2GridIndex(freq::F64, SG::SACGrid)
    @assert SG.ommin ≤ freq ≤ SG.ommax
    grid = ceil(I64, (freq - SG.ommin) / SG.grid_interval) + 1
    @assert 1 ≤ grid ≤ SG.num_grid_index
    return grid
end

function Grid2Spec(grid_index::I64, SG::SACGrid)
    @assert 1 ≤ grid_index ≤ SG.num_grid_index
    #@show (grid_index - 1) * SG.grid_interval / SG.spec_interval
    return ceil(I64, grid_index * SG.grid_interval / SG.spec_interval)
end

function init_sac(scale_factor::F64, 𝐺::GreenData, τ::ImaginaryTimeGrid, Mrot::AbstractMatrix)
    nbin = P_SAC["mc_bin_num"]
    sbin = P_SAC["mc_bin_size"]

    SG = calc_grid()

    ntau = length(τ.grid)
    Gr = Mrot * 𝐺.value
    #@show Gr
    #error()
    G1 = zeros(F64, ntau)
    G2 = zeros(F64, ntau)
    χ2 = 0.0
    χ2min = 0.0
    Θ = P_SAC["starting_theta"]
    freq = zeros(F64, SG.num_spec_index)
    spectrum = zeros(F64, SG.num_spec_index)
    SC = SACContext(Gr, G1, G2, χ2, χ2min, Θ, freq, spectrum)

    seed = rand(1:1000000)#;  seed = 840443
    rng = MersenneTwister(seed)
    @show "seed: ", seed
    acc = 0.0
    sample_acc = zeros(F64, sbin)
    sample_chi2 = zeros(F64, sbin)
    bin_acc = zeros(F64, nbin)
    bin_chi2 = zeros(F64, nbin)
    MC = SACMonteCarlo(rng, acc, sample_acc, sample_chi2, bin_acc, bin_chi2)

    SE = init_spectrum(scale_factor, SG, 𝐺, τ)
    #@show SE
    #error()

    kernel = init_kernel(τ, SG, Mrot)

    compute_corr_from_spec(kernel, SE, SC)

    χ = compute_goodness(SC.G1, SC.Gr, 𝐺.covar)
    SC.χ2 = χ
    SC.χ2min = χ
    #@show χ
    #error()
    #@show SG.num_spec_index

    return SG, SE, SC, MC, kernel
end

function sac_run(scale_factor::F64, 𝐺::GreenData, τ::ImaginaryTimeGrid, Mrot::AbstractMatrix)
    SG, SE, SC, MC, kernel = init_sac(scale_factor, 𝐺, τ, Mrot)
    anneal = perform_annealing(MC, SE, SC, SG, kernel, 𝐺)
    decide_sampling_theta(anneal, SC, SE, kernel, 𝐺)
    sample_and_collect(scale_factor, MC, SE, SC, SG, kernel, 𝐺)
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
    #@show size(τ.grid)
    #@show SG.num_grid_index
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
    #error()

    return kernel
end

function compute_corr_from_spec(kernel::AbstractMatrix, SE::SACElement, SC::SACContext)
    ndelta = P_SAC["ndelta"]
    #@show size(kernel)
    tmp_kernel = kernel[:, SE.C]
    #@show size(tmp_kernel), typeof(tmp_kernel)
    amplitude = fill(SE.A, ndelta)
    SC.G1 = tmp_kernel * amplitude
    #@show amplitude
    #@show SC.G1
    #error()
end

function compute_goodness(G::Vector{F64,}, Gr::Vector{F64}, Sigma::Vector{N64})
    #@show size(G), size(Gr), size(Sigma)
    χ = sum(((G .- Gr) .* Sigma) .^ 2.0)
    #@show χ
    return χ
end

function perform_annealing(MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::Matrix{F64}, 𝐺::GreenData)
    anneal_length = P_SAC["anneal_length"]
    #@show anneal_length

    #update_deltas_1step_single(MC, SE, SC, SG, kernel, 𝐺)

    Conf = SACElement[]
    Theta = F64[]
    Chi2 = F64[]

    for _ = 1:anneal_length
        update_fixed_theta(MC, SE, SC, SG, kernel, 𝐺)

        SC.χ2 = mean(MC.bin_chi2)

        push!(Conf, SE)
        push!(Theta, SC.Θ)
        push!(Chi2, SC.χ2)

        @show SC.χ2, SC.χ2min, SC.χ2 - SC.χ2min
        if SC.χ2 - SC.χ2min < 1e-3
            break
        end

        SC.Θ = SC.Θ / 1.1
    end

    return SACAnnealing(Conf, Theta, Chi2)
end

function update_fixed_theta(MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::Matrix{F64}, 𝐺::GreenData)
    nbin = P_SAC["mc_bin_num"]
    sbin = P_SAC["mc_bin_size"]
    ntau = length(𝐺.value)
    #@show nbin, sbin

    for n = 1:nbin
        for s = 1:sbin

            if s % 10 == 1
                SC.χ2 = compute_goodness(SC.G1, SC.Gr, 𝐺.covar)
            end

            update_deltas_1step_single(MC, SE, SC, SG, kernel, 𝐺)
            #@show n, s, SC.χ2, SC.χ2min

            MC.sample_chi2[s] = SC.χ2
            MC.sample_acc[s] = MC.acc
        end
        #error()

        MC.bin_chi2[n] = sum(MC.sample_chi2) / sbin
        MC.bin_acc[n] = sum(MC.sample_acc) / sbin

        # write log

        if MC.bin_acc[n] > 0.5
            r = SE.W * 1.5
            if ceil(I64, r) < SG.num_grid_index
                SE.W = ceil(I64, r)
            else
                SE.W = SG.num_grid_index
            end
        end

        if MC.bin_acc[n] < 0.4
            SE.W = ceil(I64, SE.W / 1.5)
        end

        @show n, MC.bin_chi2[n] / ntau, SC.χ2min / ntau, MC.bin_chi2[n] - SC.χ2min, SC.Θ
    end

    #error()
end

function update_deltas_1step_single(MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::Matrix{F64}, 𝐺::GreenData)
    ndelta = P_SAC["ndelta"]
    accept_count = 0.0

    #@show SE
    #@show SC.Θ
    #@show SG
    #@show SC.Gr

    #error()

    #@show SC.χ2
#    error()
    for i = 1:ndelta
        select_delta = rand(MC.rng, 1:ndelta)
#        select_delta = 384 # debug
        location_current = SE.C[select_delta]
#        @show location_current
        #error()

        if 1 < SE.W < SG.num_grid_index
            move_width = rand(MC.rng, 1:SE.W)
#            move_width = 5897

            if rand(MC.rng) > 0.5
                location_updated = location_current + move_width
            else
                location_updated = location_current - move_width
            end

            if location_updated < 1 || location_updated > SG.num_grid_index
                continue
            end

        elseif SE.W == SG.num_grid_index
            location_updated = rand(MC.rng, 1:SG.num_grid_index)
        else
            error("BIG PROBLEM")
        end

        SC.G2 = SC.G1 + SE.A .* (kernel[:,location_updated] .- kernel[:,location_current])
        #@show SC.G1
        #error()

        chi2_updated = compute_goodness(SC.G2, SC.Gr, 𝐺.covar)

        p = exp( (SC.χ2 - chi2_updated) / (2.0 * SC.Θ) )
        #@show p
        #error()

        if rand(MC.rng) < min(p, 1.0)
            SE.C[select_delta] = location_updated
            SC.G1 = deepcopy(SC.G2)
            SC.χ2 = chi2_updated
            if SC.χ2 < SC.χ2min
                SC.χ2min = SC.χ2
            end

            accept_count = accept_count + 1.0
        end
        #@show i, SC.χ2, SC.χ2min
    end

    MC.acc = accept_count / ndelta
    #@show MC.acc
    #error()
end

function decide_sampling_theta(anneal::SACAnnealing, SC::SACContext, SE::SACElement, kernel::AbstractMatrix, 𝐺::GreenData)
    num_anneal = length(anneal.chi2)

    c = num_anneal
    while c ≥ 1
        if anneal.chi2[c] > SC.χ2min + 2.0 * sqrt(SC.χ2min)
            break
        end
        c = c - 1
    end
    @assert 1 ≤ c ≤ num_anneal

    SE = deepcopy(anneal.Conf[c])
    SC.Θ = anneal.Theta[c]
    compute_corr_from_spec(kernel, SE, SC)
    SC.χ2 = compute_goodness(SC.G1, SC.Gr, 𝐺.covar)
end

function sample_and_collect(scale_factor::F64, MC::SACMonteCarlo, SE::SACElement, SC::SACContext, SG::SACGrid, kernel::AbstractMatrix, 𝐺::GreenData)
    ndelta = P_SAC["ndelta"]
    update_fixed_theta(MC, SE, SC, SG, kernel, 𝐺)

    for n = 1:SG.num_spec_index
        SC.freq[n] = SpecIndex2Freq(n, SG)
    end

    collecting_steps = 100000
    for i = 1:collecting_steps
        if i % 10 == 1
            SC.χ2 = compute_goodness(SC.G1, SC.Gr, 𝐺.covar)
            @show i, SC.χ2, SC.χ2min
        end

        # update
        update_deltas_1step_single(MC, SE, SC, SG, kernel, 𝐺)

        # collecting
        for j = 1:ndelta
            d_pos = SE.C[j]
            s_pos = Grid2Spec(d_pos, SG)
            #@show d_pos, s_pos
            SC.spectrum[s_pos] = SC.spectrum[s_pos] + SE.A
        end
    end

    factor = scale_factor / (collecting_steps * SG.spec_interval)
    SC.spectrum = SC.spectrum * factor

    open("SAC.spectrum", "w") do fout
        for i = 1:SG.num_spec_index
            println(fout, SC.freq[i], " ", SC.spectrum[i])
        end
    end
end
