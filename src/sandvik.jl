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

function sac_run(scale_factor::F64, 𝐺::GreenData, τ::ImaginaryTimeGrid, Mrot::AbstractMatrix)
    SG, SE, SC, MC, kernel = init_sac(scale_factor, 𝐺, τ, Mrot)
    anneal = perform_annealing(MC, SE, SC, SG, kernel, 𝐺)
    decide_sampling_theta(anneal, SC, SE, kernel, 𝐺)
    sample_and_collect(scale_factor, MC, SE, SC, SG, kernel, 𝐺)
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
