push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")
using ACFlow

#=
println("Hello world! This is a maxent code.")
ω, 𝐺 = read_data!(FermionicMatsubaraGrid)
ωc, 𝑀 = calc_moments(ω, 𝐺)
trunc_data!(ωc, ω, 𝐺)
rfg = calc_mesh(𝑀)
default_model = calc_model(rfg, 𝑀)
𝐾, 𝐾𝑀 = calc_kernel(ω, rfg)
diag_covar(𝑀, 𝐾𝑀)

norm_DM_t = 𝐾𝑀.𝐾𝑀₀ * default_model
defalut_model = 𝑀.𝑀₀ * default_model ./ norm_DM_t[1,:]
@show default_model
=#

#=
# For SOM PART

ω, 𝐺 = read_data!(FermionicMatsubaraGrid)
Aom = sum( pmap((x) -> som_run(ω, 𝐺), 1:nworkers()) ) ./ nworkers()
som_output(Aom)

=#

g0, 𝐺, τ, Mrot = read_data!(ImaginaryTimeGrid)
#error()
#grid = calc_grid()
#@show GridIndex2Freq(1, grid)
#@show GridIndex2Freq(grid.num_grid_index, grid)
#@show GridIndex2Freq(grid.num_grid_index-1, grid)
#@show Freq2GridIndex(-10.0, grid);
#@show Freq2GridIndex(9.99999, grid)
#@show Freq2GridIndex(9.99998, grid)
#@show SpecIndex2Freq(1, grid)
#@show SpecIndex2Freq(grid.num_spec_index, grid)
#@show SpecIndex2Freq(grid.num_spec_index-1, grid)

#@show Grid2Spec(1, grid)
#@show Grid2Spec(2, grid)
#@show Grid2Spec(2000, grid)

SG, SE, SC, MC, kernel = init_sac(g0, 𝐺, τ, Mrot)
#SE = init_spectrum(g0, grid, 𝐺, τ)
#kernel = init_kernel(τ, grid, Mrot)
#compute_corr_from_spec(kernel, SE, SC)
@show typeof(kernel)

sac_run(MC, SE, SC, SG, kernel, 𝐺)