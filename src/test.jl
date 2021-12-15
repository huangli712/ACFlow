push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

println("Hello world! This is a maxent code.")
ω, 𝐺 = read_data!(FermionicMatsubaraGrid)
ωc, 𝑀 = calc_moments(ω, 𝐺)
trunc_data!(ωc, ω, 𝐺)
rfg = calc_mesh(𝑀)
default_model = calc_model(rfg, 𝑀)
𝐾, 𝐾𝑀 = calc_kernel(ω, rfg)

norm_DM_t = 𝐾𝑀.𝐾𝑀₀ * default_model
defalut_model = 𝑀.𝑀₀ * default_model ./ norm_DM_t[1,:]
@show default_model