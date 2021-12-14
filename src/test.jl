push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

println("Hello world! This is a maxent code.")
ω = FermionicMatsubaraGrid()
𝐺 = GreenData()
read_data!(ω, 𝐺)
ωc, 𝑀, V𝑀 = calc_moments(ω, 𝐺)
trunc_data!(ωc, ω, 𝐺)
rfg = calc_mesh(𝑀)
default_model = calc_model(rfg, 𝑀)
𝐾, 𝐾𝑀 = calc_kernel(ω, rfg)