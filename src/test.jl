push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

println("Hello world! This is a maxent code.")
ω = FermionicMatsubaraGrid()
𝐺 = GreenData()
read_data!(ω, 𝐺)
𝑀 = calc_moments(ω, 𝐺)
rfg = calc_mesh(𝑀)