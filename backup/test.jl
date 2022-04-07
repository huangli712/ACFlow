push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")
using ACFlow

g0, 𝐺, τ, Mrot = read_data!(ImaginaryTimeGrid)

sac_run(g0, 𝐺, τ, Mrot)
