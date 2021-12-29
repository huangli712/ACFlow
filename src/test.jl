push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

function som_output(Aom::Vector{F64})
    Ngrid = P_SOM["Ngrid"]
    ommin = P_SOM["ommin"]
    ommax = P_SOM["ommax"]

    open("Aw.data", "w") do fout
        for w = 1:Ngrid
            _omega = ommin + (w - 1) * (ommax - ommin) / (Ngrid - 1)
            println(fout, _omega, " ", Aom[w])
        end
    end
end

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

ω, 𝐺 = read_data!(FermionicMatsubaraGrid)
Aom = som_run(ω, 𝐺)
som_output(Aom)