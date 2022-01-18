#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/12/14
#

function calc_model(rfg::RealFrequencyGrid, 𝑀::MomentsData)
    #for i = 1:length(rfg.grid)
    #    @show i, rfg.grid[i]
    #end

    SC = 𝑀.𝑀₁ / 𝑀.𝑀₀
    σω = sqrt(𝑀.𝑀₂ / 𝑀.𝑀₀ - SC^2)

    model_center = 𝑀.𝑀₁ / 𝑀.𝑀₀
    model_shape = 2.0
    model_width = 2.0^(1.0/model_shape) * σω

    @show model_center
    @show model_shape
    @show model_width

    default_model = F64[]
    for i = 1:length(rfg.grid)
        m = π * model_shape / model_width / gamma(1.0 / model_shape)
        m = m * exp(-abs((rfg.grid[i] - model_center) / model_width)^model_shape)
        #@show i, m
        push!(default_model, m)
    end

    return default_model
end