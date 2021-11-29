#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/30
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
end