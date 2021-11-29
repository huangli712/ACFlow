#
# Project : Gardenia
# Source  : model.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/29
#

function calc_model(rfg::RealFrequencyGrid, 𝑀::MomentsData)
    for i = 1:length(rfg.grid)
        @show i, rfg.grid[i]
    end
end