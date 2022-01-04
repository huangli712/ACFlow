#
# Project : Gardenia
# Source  : sac.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2022/01/04
#

const P_SAC = Dict{String,Any}(
    "ommax" => 10.0,
    "ommin" => -10.0,
    "grid_interval" => 1.0e-5,
    "spec_interval" => 1.0e-2,
)

function calc_grid()
    ommax = P_SAC["ommax"]
    ommin = P_SAC["ommin"]
    grid_interval = P_SAC["grid_interval"]
    spec_interval = P_SAC["spec_interval"]

    num_grid_index = ceil((ommax - ommin) / grid_interval)
    num_spec_index = ceil((ommax - ommin) / spec_interval)
    @show num_grid_index, num_spec_index
end