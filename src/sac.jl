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

struct SACGrid
    ommax :: F64
    ommin :: F64
    grid_interval :: F64
    spec_interval :: F64
    num_grid_index :: I64
    num_spec_index :: I64
end

function calc_grid()
    ommax = P_SAC["ommax"]
    ommin = P_SAC["ommin"]
    grid_interval = P_SAC["grid_interval"]
    spec_interval = P_SAC["spec_interval"]

    num_grid_index = ceil(I64, (ommax - ommin) / grid_interval)
    num_spec_index = ceil(I64, (ommax - ommin) / spec_interval)
    #@show num_grid_index, num_spec_index
    return SACGrid(ommax, ommin, grid_interval, spec_interval, num_grid_index, num_spec_index)
end

function GridIndex2Freq(grid_index::I64, SG::SACGrid)
    @assert 1 ≤ grid_index ≤ SG.num_grid_index
    return SG.ommin + (grid_index - 1) * SG.grid_interval
end

function Freq2GridIndex(freq::F64, SG::SACGrid)
    @assert SG.ommin ≤ freq ≤ SG.ommax
    grid = ceil(I64, (freq - SG.ommin) / SG.grid_interval) + 1
    @assert 1 ≤ grid ≤ SG.num_grid_index
end

function SpecIndex2Freq(spec_index::I64, SG::SACGrid)
    @assert 1 ≤ spec_index ≤ SG.num_spec_index
    return SG.ommin + (spec_index - 1) * SG.spec_interval
end

function Grid2Spec(grid_index::I64, SG::SACGrid)
    @assert 1 ≤ grid_index ≤ SG.num_grid_index
    return floor(I64, grid_index * SG.grid_interval / SG.spec_interval)
end
