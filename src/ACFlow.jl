#
# Project : Gardenia
# Source  : ACFlow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/05
#

module ACFlow

using Distributed
using LinearAlgebra
using Printf
using Dates
using TOML

include("global.jl")
include("types.jl")
include("util.jl")
include("mesh.jl")
include("grid.jl")

include("config.jl")
include("inout.jl")
include("kernel.jl")
include("model.jl")
include("flow.jl")

include("maxent.jl")
include("sac.jl")

export welcome
export overview
export read_param
export read_data
export solve

end
