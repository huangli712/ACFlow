#
# Project : Gardenia
# Source  : ACFlow.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2022/02/03
#

module ACFlow

using LinearAlgebra
using Printf
using TOML

include("global.jl")
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

export read_param
export read_data
export solve
export StochAC

end