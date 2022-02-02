include("stoch.jl")

using Distributed
using .StochFlow

println("Stochastic Analytical Continuation")
G, tmesh = read_data()
MC, SE, SG, SC = stoch_init(tmesh, G)

stoch_run(MC, SE, SC, SG, G)