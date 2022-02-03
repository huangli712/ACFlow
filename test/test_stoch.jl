#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

println("Stochastic Analytical Continuation")
G, tmesh = StochAC.read_data()
MC, SE, SG, SC = StochAC.stoch_init(tmesh, G)

StochAC.stoch_run(MC, SE, SC, SG, G)
