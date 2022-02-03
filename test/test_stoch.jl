#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

println("Stochastic Analytical Continuation")
Gtau, Gdev, tmesh = StochAC.read_data()
MC, SE, SC = StochAC.stoch_init(tmesh, Gtau, Gdev)

StochAC.stoch_run(MC, SE, SC)
