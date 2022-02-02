#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using Distributed
using StochFlow

println("Stochastic Analytical Continuation")
G, tmesh = read_data()
MC, SE, SG, SC = stoch_init(tmesh, G)

stoch_run(MC, SE, SC, SG, G)
