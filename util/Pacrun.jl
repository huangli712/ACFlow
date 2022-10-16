#!/usr/bin/env julia

using Distributed
addprocs(24)

@everywhere push!(LOAD_PATH, ENV["ACFLOW_HOME"])

@everywhere using ACFlow

welcome()
overview()
read_param()
solve(read_data())
