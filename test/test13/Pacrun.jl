#!/usr/bin/env julia -p 8

@everywhere push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

@everywhere using ACFlow

welcome()
overview()
read_param()
solve(read_data(5))
