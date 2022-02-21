#!/usr/bin/env julia -p 4

@everywhere push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

@everywhere using ACFlow

welcome()
overview()
read_param()
solve(read_data())
