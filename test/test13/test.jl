#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

welcome()
overview()
read_param()
solve(read_data(5))
