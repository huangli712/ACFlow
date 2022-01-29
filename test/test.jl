#!/usr/bin/env julia

push!(LOAD_PATH, "/Users/lihuang/Working/devel/acflow/src")

using ACFlow

println("hello world")

read_param()
raw = read_data()
solve(raw)
