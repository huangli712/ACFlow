#!/usr/bin/env julia

#
# Usage:
#
#     $ ppole.jl ac.toml
#

push!(LOAD_PATH, "/Users/lihuang/Working/devel/ACFlow/src/")

using Printf
using ACFlow

welcome()
overview()
read_param()
goodbye()
