#!/usr/bin/env julia

#
# This script is used to generate non-uniform mesh according to the
# given spectral function. It will launch only 1 process.
#
# Usage:
#
#     $ gmesh.jl ac.toml
#

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using ACFlow

welcome()
overview()
read_param()

goodbye()
