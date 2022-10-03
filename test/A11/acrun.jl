#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])

using ACFlow

welcome()
overview()
read_param()
#solve(read_data())
san_run()
