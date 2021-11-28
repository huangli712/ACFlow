#
# Project : Gardenia
# Source  : util.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/28
#

"""
    sorry()

Print an error message to the screen.
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    line_to_array(io::IOStream)

Convert a line (reading from an IOStream) to a string array.
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array.
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
end
