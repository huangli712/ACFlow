#
# Project : Gardenia
# Source  : global.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/10/01
#

#=
### *Global Constants* : *Numerical Types*
=#

"""
    I32

Alias of Integer type (32 bit).

See also: [`R32`](@ref), [`N32`](@ref).
"""
const I32 = Int32

"""
    I64

Alias of Integer type (64 bit).

See also: [`R64`](@ref), [`N64`](@ref).
"""
const I64 = Int64

"""
    API

Alias of Integer type (Arbitrary Precision Integer).

See also: [`APF`](@ref), [`APC`](@ref).
"""
const API = BigInt

"""
    F32

Alias of Float type (32 bit).

See also: [`R32`](@ref), [`N32`](@ref).
"""
const F32 = Float32

"""
    F64

Alias of Float type (64 bit).

See also: [`R64`](@ref), [`N64`](@ref).
"""
const F64 = Float64

"""
    APF

Alias of Float type (Arbitrary Precision Float).

See also: [`API`](@ref), [`APC`](@ref).
"""
const APF = BigFloat

"""
    C32

Alias of Complex type (32 bit).

See also: [`R32`](@ref), [`N32`](@ref).
"""
const C32 = ComplexF32

"""
    C64

Alias of Complex type (64 bit).

See also: [`R64`](@ref), [`N64`](@ref).
"""
const C64 = ComplexF64

"""
    APC

Alias of Complex type (Arbitrary Precision Complex).

See also: [`API`](@ref), [`APF`](@ref).
"""
const APC = Complex{BigFloat}

"""
    R32

Alias of Integer and Float types (32 bit).
Here `R` means Real.

See also: [`N32`](@ref), [`N64`](@ref), [`APN`](@ref).
"""
const R32 = Union{I32,F32}

"""
    R64

Alias of Integer and Float types (64 bit).
Here `R` means Real.

See also: [`N32`](@ref), [`N64`](@ref), [`APN`](@ref).
"""
const R64 = Union{I64,F64}

"""
    APR

Alias of Integer and Float types (Arbitrary Precision).
Here `R` means Real.

See also: [`N32`](@ref), [`N64`](@ref), [`APN`](@ref).
"""
const APR = Union{API,APF}

"""
    N32

Alias of Integer, Float, and Complex types (32 bit).
Here `N` means Number.

See also: [`R32`](@ref), [`R64`](@ref), [`APR`](@ref).
"""
const N32 = Union{I32,F32,C32}

"""
    N64

Alias of Integer, Float, and Complex types (64 bit).
Here `N` means Number.

See also: [`R32`](@ref), [`R64`](@ref), [`APR`](@ref).
"""
const N64 = Union{I64,F64,C64}

"""
    APN

Alias of Integer, Float, and Complex types (Arbitrary Precision).
Here `N` means Number.

See also: [`R32`](@ref), [`R64`](@ref), [`APR`](@ref).
"""
const APN = Union{API,APF,APC}

#=
### *Global Constants* : *Literal Strings*
=#

"""
    __LIBNAME__

Name of this julia toolkit.

See also: [`__VERSION__`](@ref).
"""
const __LIBNAME__ = "ACFlow"

"""
    __VERSION__

Version of this julia toolkit.

See also: [`__RELEASE__`](@ref).
"""
const __VERSION__ = v"2.1.1-devel.241001"

"""
    __RELEASE__

Release date of this julia toolkit.

See also: [`__AUTHORS__`](@ref).
"""
const __RELEASE__ = "2024/10"

#=
*Remarks* :

The elements of the Array `__AUTHORS__` should be a `NamedTuple` object,
such as:

```julia
(name = "author's name", email = "author's email")
```
=#

"""
    __AUTHORS__

Core authors of this julia toolkit.

See also: [`__LIBNAME__`](@ref).
"""
const __AUTHORS__ = [(name = "Li Huang", email = "huangli@caep.cn")]

"""
    authors()

Print authors / contributors of the `ACFlow` toolkit.

See also: [`__AUTHORS__`](@ref).
"""
function authors()
    println("Authors (Until $__RELEASE__):")
    for a in __AUTHORS__
        println("  $(a.name) (email: $(a.email))")
    end
end
