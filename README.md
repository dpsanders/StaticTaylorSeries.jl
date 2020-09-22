# StaticTaylorSeries.jl
A [Julia](http://julialang.org) package used to compute static 1D Taylor
polynomial expansions

[![Build Status](https://api.travis-ci.org/dpsanders/TaylorSeries.jl.svg?branch=master)](https://travis-ci.org/JuliaDiff/TaylorSeries.jl)
[![Coverage Status](https://coveralls.io/repos/dpsanders/TaylorSeries.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiff/TaylorSeries.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/dpsanders/TaylorSeries.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/dpsanders/TaylorSeries.jl/latest)

This package contains a static implementation of the Taylor1 structure in TaylorSeries.jl.

## Key Differences

This package introduces the `STaylor1{N,T}` structure which can be used as a
replacement for the `Taylor1{T}` structure of TaylorSeries.jl.  Arithmetic operators are defined via generated functions and do not allocate.

Key differences are:
- Arithmetic with `STaylor1{N,T}` structures use entirely immutable storage types
for calculations whereas `Taylor1{T}` calculations make use of a number of
arrays for intermediate storage. As a consequence, the `STaylor1{N,T}` implementation
will be significant faster for calculations involving low-order Taylor series.
- The `STaylor1{N,T}` structure stores coefficients as an `NTuple` rather than
an array.
- Constructors: Most constructors are similar to those used for Taylor1. The
constructor `STaylor1(x::T, v::Val{N})` is used in place of `Taylor1(x::T, order::Int)`.
- In place functions (e.g. `tan!(a,c,k)`) are not supported.
- Currently, **tan**, **tanh**, **asin**, **acos**, **atan** are
unsupported.

#### License

`StaticTaylorSeries` is licensed under the [MIT "Expat" license](./LICENSE.md).
