# Numerical Tests in Computing invariant manifolds in non-smooth dynamical systems
 
To run the codes in this repository, first install [julia](https://julialang.org/) and then install the following packages:
- InvariantManifolds
- LinearAlgebra
- StaticArrays
- OrdinaryDiffEq
- GLMakie
- DataInterpolations

You can install them using the following command:
```julia
using Pkg
Pkg.add(url="https://github.com/Xiaomingzzhang/InvariantManifolds.jl")
Pkg.add("StaticArrays")
Pkg.add("OrdinaryDiffEq")
Pkg.add("GLMakie")
Pkg.add("DataInterpolations")
```

You can run the codes in the following way:
```julia
include("src/lorenzimpact.jl")
```