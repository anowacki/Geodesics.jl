# Algorithms

Geodesics.jl can calculate distances and angles using multiple different
algorithms.  By default, Vincenty's (1975) algorithm is used.

Algorithms are optionally passed as the first argument to the calculation
functions such as [`surface_distance`](@ref) or [`Geodesics.forward`](@ref).
They include information on the ellipsoid on which the computation should be
performed.

## Choosing an algorithm

Currently two algorithms are available:

- Vincenty's algorithm ([`Vincenty`](@ref))
- The Haversine algorithm ([`Haversine`](@ref))

Haversine can only be used on spheres (i.e., when flattening is 0).

Algorithms are constructed by passing an [`Ellipsoid`](@ref Ellipsoid{A,T}) to
one of the algorithm constructors, like so:

```@repl example
using Geodesics
alg = Vincenty(Ellipsoid())
```

Note that without any arguments, an algorithm uses the default ellipsoid,
WGS84, so the above is equivalent to just `Vincenty()`:

```@repl example
Vincenty() == Vincenty(Ellipsoid())
```

