# Type index

This page lists the docstrings for each public type exported
by Geodesics.jl

## Ellipsoids

```@docs
Geodesics.Ellipsoid
Geodesics.Ellipsoid()
Geodesics.Sphere()
```

## Computation algorithms

```@docs
Geodesics.Vincenty
Geodesics.Vincenty()
Geodesics.Haversine
Geodesics.Haversine()
```

Geodesics.jl also defines an abstract supertype of all algorithms
`Geodesics.GeoDistance`.  Developers wishing to add new algorithms
to Geodesics.jl or their own packages should subtype this and
implement `forward` and `inverse` for their type.  These will then
be called by the convenience functions.