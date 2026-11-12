# Unitful ellipsoids

Geodesics.jl is written to be as agnostic as possible to the types
used to construct `Ellipsoid`s, and to accept any types for the
arguments to the calculation functions.  Consequently, Geodesics.jl
handles ellipsoids with units.  In particular, Geodesics.jl tests
that units from the [Unitful.jl](https://github.com/PainterQubits/Unitful.jl)
package work.

## Ellipsoids with unitful axes

To construct an ellipsoid with units for the axes, simply pass a value
with a unit as the equatorial semi-axis `a` to the `Ellipsoid` constructor
or the constructor for one of the algorithms:

```@repl example
using Geodesics, Unitful
Ellipsoid(a=1_000_000*u"m", f=0.001)
```

If defining the ellipsoid with the semi-major and -minor axes, the units of
the equatorial 'radius' are used and the polar 'radius' values are converted
to the same units as those of ``a``:

```@repl  example
Ellipsoid(a=100*u"m", b=0.1u"km")
```

### Returned units
By default, the units returned for distances are the same as that initially
passed in to the ellipsoid.

For example, when calculating distances along the surface of the ellipsoid:
```@repl example
alg = Vincenty(a=6378.0*u"km", f=1/298)
surface_distance(alg, 0, 0, 10, 10)
```

### Input units
In functions which take distances as an input (e.g., [`surface_step`](@ref)),
a unitless value is assumed to have the same units as the ellipsoid:

```@repl example
surface_step(alg, 0, 0, 45, 10)
```

However, you can use units as input distances quite happily:

```@repl example
surface_step(alg, 0, 0, 45, 10*u"km")
surface_step(alg, 0, 0, 45, 10_000*u"m")
```
