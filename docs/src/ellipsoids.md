# Ellipsoids

By default, all the functions provided as part of Geodesics.jl
use the WGS84 ellipsoid.  However, users (i.e., you) can
create new ellipsoids for any spherical or uniaxial elliptical
body using the `Ellipsoid` class.

## `Ellipsoid{A,T}`

An ellipsoid describes the shape of a body.  In Geodesics.jl, such
a body can be an 'ellipsoid of rotation', defined by two of three
related parameters:

- the equatorial radius or semimajor axis ``a``;
- the polar radius or semiminor axis ``b``; and
- the flattening ``f``.

Flattening is defined as

``
    f = \frac{a - b}{a} \, .
``

(Here we often refer to the 'radius' of an ellipsoid, but this is not
technically correct, as only spheres have a radius.  But it is often
convenient to talk, for example for an oblate ellipsoid, of the
'polar radius' rather than the 'semi-minor axis length' for brevity.

An instance of the [`Ellipsoid`](@ref) type can be created in two ways:

1. by passing in the equatorial 'radius' ``a`` and flattening ``f``, or
2. passing in the equatorial 'radius' ``a`` and polar 'radius' ``b``.

!!! note
    `Ellipsoid` only supports uniaxial ellipsoids; that is ellipsoids
    where the equatorial radius is the same everywhere, but the polar
    radius may be larger or smaller than the equatorial radius.
    *(Triaxial ellipsoids or more complicated surfaces could be added in
    future if someone wanted to contribute that.
    [Open an issue](https://github.com/anowacki/Geodesics.jl/issues/new)
    if you are interested in helping.)*

As an example, we can use the semimajor axis and flattening for the
World Geodetic System ellipsoid like so:

```@repl example
using Geodesics
a = 6_378_137
f = 0.0033528106647474805
wgs84 = Ellipsoid(a, f)
```

It so happens that the default parameters for an `Ellipsoid` are indeed
the WGS84 parameters, so:

```@repl example
wgs84 == Ellipsoid()
```

We can also pass the semimajor and semiminor axes.
For example, to create an ellipsoid for Mars using the parameters of
[Ardalan et al., (2010)](https://doi.org/10.1007/s11038-009-9342-7):

```@repl example
mars = Ellipsoid(a=3_395_428, b=3_377_678)
```

## Spheres

Sometimes a spherical approximation of a body is good enough, or you are
only interested in truly spherical distances.  In this case, you can
construct an ellipsoid with no flattening (i.e., where ``a = b`` and ``f = 0``)
either by passing `f=0` to the `Ellipsoid` constructor, or you can use
the convenience constructor [`Sphere`](@ref):

```@repl example
sphere = Sphere(1) # The unit sphere
```

## Type parameters

`Ellipsoid`s have two type parameters: `A` and `T`:
- `A` is the type of the semimajor axis.
- `T` is the type of the flattening.

Usually (e.g., for the default ellipsoid), `A` and `T` are the same
and just define the precision of the type used.  (Note that the precision
of `A` and `T` are enforced to be the same).  However, `A` can have the
same precision but be of a different type since it can have units.
See [Unitful ellipsoids](@ref) for more information.

Likewise, `A` or `T` (or both) could be uncertain.
[Measurements](https://github.com/JuliaPhysics/Measurements.jl)
package is one Julia package which provides numbers which have associated
uncertainties, which can be propagated through calculations as expected:

In this example we have a sphere with a radius known to 10%:
```@repl example
using Measurements
ell = Sphere(1 ± 0.1) # Sphere with radius (1 ± 0.1) m
alg = Haversine(ell)
surface_distance(alg, 0, 0, 10, 10) # Distance from (lon, lat = 0, 0) to (10, 10)
```
