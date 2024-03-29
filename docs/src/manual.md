# User manual

## Importing the module

Once you have [installed Geodesics.jl](#Install), you can load the module
by doing

```@repl example
using Geodesics
```

## Basic functions

### [`Geodesics.forward`](@ref)

A standard forward Vincenty computation finds the end point on a flattened sphere,
given a starting location at longitude-latitude `(lon,lat)`, a forward azimuth `az`
and a distance `dist` travelled along the surface.  You also need to specify the
semimajor (equatorial) radius of the ellipsoid `a` and the flattening `f`.

Let's travel 30 km northeast of
[Trafalgar Square](https://www.google.co.uk/maps/place/Trafalgar+Square/@51.50809,-0.1302376,17z/),
london, and see where we end up:

```@repl example
using Geodesics
lon, lat, az, dist, a, f = deg2rad(0.1281), deg2rad(51.5080), deg2rad(45), 30_000, Geodesics.EARTH_R_MAJOR_WGS84, Geodesics.F_WGS84
Geodesics.forward(lon, lat, az, dist, a, f)
lon1, lat1, backazimuth = rad2deg.(ans)
```

So we get to [here](https://www.google.com/maps/place/51.69826376194629,0.4349208715005321),
which is most of the way to Chelmsford in Essex.  Exciting!

Note that here, we've use the semimajor Earth radius from the WGS84 ellipsoid,
and WGS84's flattening, avaialable in the geodesics package.  If we'd assumed a
perfectly spherical Earth, then we end up very slightly
[elsewhere](https://www.google.com/maps/place/51.698479560055226,0.4360730564880621).

### [`Geodesics.inverse`](@ref)

A standard inverse Vincenty calculation finds the distance, azimuth and backazimuth
between two known points on the ellipsoid, (`lon1,lat1`) and (`lon2,lat2`).  Again,
we need to specify the semimajor radius and flattening.

Let's see how far it is (as the crow flies) between
[Jarrow](https://en.wikipedia.org/wiki/Jarrow) and
[Parliament](https://en.wikipedia.org/wiki/Palace_of_Westminster):

```@repl example
lon1, lat1 = -1.4951547, 54.967618; # Jarrow
lon2, lat2 = -0.1270032, 51.4994794; # Palace of Westminster
a, f = Geodesics.EARTH_R_MAJOR_WGS84, Geodesics.F_WGS84;
dist, az, baz = Geodesics.inverse(deg2rad.((lon1, lat1, lon2, lat2))..., a, f)
```

So the [Jarrow Marchers](https://en.wikipedia.org/wiki/Jarrow_March) had at the very least
(ignoring topography) 396 km to walk!


## Convenience functions

`Geodesics.jl` comes with convenience functions which use angles in either degrees
or radians, as you prefer.  It is also often true that you are only interested in
one quantity at the time, such only wanting to find out the forward azimuth between
two known points.  In this case, use `Geodesics.azimuth` to find which azimuth you should
sail along (on a perfectly calm day) to reach St Malo from Jersey:

```@repl example
lon0, lat0 = -2.117641, 49.176924; # St. Helier, Jersey
lon1, lat1 = -2.032614, 48.641570; # St. Malo, France
Geodesics.azimuth(lon0, lat0, lon1, lat1) # Uses degrees and WGS84 ellipsoid by default
```

Note that I didn't need to specify the flattening (which defaults to that of
the WGS84 ellipsoid), but this can be overridden using the `f` keyword argument.

The full list of convenience functions:

- [`angular_distance`](@ref Geodesics.angular_distance): Find the surface distance between two points in terms of
  an angle measured from the centre of the ellipsoid.
- [`surface_distance`](@ref Geodesics.surface_distance): The distance between two points.
- [`angular_step`](@ref Geodesics.angular_step): Find the end point and backazimuth from one point when travelling
  a set angular distance along a defined azimuth.
- [`azimuth`](@ref Geodesics.azimuth): The forward azimuth between two points.


## Choice of geodesic calculation

This package so far only implements Vincenty's methods, but others are available.
Pull requests to add these are welcome.

## References

Thaddeus Vincenty published the forward and inverse methods used in this package in
the following paper:

1. Vincenty, T. (1975). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
   application of nested equations" (PDF). Survey Review. XXIII (176): 88–93. [doi:10.1179/sre.1975.23.176.88](https://doi.org/10.1179/sre.1975.23.176.88).
