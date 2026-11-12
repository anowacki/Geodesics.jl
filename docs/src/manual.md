# User manual

## Importing the module

Once you have [installed Geodesics.jl](@ref Install), you can load the module
by doing

```@repl example
using Geodesics
```

## High-level functions

`Geodesics.jl` comes with high-level functions which use angles in either degrees
or radians, as you prefer.  It is also often true that you are only interested in
one quantity at the time, such only wanting to find out the forward azimuth between
two known points.  In this case, use `azimuth` to find which azimuth you should
sail along (on a perfectly calm day) to reach St Malo from Jersey:

```@repl example
lon0, lat0 = -2.117641, 49.176924; # St. Helier, Jersey
lon1, lat1 = -2.032614, 48.641570; # St. Malo, France
azimuth(lon0, lat0, lon1, lat1) # Uses degrees and WGS84 ellipsoid by default
```

Note that I didn't need to specify the flattening (which defaults to that of
the WGS84 ellipsoid), but this can be overridden by explicitly passing a
method as the first argument.

The full list of convenience functions:

- [`angular_distance`](@ref Geodesics.angular_distance): Find the surface distance between two points in terms of
  an angle measured from the centre of the ellipsoid.
- [`surface_distance`](@ref Geodesics.surface_distance): The distance between two points on the surface of the body.
- [`angular_step`](@ref Geodesics.angular_step): Find the end point and backazimuth from one point when travelling
  a set angular distance along a defined azimuth.
- [`surface_step`](@ref Geodesics.surface_step): Find the end point and backazimuth
  from one point when travelling a set distance along a defined azimuth on the
  surface of the ellipsoid.
- [`azimuth`](@ref Geodesics.azimuth): The forward azimuth between two points.

All high-level functions take as an optional first argument the
[algorithm](#Algorithms) you want to use.  This defaults to `Vincenty()`,
which in turn defaults to using the WGS84 ellipsoid.

### `angular_distance` and `surface_distance`
A common problem on a sphere is to find out how far it is between two known
points.  Let's see how far it is (as the crow flies) between
[Jarrow](https://en.wikipedia.org/wiki/Jarrow) and
[Parliament](https://en.wikipedia.org/wiki/Palace_of_Westminster):

```@repl example
lon1, lat1 = -1.4951547, 54.967618; # Jarrow
lon2, lat2 = -0.1270032, 51.4994794; # Palace of Westminster
surface_distance(lon1, lat1, lon2, lat2)
```

So the [Jarrow Marchers](https://en.wikipedia.org/wiki/Jarrow_March) had at the very least
(ignoring topography) 396 km to walk!

Expressing this in terms of angle at the centre of the Earth:

```@repl example
angular_distance(lon1, lat1, lon2, lat2)
```

### `surface_step` and `angular_step`
Another common problem on a sphere is to calculate where one ends up after
travelling some set distance along a starting azimuth.

Let's travel 30 km northeast of
[Trafalgar Square](https://www.google.co.uk/maps/place/Trafalgar+Square/@51.50809,-0.1302376,17z/),
london, and see where we end up:

```@repl example
lon, lat, az, dist = 0.1281, 51.5080, 45, 30_000;
lon_end, lat_end, baz = surface_step(lon, lat, az, dist)
```

So we get to [here](https://www.google.com/maps/place/51.698263761946066,0.43492087150187003),
which is most of the way to Chelmsford in Essex.  Exciting!

You can see that `surface_step` and `angular_step` return the longitude
and latitude of the final point, and the azimuth from the final point
to the first (known as the backazimuth).

Above we've used the default WGS84 ellipsoid and the `Vincenty`
algorithm.  If we'd assumed a
perfectly spherical Earth, then we end up very slightly
[elsewhere](https://www.google.com/maps/place/51.698373226487455,0.43590004438265734):

```@repl example
surface_step(Haversine(), lon, lat, az, dist)
```


## Low-level functions

As well as the high-level functions above, you can also use the low-level
functions [`Geodesics.forward`](@ref) and [`Geodesics.inverse`](@ref).  These are
not exported and their input and output is in radians for angular quantities.

### [`Geodesics.forward`](@ref)

A standard forward geodesic computation finds the end point on a flattened sphere,
given a starting location at a longitude-latitude point `(lon,lat)`, a forward
azimuth `az` and a distance `dist` travelled along the surface.  It als computes
the backazimuth from the end point back to the first.  High-level functions
`angular_step` and `surface_step` rely on this calculation.

Our trafalgar square example above now has to be written like this:

```@repl example
lon, lat, az, dist = deg2rad(0.1281), deg2rad(51.5080), deg2rad(45), 30_000;
alg = Vincenty()
out = Geodesics.forward(alg, lon, lat, az, dist)
lon1, lat1, backazimuth = rad2deg.((out.lon, out.lat, out.backazimuth))
```

As before, we can change our algorithm by passing a different one as the first
argument:
```@repl example
Geodesics.forward(Haversine(), lon, lat, az, dist)
```

### [`Geodesics.inverse`](@ref)

A standard inverse calculation finds the distance, azimuth and backazimuth
between two known points on the ellipsoid, (`lon1,lat1`) and (`lon2,lat2`).
Again, this gives us the building block to enable the high-level functions
`azimuth`, `surface_distance` and `angular_distance`.

Taking the same example as before:

```@repl example
lon1, lat1 = -1.4951547, 54.967618; # Jarrow
lon2, lat2 = -0.1270032, 51.4994794; # Palace of Westminster
Geodesics.inverse(deg2rad.((lon1, lat1, lon2, lat2))...)
```

Note that the default argument to `Geodesics.forward` and `Geodesics.inverse`
is `Vincenty()`, so we didn't need to pass it in the examples above.
