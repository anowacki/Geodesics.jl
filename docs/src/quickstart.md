# Quickstart

Geodesics.jl allows you to calculate:

- the azimuth from one point on a flattened sphere to another;
- the distance along the surface between two points, either in units of
  distance or the angle subtended at the centre of the body; and
- the ending point reached by travelling from a starting point along a given
  direction for a given distance.

By default it does this using Vincenty's algorithm on the WGS84 ellipsoid,
giving distances accurate to mm.

## Preliminaries

Bring the functions needed into scope before running any of the other examples.

```@repl quickstart
using Geodesics
```

## Azimuth between two points

```@repl quickstart
lon0, lat0 = -2.117641, 49.176924; # St. Helier, Jersey
lon1, lat1 = -2.032614, 48.641570; # St. Malo, France
azimuth(lon0, lat0, lon1, lat1) # Uses degrees and WGS84 ellipsoid by default
```

## Distance between two points

```@repl quickstart
surface_distance(lon0, lat0, lon1, lat1) # Returns distance in m
angular_distance(lon0, lat0, lon1, lat1) # Angle at centre of Earth between points
```

## End point
Travelling a certain distance along the surface from a starting point:
```@repl quickstart
distance = 100_000; # Surface distance in m
azi = 45; # Azimuth from starting point in degrees
surface_step(lon0, lat0, azi, distance)
```

Travelling a certain angular distance about the centre of the Earth from a
starting point:
```@repl quickstart
ang_dist = 1; # degrees
angular_step(lon0, lat0, azi, ang_dist)
```
