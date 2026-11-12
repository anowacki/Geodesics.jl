# Function index

This page lists the docstrings for each public function.  You can
also see this information in the Julia REPL by typing `?` then
the name of the function.

## High-level functions
These functions make use of the low-level functions to allow you to
calculate specific quantities easily.

```@docs
Geodesics.angular_distance
Geodesics.angular_step
Geodesics.azimuth
Geodesics.surface_distance
Geodesics.surface_step
```

## Low-level functions
These functions solve the geodesic inverse or forward problems,
returning all available information.

```@docs
Geodesics.forward
Geodesics.inverse
```
