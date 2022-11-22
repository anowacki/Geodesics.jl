"""
# `Geodesics`

Compute geodesics, the shortest path between two points, on a flattened
sphere.


## User functions
### Exported functions
The following functions are part of the public interface to the module:

- [`angular_distance`](@ref): Compute the distance between two points
  in terms of spherical angle.
- [`angular_step`](@ref): Compute the final point reached by travelling
  a spherical angle distance along a known azimuth from a starting point,
  returning also the backazimuth from the final point to the first
- [`azimuth`](@ref): Compute the azimuth between two points
- [`surface_distance`](@ref): Compute the distance between two points in
  terms of physical surface distance.
- [`surface_step`](@ref): Compute the final point reached by travelling
  a physical surface distance along a known azimuth from a starting point,
  returning also the backazimuth from the final point to the first

### Unexported functions
Two functions ([`Geodesics.forward`](@ref) and [`Geodesics.inverse`](@ref))
are not exported, but can be used by users.

`Geodesics.forward` solves the 'forward problem', meaning it finds
the final point reached when travelling from a starting point along a known
azimuth for a set surface distance.  It returns the final coordinates,
backazimuth and spherical angle distance.

`Geodesics.inverse` solves the 'inverse problem', meaning it finds
the azimuth, backazimuth, spherical angle and surface distance between
two points on the ellipsoid.

### Default ellipsoid and algorithm
By default, all of the above solve the problems using Vincenty's algorithm
on the WGS84 ellipsoid.


## Types
### Ellipsoids
The parameters of an ellipsoid of rotation are the major radius ``a`` and
flattening ``f``.  You can construct an ellipsoid with [`Ellipsoid`](@ref)
and pass it to any of the following algorithms.

### Algorithms
You can control how azimuths, distances, etc., are calculated by passing
different algorithms to methods exported by this module.  The following
algorithms are available:

- [`Haversine`](@ref): Use the Haversine formulae.  Only possible for
  'ellipsoids' which are spherical (i.e., where ``f = 0``.
- [`Vincenty`](@ref): Use Vincenty's algorithm.  Accurate with `Float64`s to
  about 0.1 mm on the (WGS84) Earth, but does not converge for antipodal points.
  Default algorithm.

Algorithms can either be passed an `Ellipsoid` or raw ellipsoid parameters
upon construction.
"""
module Geodesics

export
    # Types
    Ellipsoid,
    Haversine,
    Sphere,
    Vincenty,

    # Methods
    angular_distance,
    angular_step,
    azimuth,
    surface_distance,
    surface_step

# Utility functions
include("util.jl")

# Constants
include("constants.jl")

# Generic and abstract types
include("types.jl")

# Algorithms
include("vincenty.jl")
include("haversine.jl")

# Default methods
include("default_methods.jl")

end # module
