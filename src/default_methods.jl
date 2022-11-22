# Methods defining the default algorithm to use for each and docstrings

const DEFAULT_METHOD = Vincenty()

"""
    angular_distance([method::Vincenty,] lon0, lat0, lon1, lat1, degrees=true; check=true) -> Δ

Return the angular distance between points (`lon0`,`lat0`) and (`lon1`,`lat1`)
on a flattened sphere using `method`.  By default,
input and output are in degrees, but specify `degrees` as `false` to use radians.

If `check` is `false`, then no checks are made that the inputs are correct.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.
"""
angular_distance(lon0, lat0, lon1, lat1; degrees::Bool=true, check::Bool=true) =
    angular_distance(DEFAULT_METHOD, lon0, lat0, lon1, lat1; degrees=degrees, check=check)

"""
    angular_step([method::Vincenty,] lon, lat, azimuth, distance, degrees=true; check=true) -> lon′, lat′, backazimuth

Return the longitude `lon′`, latitude `lat′` and backazimuth `baz` reached by
travelling an angular `distance` along `azimuth` from the starting point at
(`lon`,`lat`), using `method`.  By default,
input and output are in degrees, but specify `degrees` as `false` to use radians.

If the ellipsoid used has unitful dimensions, then `distance` must be dimensionally
compatible if it also has units.  If it is unitless, then the units of the input
are assumed to be the same as that of the ellipsoid.

If `check` is `false`, then no checks are made that the inputs are correct.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.
"""
angular_step(lon, lat, azimuth, distance; degrees::Bool=true, check::Bool=true) =
    angular_step(DEFAULT_METHOD, lon, lat, azimuth, distance; degrees=degrees, check=check)

"""
    azimuth([method::Vincenty,] lon0, lat0, lon1, lat1, degrees=true; check=true)

Return the azimuth from point (`lon0`,`lat0`) to point (`lon1`,`lat1`) on a
flattened sphere using `method`.  By default,
input and output are in degrees, but specify `degrees` as `false` to use radians.

If `check` is `false`, then no checks are made that the inputs are correct.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.
"""
azimuth(lon0, lat0, lon1, lat1; degrees::Bool=true, check::Bool=true) = 
    azimuth(DEFAULT_METHOD, lon0, lat0, lon1, lat1; degrees=degrees, check=check)

"""
    surface_distance([method::Vincenty,] lon0, lat0, lon1, lat1, degrees::Bool=true; check=true)

Return the physical distance between points (`lon0`,`lat0`) and (`lon1`,`lat1`) on
the flattened sphere.  Distance is given
in the same units as the equatorial radius.  By default, input angles are in degrees, but specify `degrees`
as `false` to use radians.

If `check` is `false`, then no checks are made that the inputs are correct.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.
"""
surface_distance(lon0, lat0, lon1, lat1; degrees::Bool=true, check::Bool=true) =
    surface_distance(DEFAULT_METHOD, lon0, lat0, lon1, lat1; degrees=degrees, check=check)

"""
    surface_step([method::Vincenty,] lon, lat, azimuth, distance, degrees=true; check=true) -> lon′, lat′, backazimuth

Return the longitude `lon′`, latitude `lat′` and backazimuth `baz` reached by
travelling a `distance` along the surface of an ellipsoid, taking an azimuth
`azimuth` from the starting point at (`lon`,`lat`), using `method`.  By default,
input and output are in degrees, but specify `degrees` as `false` to use radians.

If the ellipsoid used has unitful dimensions, then `distance` must be dimensionally
compatible if it also has units.  If it is unitless, then the units of the input
are assumed to be the same as that of the ellipsoid.

If `check` is `false`, then no checks are made that the inputs are correct.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.
"""
surface_step(lon, lat, azimuth, distance; degrees::Bool=true, check::Bool=true) =
    surface_step(DEFAULT_METHOD, lon, lat, azimuth, distance; degrees=degrees, check=check)

"""
    forward([method=Vincenty(),] lon, lat, azimuth, distance; check=true) -> lon′, lat′, backazimuth

Return the longitude `lon′` and latitude `lat′` and `backazimuth` of a projected
point, reached by travelling along an `azimuth` for `distance` from an original
point at (`lon`, `lat`).

Coordinates and azimuth are in radians.

If the ellipsoid used has unitful dimensions, then `distance` must be dimensionally
compatible if it also has units.  If it is unitless, then the units of the input
are assumed to be the same as that of the ellipsoid.

If `check` is `false`, then no checks are made that the inputs are correct.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.

#### References

1. Vincenty, T. (1975). "Direct and Inverse Solutions of Geodesics on the Ellipsoid
   with application of nested equations" (PDF). Survey Review. XXIII (176): 88–93.
   doi:10.1179/sre.1975.23.176.88
"""
forward(lon, lat, azimuth, distance; check=true) =
    forward(DEFAULT_METHOD, lon, lat, azimuth, distance; check=check)

"""
    forward_angle([method=Vincenty(),] lon, lat, azimuth, angular_distance) -> lon, lat, backazimuth

Return the longitude `lon′` and latitude `lat′` and `backazimuth` of a projected
point, reached by travelling along an `azimuth` for a spherical angle of
`angular_distance` from an original point at (`lon`, `lat`).

Coordinates, azimuth and spherical distance are in radians.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.

If `check` is `false`, then no checks are made that the inputs are correct.
"""
forward_angle(lon, lat, azimuth, angular_distance; check=true) =
    forward(DEFAULT_METHOD, lon, lat, azimuth, angular_distance; check=check)

"""
    inverse([method=Vincenty(),] lon1, lat1, lon2, lat2; check=true) -> surface_distance, angular_distance, azimuth, backazimuth

Return the `surface_distance`, `angular_distance`, `azimuth` and `backazimuth`
between two points with longitudes `lon1` and `lon2`, and latitudes `lat1`
and `lat2`, using `method`.

Coordinates and angles are in radians, whilst `distance` is in the same units as
the equatorial radius in `method`.

If no `method` is supplied, then computations are performed using Vincenty's (1975)
method on the WGS84 ellipsoid.

If `check` is `false`, then no checks are made that the inputs are correct.
"""
inverse(lon1, lat1, lon2, lat2; check=true) =
    inverse(DEFAULT_METHOD, lon1, lat1, lon2, lat2; check=check)
