# Haversine formula

"""
    Haversine <: Geodesics.GeoDistance

Struct representing the calculation of geodesics using the Haversine
formula.

!!! note
    The Haversine formula only applies to spheres; that is, it only
    works for ellipsoids with no flattening (``f = 0``).  If you try
    to use `Haversine` with a non-spherical ellipsoid, an error
    will be thrown.

# References
1. https://en.wikipedia.org/wiki/Haversine_formula
"""
struct Haversine{A,T} <: GeoDistance
    ell::Ellipsoid{A,T}

    function Haversine{A,T}(ell) where {A,T}
        iszero(ell.f) ||
            throw(ArgumentError(
                """non-zero flattening not allowed for Haversine, since computations
                   require a sphere with the haversine method"""))
        new{A,T}(ell)
    end

    Haversine(ell::Ellipsoid{A,T}) where {A,T} = Haversine{A,T}(ell)
end

"""
    Haversine(a=$(EARTH_R_MEAN_WGS84))
    Haversine(; a=$(EARTH_R_MEAN_WGS84))

Produce a new spherical ellipsoid from which geodesic calculations can be
made using the Haversine formula.  By default, a sphere with the average
WGS84 radius is used if no parameters are passed in.

# Examples
```
julia> method = Haversine()
Haversine{Float64, Float64}(Ellipsoid{Float64, Float64}(6.371008771415059e6, 0.0, 6.371008771415059e6))

julia> Haversine(1.0)
Haversine{Float64, Float64}(Ellipsoid{Float64, Float64}(1.0, 0.0, 1.0))
```

See also: [`Ellipsoid`](@ref)

---

    Haversine(ellipsoid::Ellipsoid)

Construct `Haversine` using an ellipsoid.
"""
Haversine(; kwargs...) =
    Haversine(Ellipsoid(; a=EARTH_R_MEAN_WGS84, f=zero(EARTH_R_MEAN_WGS84), kwargs...))
Haversine(a) = Haversine(Ellipsoid(a, zero(a)/oneunit(a)))

function angular_distance(::Haversine{A,T}, lon1, lat1, lon2, lat2;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    if check
        _check_latitude(lat1, degrees)
        _check_latitude(lat2, degrees)
    end

    lon1, lat1, lon2, lat2 = T.((lon1, lat1, lon2, lat2))
    if degrees
        lon1, lat1, lon2, lat2 = deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2)
    end
    Δlon = lon2 - lon1
    sinΔlon, cosΔlon = sincos(Δlon)
    sin_lat1, cos_lat1 = sincos(lat1)
    sin_lat2, cos_lat2 = sincos(lat2)
    d = atan(sqrt(
               (cos_lat2*sinΔlon)^2 + (cos_lat1*sin_lat2 -
                sin_lat1*cos_lat2*cosΔlon)^2),
               sin_lat1*sin_lat2 + cos_lat1*cos_lat2*cosΔlon
              )
    degrees ? rad2deg(d) : d
end

function angular_step(method::Haversine{A,T}, lon, lat, az, delta;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    check && _check_latitude(lat, degrees)

    lon, lat, az, delta = T.((lon, lat, az, delta))

    if iszero(delta)
        return lon, lat, (degrees ? mod(az + 180, 360) : mod2pi(az + π))
    end

    if degrees
        lon, lat, az, delta = deg2rad(lon), deg2rad(lat), deg2rad(az), deg2rad(delta)
    end

    sin_lat, cos_lat = sincos(lat)
    sinΔ, cosΔ = sincos(delta)
    sin_az, cos_az = sincos(az)

    lat2 = asin(sin_lat*cosΔ + cos_lat*sinΔ*cos_az)
    lon2 = lon + atan(sin_az*sinΔ*cos_lat,
                      cosΔ - sin_lat*sin(lat2))
    # This includes a check for latitude range
    baz = azimuth(method, lon2, lat2, lon, lat; degrees=false, check=check)

    if degrees
        lon2, lat2, baz = rad2deg(lon2), rad2deg(lat2), rad2deg(baz)
    end

    lon2, lat2, baz
end

function azimuth(::Haversine{A,T}, lon1, lat1, lon2, lat2;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    if check
        _check_latitude(lat1, degrees)
        _check_latitude(lat2, degrees)
    end

    lon1, lat1, lon2, lat2 = T.((lon1, lat1, lon2, lat2))
    if degrees
        lon1, lat1, lon2, lat2 = deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2)
    end
    Δlon = lon2 - lon1
    sinΔlon, cosΔlon = sincos(Δlon)
    sin_lat1, cos_lat1 = sincos(lat1)
    sin_lat2, cos_lat2 = sincos(lat2)
    azimuth = atan(sinΔlon*cos_lat2,
                   cos_lat1*sin_lat2 - sin_lat1*cos_lat2*cosΔlon)
    azimuth = mod2pi(azimuth)
    degrees ? rad2deg(azimuth) : azimuth
end

function surface_distance(method::Haversine{A,T}, lon1, lat1, lon2, lat2;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    Δ = angular_distance(method, lon1, lat1, lon2, lat2; degrees=degrees, check=check)
    method.ell.a*(degrees ? deg2rad(Δ) : Δ)
end

function surface_step(method::Haversine{A,T}, lon, lat, azimuth, distance;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    lon, lat, azimuth = T.((lon, lat, azimuth))

    if iszero(distance)
        return lon, lat
    end

    if degrees
        lon, lat, azimuth = deg2rad.((lon, lat, azimuth))
    end

    delta = distance/method.ell.a

    lon2, lat2, baz = angular_step(method, lon, lat, azimuth, delta;
        degrees=false, check=check)
    if degrees
        lon2, lat2, baz = rad2deg.((lon2, lat2, baz))
    end

    lon2, lat2, baz
end

function forward(method::Haversine{A,T}, lon, lat, azimuth, distance; check=true) where {A,T}
    check && _check_latitude(lat, false)
    lon, lat, azimuth = T.((lon, lat, azimuth))
    if iszero(distance)
        (lon=lon, lat=lat, backazimuth=mod2pi(azimuth + π))
    else
        lon′, lat′ = surface_step(method, lon, lat, azimuth, distance;
            degrees=false, check=false)
        baz = Geodesics.azimuth(method, lon′, lat′, lon, lat; degrees=false, check=false)

        (lon=lon′, lat=lat′, backazimuth=baz)
    end
end

function forward_angle(method::Haversine{A,T}, lon, lat, azimuth, angular_distance;
    check=true
) where {A,T}
    check && _check_latitude(lat, false)
    lon′, lat′, baz = angular_step(method, lon, lat, azimuth, angular_distance;
        degrees=false, check=false)
    (lon=lon′, lat=lat′, backazimuth=baz)
end

function inverse(method::Haversine{A,T}, lon1, lat1, lon2, lat2; check=true) where {A,T}
    if check
        _check_latitude(lat1, false)
        _check_latitude(lat2, false)
    end
    azi = azimuth(method, lon1, lat1, lon2, lat2; degrees=false, check=false)
    baz = azimuth(method, lon2, lat2, lon1, lat1; degrees=false, check=false)
    delta = angular_distance(method, lon1, lat1, lon2, lat2; degrees=false, check=false)
    distance = method.ell.a*delta

    (surface_distance=distance, angular_distance=delta, azimuth=azi, backazimuth=baz)
end
