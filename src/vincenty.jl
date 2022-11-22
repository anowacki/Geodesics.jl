# Vincenty's method

"""
    Vincenty{A,T} <: Geodesics.GeoDistance

Struct representing the calculation of geodesics using Vincenty's (1975)
algorithm.  For `Float64`-precision input, this algorithm calculates
distances on the Earth to within about 0.1 mm (Vincenty, 1975).

!!! note
    For points which are antipodal or nearly antipodal, Vincenty's algorithm
    can fail to converge, and thus calculations can become inaccurate.  Avoid
    using this method if you expect to encounter such points.

# References

Vincenty, T. (1975). "Direct and Inverse Solutions of Geodesics on the Ellipsoid
with application of nested equations" (PDF). _Survey Review_. XXIII (176): 88–93.
[doi:10.1179/sre.1975.23.176.88](https://doi.org/10.1179/sre.1975.23.176.88)

"""
struct Vincenty{A,T} <: GeoDistance
    ell::Ellipsoid{A,T}
end

"""
    Vincenty(a=$(EARTH_R_MAJOR_WGS84), f=$(EARTH_F_WGS84))
    Vincenty(; a=$(EARTH_R_MAJOR_WGS84), b=nothing, f=$(EARTH_F_WGS84)))

Produce a new ellipsoid from which geodesic calculations can be made using
Vincenty's algorithm.  By default, the WGS84 ellipsoid is used if no
parameters are passed in.

The equatorial radius `a` must be given.  In the first form, the flattening
`f` is given as the second argument.  In the second form, the polar radiud `b`
may be given, in which case `f` is calculated from `a` and `b`.

If no arguments are given, the default (WGS84) ellipsoid is used.

# Examples
```jldoctest
julia> method = Vincenty()
Vincenty{Float64, Float64}(Ellipsoid{Float64, Float64}(6.378137e6, 0.0033528106647474805, 6.356752314245179e6))

julia> Vincenty(1.f0, 0.003f0)
Vincenty{Float32, Float32}(Ellipsoid{Float32, Float32}(1.0f0, 0.003f0, 0.997f0))
```

See also: [`Ellipsoid`](@ref)

---

    Vincenty(ellipsoid::Ellipsoid)

Construct `Vincenty` using an ellipsoid.
"""
Vincenty(; kwargs...) = Vincenty(Ellipsoid(; kwargs...))

Vincenty(a, f) = Vincenty(Ellipsoid(a, f))

function angular_distance(method::Vincenty{A,T}, lon0, lat0, lon1, lat1;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    lon0, lat0, lon1, lat1 = T.((lon0, lat0, lon1, lat1))
    if degrees
        lon0, lat0, lon1, lat1 = deg2rad(lon0), deg2rad(lat0), deg2rad(lon1), deg2rad(lat1)
    end
    result = inverse(method, lon0, lat0, lon1, lat1; check=check)
    degrees ? rad2deg(result.angular_distance) : result.angular_distance
end

function azimuth(method::Vincenty{A,T}, lon0, lat0, lon1, lat1;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    lon0, lat0, lon1, lat1 = T.((lon0, lat0, lon1, lat1))
    if degrees
        lon0, lat0, lon1, lat1 = deg2rad(lon0), deg2rad(lat0), deg2rad(lon1), deg2rad(lat1)
    end
    result = inverse(method, lon0, lat0, lon1, lat1; check=check)
    degrees ? rad2deg(result.azimuth) : result.azimuth
end

function angular_step(method::Vincenty{A,T}, lon, lat, azimuth, distance;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    lon, lat, azimuth, distance = T.((lon, lat, azimuth, distance))
    if degrees
        lon, lat, azimuth, distance = deg2rad(lon), deg2rad(lat), deg2rad(azimuth), deg2rad(distance)
    end
    result = forward_angle(method, lon, lat, azimuth, distance; check=check)
    if degrees
        rad2deg(result.lon), rad2deg(result.lat), rad2deg(result.backazimuth)
    else
        result.lon, result.lat, result.backazimuth
    end
end

function surface_distance(method::Vincenty{A,T}, lon0, lat0, lon1, lat1;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    lon0, lat0, lon1, lat1 = T.((lon0, lat0, lon1, lat1))
    if degrees
        lon0, lat0, lon1, lat1 = deg2rad(lon0), deg2rad(lat0), deg2rad(lon1), deg2rad(lat1)
    end
    result = inverse(method, lon0, lat0, lon1, lat1; check=check)
    result.surface_distance
end

function surface_step(method::Vincenty{A,T}, lon, lat, azimuth, distance;
    degrees::Bool=true, check::Bool=true
) where {A,T}
    # Don't convert distance in case it's unitful; it's done in forward anyway
    lon, lat, azimuth = T.((lon, lat, azimuth))
    if degrees
        lon, lat, azimuth = deg2rad(lon), deg2rad(lat), deg2rad(azimuth)
    end
    result = forward(method, lon, lat, azimuth, distance; check=check)
    lon′ = degrees ? rad2deg(result.lon) : result.lon
    lat′ = degrees ? rad2deg(result.lat) : result.lat
    baz = degrees ? rad2deg(result.backazimuth) : result.backazimuth
    lon′, lat′, baz
end

function forward(method::Vincenty{A,T}, lon, lat, azimuth, distance; check=true)::NamedTuple{(:lon,:lat,:backazimuth),Tuple{T,T,T}} where {A,T}
    check && _check_latitude(lat, false)

    a_units = method.ell.a
    f = method.ell.f
    b_units = method.ell._b

    # Make unitless and get units
    a, a_unit = _remove_units(a_units)
    b, b_unit = _remove_units(b_units)

    # Convert to same units as A and remove units
    distance_nounit, distance_unit = _remove_units(A(distance))

    # Calculations are done with Float64s internally as the tolerances are hard-wired
    lambda1, phi1, alpha12, s = T(lon), T(lat), T(azimuth), T(distance_nounit)
    alpha12 = mod2pi(alpha12)

    TanU1 = (1 - f)*tan(phi1)
    U1 = atan(TanU1)
    sigma1 = atan(TanU1, cos(alpha12))
    Sinalpha = cos(U1)*sin(alpha12)
    cosalpha_sq = 1 - Sinalpha^2

    u2 = cosalpha_sq*(a^2 - b^2)/b^2
    # `A` is already taken as the equatorial radius type parameter
    AA = 1 + (u2/16384)*(4096 + u2*(-768 + u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))

    # Starting with the approximation
    sigma = s/(b*AA)

    # Not moving anywhere. We can return the location that was passed in.
    if sigma == 0
        return (lon=lambda1, lat=phi1, backazimuth=mod2pi(alpha12 + π))
    end

    last_sigma = 2*sigma + 2 # something impossible

    # Iterate the following three equations
    # until there is no significant change in sigma, two_sigma_m and delta_sigma
    local two_sigma_m
    while abs((last_sigma - sigma)/sigma) > #=1.0e-9=# √eps(T)
        two_sigma_m = 2*sigma1 + sigma
        delta_sigma = B*sin(sigma)*(cos(two_sigma_m) + (B/4)*(cos(sigma)*(-1 + 2*cos(two_sigma_m)^2 - (B/6)*cos(two_sigma_m)*(-3 + 4*sin(sigma)^2)*(-3 + 4*cos(two_sigma_m)^2 ))))
        last_sigma = sigma
        sigma = (s/(b*AA)) + delta_sigma
    end

    phi2 = atan((sin(U1)*cos(sigma) + cos(U1)*sin(sigma)*cos(alpha12)),
        ((1-f)*sqrt(Sinalpha^2 + (sin(U1)*sin(sigma) - cos(U1)*cos(sigma)*cos(alpha12))^2)))

    lambda = atan((sin(sigma)*sin(alpha12)),
        (cos(U1)*cos(sigma) - sin(U1)*sin(sigma)*cos(alpha12)))

    C = (f/16)*cosalpha_sq*(4 + f*(4 - 3*cosalpha_sq))

    omega = lambda - (1-C)*f*Sinalpha*(sigma + C*sin(sigma)*(
        cos(two_sigma_m) + C*cos(sigma)*(-1 + 2*cos(two_sigma_m)^2)))

    lambda2 = lambda1 + omega

    alpha21 = atan(Sinalpha, (-sin(U1)*sin(sigma) + cos(U1)*cos(sigma)*cos(alpha12)))
    alpha21 = mod2pi(alpha21 + π)

    return (lon=lambda2, lat=phi2, backazimuth=alpha21)
end

function forward_angle(method::Vincenty{A,T}, lon, lat, azimuth, angular_distance;
    check=true
) where {A,T}
    check && _check_latitude(lat, false)

    f = method.ell.f

    # Notation follows https://en.wikipedia.org/wiki/Vincenty%27s_formulae
    L₁, ϕ₁, α₁, σ = T.((lon, lat, azimuth, angular_distance))

    tanU₁ = (1 - f)*tan(ϕ₁)
    U₁ = atan(tanU₁)

    sinU₁, cosU₁ = sincos(U₁)
    sinα₁, cosα₁ = sincos(α₁)
    sinσ, cosσ = sincos(σ)

    σ₁ = atan(tanU₁, cosα₁)
    sinσ₁, cosσ₁ = sincos(σ₁)

    sinα = cosU₁*sinα₁
    sin²α = sinα^2
    cos²α = 1 - sin²α

    twoσₘ = 2*σ₁ + σ

    ϕ₂ = atan(sinU₁*cosσ + cosU₁*sinσ*cosα₁,
        (1 - f)*sqrt(sin²α + (sinU₁*sinσ - cosU₁*cosσ*cosα₁)^2))
    λ = atan(sinσ*sinα₁, cosU₁*cosσ - sinU₁*sinσ*cosα₁)
    C = f/16*cos²α*(4 + f*(4 - 3*cos²α))
    L = λ - (1 - C)*f*sinα*(σ + C*sinσ*(cos(twoσₘ) + C*cosσ*(-1 + 2*cos(twoσₘ)^2)))
    L₂ = L + L₁
    α₂ = atan(sinα, -sinU₁*sinσ + cosU₁*cosσ*cosα₁)

    return (lon=L₂, lat=ϕ₂, backazimuth=mod2pi(α₂ + π))
end

function inverse(method::Vincenty{A,T}, lon1, lat1, lon2, lat2; check=true) where {A,T}
    if check
        _check_latitude(lat1, false)
        _check_latitude(lat2, false)
    end

    a_units = method.ell.a
    f = method.ell.f
    b_units = method.ell._b

    a, unit = _remove_units(a_units)
    b, _ = _remove_units(b_units)

    lambda1, phi1, lambda2, phi2 = T(lon1), T(lat1), T(lon2), T(lat2)
    tol = sqrt(eps(T))
    if (abs(phi2 - phi1) < tol) && (abs(lambda2 - lambda1) < tol)
        return (surface_distance=zero(T)*unit, angular_distance=zero(T),
            azimuth=zero(T), backazimuth=zero(T))
    end

    TanU1 = (1 - f)*tan(phi1)
    TanU2 = (1 - f)*tan(phi2)

    U1 = atan(TanU1)
    U2 = atan(TanU2)

    sinU1, cosU1 = sincos(U1)
    sinU2, cosU2 = sincos(U2)

    lambda = lambda2 - lambda1
    last_lambda = #=-4000000.0=# T(-4000000) # an impossibe value
    omega = lambda

    # Iterate the following equations until there is no significant change in lambda
    alpha = sigma = sin_sigma = Cos2sigma_m = cos_sigma = sin²σ =
        typemin(T)

    while ((last_lambda < #=-3000000.0=# T(-3000000)) || (lambda != 0)) &&
            (abs((last_lambda - lambda)/lambda) > #=1.0e-9=# √eps(T))
        sin²σ = (cosU2*sin(lambda))^2 +
                         ((cosU1*sinU2 - sinU1*cosU2*cos(lambda)))^2
        sin_sigma = sqrt(sin²σ)
        cos_sigma = sinU1*sinU2 + cosU1*cosU2*cos(lambda)
        sigma = atan(sin_sigma, cos_sigma)

        Sin_alpha = cosU1*cosU2*sin(lambda)/sin(sigma)

        if Sin_alpha >= 1
            Sin_alpha = one(T)
        elseif Sin_alpha <= -1
            Sin_alpha = -one(T)
        end

        alpha = asin(Sin_alpha)
        Cos2sigma_m = cos(sigma) - 2*sinU1*sinU2/cos(alpha)^2
        C = (f/16)*cos(alpha)^2*(4 + f*(4 - 3*cos(alpha)^2))
        last_lambda = lambda
        lambda = omega + (1 - C)*f*sin(alpha)*(sigma +
            C*sin(sigma)*(Cos2sigma_m + C*cos(sigma)*(-1 + 2*Cos2sigma_m^2)))
    end

    u2 = cos(alpha)^2*(a^2 - b^2)/b^2
    # The name `A` is taken by the major radius's type parameter
    AA = 1 + (u2/16384)*(4096 + u2*(-768 + u2*(320 - 175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
    delta_sigma = B*sin_sigma*(Cos2sigma_m + (B/4)*(
        cos_sigma*(-1 + 2*Cos2sigma_m^2) -
        (B/6)*Cos2sigma_m*(-3 + 4*sin²σ)*(-3 + 4*Cos2sigma_m^2)))
    s = b*AA*(sigma - delta_sigma)

    alpha12 = atan((cosU2*sin(lambda)), ( cosU1*sinU2 - sinU1*cosU2*cos(lambda)))
    alpha21 = atan((cosU1*sin(lambda)), (-sinU1*cosU2 + cosU1*sinU2*cos(lambda)))

    alpha12 = mod2pi(alpha12)
    alpha21 = mod2pi(alpha21 + π)

    return (surface_distance=s*unit, angular_distance=sigma,
        azimuth=alpha12, backazimuth=alpha21)
end
