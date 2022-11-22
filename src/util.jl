# Utility functions

"""
    _flattening(a, b) -> f

Compute the flattening `f` from the equatorial radius `a`
and the polar radius `b` of an ellipsoid of rotation whose
symmetry axis passes through the pole.
"""
_flattening(a, b) = (a - b)/a

"""
    _polar_radius(a, f) -> b

Return the polar radius `b` for an ellipsoid of rotation given
the equatorial radius `a` and flattening `f`.
"""
_polar_radius(a, f) = a*(1 - f)

"""
    _remove_units(x) -> (x=x̃, unit=unit)

For a quantity `x` which may or may not have units, return the value
of `x` stripped of its units `x̃`, and the `unit`.
"""
_remove_units(x) = (x=x/oneunit(x), unit=oneunit(x))

"""
    _promote_radius_and_flattening(a, f)

Return a pair of `a` and `f` which retain units if any on `a`, but
promoted to a common floating-point type.
"""
_promote_radius_and_flattening(a, f) = promote(float(a), float(f))

"""
    _check_latitude(lat, degrees)

Throw an `ArgumentError` if `lat` is outside the range [-π/2, π/2]
when `degrees` is `false`, or [-90°, 90°] when `degrees` is `true`.
"""
function _check_latitude(lat, degrees)
    if degrees
        abs(lat) <= 90 || throw(ArgumentError("latitude not in range [-90°, 90]"))
    else
        abs(lat) <= π/2 || throw(ArgumentError("latitude not in range [-π/2, π/2]"))
    end
    nothing
end
