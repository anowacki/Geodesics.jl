# Types used by the rest of the module, and abstract types

"""
    GeoDistance

Abstract type from which new types to calculate geodesics should subtype.
"""
abstract type GeoDistance end

"""
    Ellipsoid{A,T}

Struct representing an ellipsoid of rotation; that is, a flattened
sphere.

The type parameter `A` is the type of the equatorial radius.  `T` is
the type of the flattening.  The precision of calculations using
this ellipsoid will always be the greater of `A` and `T`.  Hence to
use (for example) lesser precision for calculations, create your
ellipsoid with smaller floats.

# Example
```jldoctest
julia> ell = Ellipsoid(6378_000.f0, Float32(1/298))
Ellipsoid{Float32, Float32}(6.378f6, 0.0033557047f0)

julia> Geodesics.azimuth(Vincenty(ell), 0.0, 0.0, 10.0, 10.0, true) # 32-bit precision
44.7519101709054

julia> ell_big = Ellipsoid(6378_000, 1/BigFloat(298))
Ellipsoid{BigFloat, BigFloat}(6.378e+06, 0.003355704697986577181208053691275167785234899328859060402684563758389261744966452, 6.356597315436241610738255033557046979865771812080536912751677852348993288590605e+06)

julia> Geodesics.azimuth(Vincenty(ell_big), 0, 0, 10, 10, true) # BigFloat precision
44.75207485058905713999544689782065320744100680363085922256125906673345187497678
```

!!! note
    Integer types (e.g., `Int32`, `Int64`) will be promoted to their
    equivalent floating-point types, since geodesic computations
    are performed in floating-point.  Hence if you want to preserve
    the full precision of any integer values, convert them to an adequate
    float type before constructing an `Ellipsoid`.
---

    Ellipsoid(a, f)
    Ellipsoid{A,T}(a::A, f::T)

Create a new ellipsoid with equatorial radius `a` and flattening `f`.
It must be an ellipsoid of rotation; that is, it is rotationally symmetric
about the pole.

By convention, `a` is assumed to be in m if `a` has no units, but
the distances output by [`surface_distance`](@ref) will be in
whatever units are used for `a` (i.e., of type `A`).

Flattening is defined as ``f = (a - b)/a``, where `b` is the polar
radius.  When ``f`` is positive, the ellipsoid is oblate, meaning
it is wider at the equator and 'hamburger-shaped'.
A negative ``f`` gives a prolate ellipsoid like an American football
or rugby ball.

See also: [`Sphere`](@ref)

# Example
```jldoctest
julia> ell = Ellipsoid(6378, 1/298)
Ellipsoid{Float64, Float64}(6378.0, 0.003355704697986577)
```
"""
struct Ellipsoid{A,T}
    "Equatorial radius"
    a::A
    "Flattening"
    f::T
    "Polar radius (private field); stored to avoid recalculating unnecessarily"
    _b::A

    #=
        'Public' internal constructors
    =#
    function Ellipsoid(a, f)
        a′, f′ = promote(float(a), float(f))
        f″, _ = _remove_units(f′)
        _b = _polar_radius(a′, f″)
        Ellipsoid{typeof(a′),typeof(f″)}(a′, f″, _b)
    end


    #=
        'Private' internal constructors
    =#
    function Ellipsoid{A,T}(a, f) where {A,T}
        _b = _polar_radius(a, f)
        new{A,T}(a, f, _b)
    end

    function Ellipsoid{A,T}(a, f, _b) where {A,T}
        a > zero(a) || throw(ArgumentError("equatorial radius must be greater than 0"))
        _b > zero(a) || throw(ArgumentError("polar radius must be greater than 0"))
        abs(f) < one(f) || throw(ArgumentError("flattening must be in range -1 < f < 1"))
        new{A,T}(a, f, _b)
    end
end

"""
    Ellipsoid(; a=$(EARTH_R_MAJOR_WGS84), b=nothing, f=$(EARTH_F_WGS84))

Create a new ellipsoid by providing the equatorial radius `a`, and
one of: the flattening `f`; or polar radius `b`.
"""
function Ellipsoid(; a=EARTH_R_MAJOR_WGS84, b=nothing, f=EARTH_F_WGS84)
    if b !== nothing
        f = _flattening(a, b)
    end
    Ellipsoid(a, f)
end

Base.eltype(::Ellipsoid{A,T}) where {A,T} = T

"""
    Sphere(a=1)

Create a new `Ellipsoid` which represents a sphere with radius `a`.
The type of `a` will be used to construct the `Ellipsoid`.

---
    Sphere(::Type{T}) where T

Create a spherical `Ellipsoid` with unit radius and floating point
type `T`.

# Examples
```jldoctest
julia> Sphere()
Ellipsoid{Float64, Float64}(1.0, 0.0)

julia> Sphere(6371.f0)
Ellipsoid{Float32, Float32}(6371.0f0, 0.0f0)
```
"""
Sphere(a::T=1.0) where T = Ellipsoid(a, zero(T)/oneunit(T))
Sphere(::Type{T}) where T = Ellipsoid(oneunit(T), zero(T)/oneunit(T))
