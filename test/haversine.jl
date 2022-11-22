# Haversine method

using Geodesics
using Test
using Unitful: @u_str

float_types = (Float16, Float32, Float64, BigFloat)
a_float_types = (float_types..., typeof(1.0u"m"))

"""
Check that the type of the ellipsoid created with the arguments passed to `Haversine`
are the same as for the constructor itself.
"""
function test_type_parameters(::Haversine{Av,Tv}, ::Ellipsoid{Ae,Te}) where {Av,Tv,Ae,Te}
    @test Av == Ae
    @test Tv == Te
end

function random_points(; radians=false)
    lon0 = 360*rand() - 180
    lat0 = 180*rand() - 90
    lon1 = 360*rand() - 180
    lat1 = 180*rand() - 90
    if radians
        lon0, lat0, lon1, lat1 = deg2rad.((lon0, lat0, lon1, lat1))
    end
    (; lon0, lat0, lon1, lat1)
end


@testset "Haversine" begin
    # Default method to use
    alg = Haversine()

    @testset "Constructors" begin
        @test Haversine() isa Haversine
        @test Haversine(2.f0).ell.a == 2.f0

        @testset "A: $A" for A in a_float_types
            @testset "T: $T" for T in float_types
                test_type_parameters(Haversine(; a=one(A), f=zero(T)),
                    Ellipsoid(one(A), zero(T)))
                test_type_parameters(Haversine(; a=oneunit(A), f=zero(T)),
                    Ellipsoid(a=oneunit(A), f=zero(T)))
                test_type_parameters(Haversine(; a=oneunit(A), b=oneunit(A)),
                    Ellipsoid(a=oneunit(A), b=oneunit(A)))
            end
        end

        @testset "Errors" begin
            @test_throws ArgumentError Haversine(-1)
            @test_throws ArgumentError Haversine(0)
            @test_throws ArgumentError Haversine(; a=1, f=1)
        end
    end

    @testset "forward" begin
        @testset "Wrong latitude" begin
            @test_throws ArgumentError Geodesics.forward(alg, 0, 2π, 0, 0)
        end

        @testset "No check on latitude" begin
            @test Geodesics.forward(alg, 0, 2π, 0, 0; check=false) isa NamedTuple
        end

        # Compare to known values on the unit sphere
        @testset "Known answer" begin
            result = Geodesics.forward(Haversine(1), 0, 0, 1, 1)
            @test result.lon ≈ 0.9189892552937267
            @test result.lat ≈ 0.4719777676633856
            @test result.backazimuth ≈ mod2pi(-1.9047280216911358) atol=1e-6
        end

        # Make sure if we don't move, we end up in the same place
        @testset "Zero distance" begin
            lon0 = 0.1
            lat0 = 0.2
            azimuth = 1
            distance = 0
            result = Geodesics.forward(Haversine(1), lon0, lat0, azimuth, distance)
            @test result.lon == lon0
            @test result.lat == lat0
            @test result.backazimuth ≈ mod2pi(azimuth + π)
        end
    end

    @testset "inverse" begin
        @testset "Wrong latitude" begin
            @test_throws ArgumentError Geodesics.inverse(Haversine(), 0, 0, 0, -3π)
            @test_throws ArgumentError Geodesics.inverse(Haversine(), 0, 2π, 0, 0)
        end

        @testset "Versus Vincenty on sphere" begin
            # Compare to Vincenty computation on the unit sphere
            result = Geodesics.inverse(Haversine(1), 0, 0, 1, 1)
            atol = 1e-6
            @test result.surface_distance ≈ 1.2745557823062943 atol=atol
            @test result.azimuth ≈ 0.49536728921867335 atol=atol
            @test result.backazimuth ≈ 4.217021691166017 atol=atol
        end
    end

    @testset "angular_distance" begin
        @testset "Degrees" begin
            lon0, lat0, lon1, lat1 = random_points()
            atol = 1e-6

            @testset "Sphere" begin
                method = Haversine()
                Δ₁ = Geodesics.angular_distance(method, lon0, lat0, lon1, lat1)
                Δ₂ = rad2deg(Geodesics.inverse(method, deg2rad.((lon0, lat0, lon1, lat1))...).angular_distance)
                @test Δ₁ ≈ Δ₂ atol=atol
            end
        end

        # Radians
        @testset "Radians" begin
            lon0, lat0, lon1, lat1 = random_points(radians=true)
            atol = 1e-6

            @testset "Sphere" begin
                @test Geodesics.angular_distance(Haversine(), lon0, lat0, lon1, lat1; degrees=false) ≈
                    Geodesics.inverse(Haversine(a=1, f=0),
                        lon0, lat0, lon1, lat1).angular_distance atol=atol
            end
        end
    end

    @testset "azimuth" begin
        @testset "Degrees" begin
            lon0, lat0, lon1, lat1 = random_points()
            atol = 1e-6

            @test Geodesics.azimuth(Haversine(f=0), lon0, lat0, lon1, lat1) ≈
                rad2deg(Geodesics.inverse(Haversine(a=1, f=0),
                    deg2rad.((lon0, lat0, lon1, lat1))...).azimuth) atol=atol
        end

        @testset "Radians" begin
            lon0, lat0, lon1, lat1 = random_points(radians=true)
            atol = 1e-6

            @test Geodesics.azimuth(Haversine(f=0), lon0, lat0, lon1, lat1; degrees=false) ≈
                Geodesics.inverse(Haversine(a=1, f=0), lon0, lat0, lon1, lat1).azimuth atol=atol
        end

        @testset "One point at pole" begin
            lon0 = 0
            lat0 = -90
            lon1 = 0
            lat1 = 0

            for T in (Int16, Int32, Int64, Float16, Float32, Float64)
                @test Geodesics.azimuth(Haversine(), T(lon0), T(lat0), T(lon1), T(lat1); degrees=true) == 0.0
            end

            lon0 = 90
            lat0 = 0
            lon1 = 90
            lat1 = 90

            for T in (Int16, Int32, Int64, Float16, Float32, Float64)
                @test Geodesics.azimuth(Haversine(), T(lon0), T(lat0), T(lon1), T(lat1); degrees=true) == 0.0
            end
        end
    end

    @testset "angular_step" begin
        @testset "Degrees" begin
            lon0, lat0, lon1, lat1 = random_points()
            dist = 360rand()
            az = 360rand() - 90
            atol = 1e-6

            @testset "Sphere" begin
                lon, lat, baz = Geodesics.angular_step(Haversine(a=1, f=0), lon0, lat0, az, dist)
                result = Geodesics.forward(Haversine(a=1, f=0), deg2rad.((lon0, lat0, az, dist))...)
                @test lon ≈ rad2deg(result.lon) atol=atol
                @test lat ≈ rad2deg(result.lat) atol=atol
                @test baz ≈ rad2deg(result.backazimuth) atol=atol
            end
        end

        @testset "Radians" begin
            lon0, lat0, lon1, lat1 = random_points(radians=true)
            dist = 2π*rand()
            az = 2π*rand() - π/2
            atol = 1e-6

            lon, lat, baz = Geodesics.angular_step(Haversine(a=1, f=0),
                lon0, lat0, az, dist; degrees=false)
            result = Geodesics.forward(Haversine(a=1, f=0), lon0, lat0, az, dist)
            @test lon ≈ result.lon atol=atol
            @test lat ≈ result.lat atol=atol
            @test baz ≈ result.backazimuth atol=atol
        end
    end

    @testset "surface_distance" begin
        @testset "Degrees" begin
            lon0, lat0, lon1, lat1 = random_points()
            a = 10_000_000*rand()
            atol = 1e-6
            @test Geodesics.surface_distance(Haversine(a=a, f=0), lon0, lat0, lon1, lat1) ≈
                Geodesics.inverse(Haversine(a=a, f=0), deg2rad.((lon0, lat0, lon1, lat1))...).surface_distance atol=atol
        end

        @testset "Radians" begin
            lon0, lat0, lon1, lat1 = random_points(radians=true)
            a = 10_000_000*rand()
            atol = 1e-6

            @test Geodesics.surface_distance(Haversine(a=a, f=0), lon0, lat0, lon1, lat1; degrees=false) ≈
                Geodesics.inverse(Haversine(a=a, f=0), lon0, lat0, lon1, lat1).surface_distance atol=atol
        end
    end

    @testset "surface_step" begin
        @testset "Degrees" begin
            lon0, lat0, lon1, lat1 = random_points()
            dist = 360rand()
            az = 360rand() - 90
            atol = 1e-6

            @testset "Sphere" begin
                lon, lat, baz = Geodesics.surface_step(Haversine(), lon0, lat0, az, dist)
                result = Geodesics.forward(Haversine(), deg2rad.((lon0, lat0, az))..., dist)
                @test lon == rad2deg(result.lon)
                @test lat == rad2deg(result.lat)
                @test baz == rad2deg(result.backazimuth)
            end
        end

        @testset "Radians" begin
            lon0, lat0, lon1, lat1 = random_points(radians=true)
            dist = rand()
            az = 2π*rand() - π/2
            atol = 1e-6

            lon, lat, baz = Geodesics.surface_step(Haversine(a=1),
                lon0, lat0, az, dist; degrees=false)
            result = Geodesics.forward(Haversine(a=1), lon0, lat0, az, dist)
            @test lon ≈ result.lon atol=atol
            @test lat ≈ result.lat atol=atol
            @test baz ≈ result.backazimuth atol=atol
        end

        @testset "Units" begin
            lon0, lat0, _, _ = random_points(radians=true)
            az = 2π*rand() - π/2
            ell = Sphere(Geodesics.EARTH_R_MAJOR_WGS84/1000*u"km")
            alg = Haversine(ell)

            @testset "km and km" begin
                dist = 10_000*u"km"
                lon, lat, baz = Geodesics.surface_step(alg, lon0, lat0, az, dist; degrees=false)
                result = Geodesics.forward(alg, lon0, lat0, az, dist)
                @test lon == result.lon
                @test lat == result.lat
                @test baz == result.backazimuth
            end

            @testset "m and km" begin
                dist = 10_000_000*u"m"
                lon, lat, baz = Geodesics.surface_step(alg, lon0, lat0, az, dist; degrees=false)
                result = Geodesics.forward(alg, lon0, lat0, az, dist)
                @test lon == result.lon
                @test lat == result.lat
                @test baz == result.backazimuth
            end
        end
    end
end
