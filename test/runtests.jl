using Geodesics
using Test

@testset "Forward" begin
    @test_throws ArgumentError Geodesics.forward(0, 2π, 0, 0, 1, 0)
    @test_throws ArgumentError Geodesics.forward(0, 2π, 0, 0, 0, 0)
    @test_throws ArgumentError Geodesics.inverse(0, 0, 1, 1, 1, 1)

    # Compare to Haversine computation on the unit sphere
    let lon0 = 0, lat0 = 0, azimuth = 1, distance = 1, a = 1, f = 0
        lon1, lat1, baz = Geodesics.forward(lon0, lat0, azimuth, distance, a, f)
        @test lon1 ≈ 0.9189892552937267
        @test lat1 ≈ 0.4719777676633856
        @test baz ≈ mod(-1.9047280216911358, 2π) atol=1e-6
    end

    # Make sure if we don't move, we end up in the same place
    let lon0 = 0.1, lat0 = 0.2, azimuth = 1, distance = 0, a = 1, f = 0
        lon1, lat1, baz = Geodesics.forward(lon0, lat0, azimuth, distance, a, f)
        @test lon1 == lon0
        @test lat1 == lat0
        @test baz ≈ mod(azimuth + π, 2π)
    end
end

@testset "Inverse" begin
    @test_throws ArgumentError Geodesics.inverse(0, 0, 1, 1, 0, 0)
    @test_throws ArgumentError Geodesics.inverse(0, 2π, 0, -3π, 1, 1)
    @test_throws ArgumentError Geodesics.inverse(0, 0, 1, 1, 1, 1)

    # Compare to Haversine computation on the unit sphere
    let lon0 = 0, lat0 = 0, lon1 = 1, lat1 = 1, a = 1, f = 0, atol = 1e-6
        distance, azimuth, backazimuth = Geodesics.inverse(lon0, lat0, lon1, lat1, a, f)
        @test distance ≈ 1.2745557823062943 atol=atol
        @test azimuth ≈ 0.49536728921867335 atol=atol
        @test backazimuth ≈ mod(-2.0661636160135703, 2π) atol=atol
    end
end

@testset "Ang dist" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, lon1 = 360*rand() - 180,
            lat1 = 180*rand() - 90, atol = 1e-6
        @test Geodesics.angular_distance(lon0, lat0, lon1, lat1, f=0) ≈
            rad2deg(Geodesics.inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, 0)[1]) atol=atol
        @test Geodesics.angular_distance(lon0, lat0, lon1, lat1) ≈
            rad2deg(Geodesics.inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, Geodesics.F_WGS84)[1]) atol=atol
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), lon1 = 2π*(rand() - 0.5),
            lat1 = π*(rand() - 0.5), atol = 1e-6
        @test Geodesics.angular_distance(lon0, lat0, lon1, lat1, false, f=0) ≈
            Geodesics.inverse(lon0, lat0, lon1, lat1, 1, 0)[1] atol=atol
        @test Geodesics.angular_distance(lon0, lat0, lon1, lat1, false) ≈
            Geodesics.inverse(lon0, lat0, lon1, lat1, 1, Geodesics.F_WGS84)[1] atol=atol
    end
end

@testset "Azimuth" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, lon1 = 360*rand() - 180,
            lat1 = 180*rand() - 90, atol = 1e-6
        @test Geodesics.azimuth(lon0, lat0, lon1, lat1, f=0) ≈
            rad2deg(Geodesics.inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, 0)[2]) atol=atol
        @test Geodesics.azimuth(lon0, lat0, lon1, lat1) ≈
            rad2deg(Geodesics.inverse(deg2rad.((lon0, lat0, lon1, lat1))..., 1, Geodesics.F_WGS84)[2]) atol=atol
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), lon1 = 2π*(rand() - 0.5),
            lat1 = π*(rand() - 0.5), atol = 1e-6
        @test Geodesics.azimuth(lon0, lat0, lon1, lat1, false, f=0) ≈
            Geodesics.inverse(lon0, lat0, lon1, lat1, 1, 0)[2] atol=atol
        @test Geodesics.azimuth(lon0, lat0, lon1, lat1, false) ≈
            Geodesics.inverse(lon0, lat0, lon1, lat1, 1, Geodesics.F_WGS84)[2] atol=atol
    end

    # One point at pole
    let lon0 = 0, lat0 = -90, lon1 = 0, lat1 = 0
        for T in (Int16, Int32, Int64, Float16, Float32, Float64)
            @test Geodesics.azimuth(T(lon0), T(lat0), T(lon1), T(lat1), true) == 0.0
        end
    end
    let lon0 = 90, lat0 = 0, lon1 = 90, lat1 = 90
        for T in (Int16, Int32, Int64, Float16, Float32, Float64)
            @test Geodesics.azimuth(T(lon0), T(lat0), T(lon1), T(lat1), true) == 0.0
        end
    end
end

@testset "Ang step" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, dist = 360*rand(),
            az = 360*rand() - 90, atol = 1e-6
        @test all(isapprox.(Geodesics.angular_step(lon0, lat0, az, dist, f=0),
            rad2deg.(Geodesics.forward(deg2rad.((lon0, lat0, az, dist))..., 1, 0)), atol=atol))
        @test all(isapprox.(Geodesics.angular_step(lon0, lat0, az, dist),
            rad2deg.(Geodesics.forward(deg2rad.((lon0, lat0, az, dist))..., 1, Geodesics.F_WGS84)), atol=atol))
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), dist = 2π*rand(),
            az = 2π*(rand() - 0.5), atol = 1e-6
        @test all(isapprox.(Geodesics.angular_step(lon0, lat0, az, dist, false, f=0),
            Geodesics.forward(lon0, lat0, az, dist, 1, 0), atol=atol))
        @test all(isapprox.(Geodesics.angular_step(lon0, lat0, az, dist, false),
            Geodesics.forward(lon0, lat0, az, dist, 1, Geodesics.F_WGS84), atol=atol))
    end
end

@testset "Surf dist" begin
    # Degrees
    let lon0 = 360*rand() - 180, lat0 = 180*rand() - 90, lon1 = 360*rand() - 180,
            lat1 = 180*rand() - 90, a = 10_000_000*rand(), atol = 1e-6
        @test Geodesics.surface_distance(lon0, lat0, lon1, lat1, a, f=0) ≈
            Geodesics.inverse(deg2rad.((lon0, lat0, lon1, lat1))..., a, 0)[1] atol=atol
        @test Geodesics.surface_distance(lon0, lat0, lon1, lat1, a) ≈
            Geodesics.inverse(deg2rad.((lon0, lat0, lon1, lat1))..., a, Geodesics.F_WGS84)[1] atol=atol
    end

    # Radians
    let lon0 = 2π*(rand() - 0.5), lat0 = π*(rand() - 0.5), lon1 = 2π*(rand() - 0.5),
            lat1 = π*(rand() - 0.5), a = 10_000_000*rand(), atol = 1e-6
        @test Geodesics.surface_distance(lon0, lat0, lon1, lat1, a, false, f=0) ≈
            Geodesics.inverse(lon0, lat0, lon1, lat1, a, 0)[1] atol=atol
        @test Geodesics.surface_distance(lon0, lat0, lon1, lat1, a, false) ≈
            Geodesics.inverse(lon0, lat0, lon1, lat1, a, Geodesics.F_WGS84)[1] atol=atol
    end
end
