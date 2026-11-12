# Test against a subset of Karney's dataset
# See tests/data/generate_test_data.jl for details of how the test
# data are obtained, and references therein.

using Test
using Geodesics

"""
Read a single line from `io`, which should contain 10 space-delimited numbers,
and parse as a `T`.

The columns and returned values are:
- latitude at point 1, φ₁ (degrees, exact)
- longitude at point 1, λ₁ (degrees, always 0)
- azimuth at point 1, α₁ (clockwise from north in degrees, exact)
- latitude at point 2, φ₂ (degrees, accurate to 10e−18 deg)
- longitude at point 2, λ₂ (degrees, accurate to 10e−18 deg)
- azimuth at point 2, α₂ (degrees, accurate to 10e−18 deg)
- geodesic distance from point 1 to point 2, s₁₂ (meters, exact)
- arc distance on the auxiliary sphere, σ₁₂ (degrees, accurate to 10e−18 deg)
- reduced length of the geodesic, m₁₂ (meters, accurate to 0.1 pm)
- the area between the geodesic and the equator, S₁₂ (m², accurate to 1 mm²)
"""
function read_line(io::IO, T=Float64)
    line = readline(io)
    ϕ₁, λ₁, α₁, ϕ₂, λ₂, α₂, s₁₂, σ₁₂, m₁₂, S₁₂ = parse.(T, split(line))
    (; ϕ₁, λ₁, α₁, ϕ₂, λ₂, α₂, s₁₂, σ₁₂, m₁₂, S₁₂)
end

"""
Test the 10 numbers on the next line of `io` agree with the values in the
file when using `method`.
"""
function test_line(io, method::Geodesics.GeoDistance; atol_azi=0.01, atol_dist=0.00001, atol_ang_dist=0.001)
    T = eltype(method.ell)
    vals = read_line(io, T)

    ≈(surface_distance(method, vals.λ₁, vals.ϕ₁, vals.λ₂, vals.ϕ₂), vals.s₁₂, atol=atol_dist) ||
        @show :surface_distance vals
    ≈(angular_distance(method, vals.λ₁, vals.ϕ₁, vals.λ₂, vals.ϕ₂), vals.σ₁₂, atol=atol_ang_dist) ||
        @show :angular_distance vals
    ≈(azimuth(method, vals.λ₁, vals.ϕ₁, vals.λ₂, vals.ϕ₂), vals.α₁, atol=atol_azi) ||
        @show :azimuth vals (vals.α₁, azimuth(method, vals.λ₁, vals.ϕ₁, vals.λ₂, vals.ϕ₂))
    ≈(mod(azimuth(method, vals.λ₂, vals.ϕ₂, vals.λ₁, vals.ϕ₁) + 180, 360), vals.α₂, atol=atol_azi) ||
        @show :backazimuth vals (vals.α₂, mod(azimuth(method, vals.λ₂, vals.ϕ₂, vals.λ₁, vals.ϕ₁) + 180, 360))

    @testset "surface_distance" begin
        @test surface_distance(method, vals.λ₁, vals.ϕ₁, vals.λ₂, vals.ϕ₂) ≈ vals.s₁₂ atol=atol_dist
    end
    @testset "angular_distance" begin
        @test angular_distance(method, vals.λ₁, vals.ϕ₁, vals.λ₂, vals.ϕ₂) ≈ vals.σ₁₂ atol=atol_ang_dist
    end
    @testset "azimuth" begin
        @test azimuth(method, vals.λ₁, vals.ϕ₁, vals.λ₂, vals.ϕ₂) ≈ vals.α₁ atol=atol_azi
    end
    @testset "backazimuth" begin
        @test mod(azimuth(method, vals.λ₂, vals.ϕ₂, vals.λ₁, vals.ϕ₁) + 180, 360) ≈ vals.α₂ atol=atol_azi
    end
end

function test_file(file, method::Geodesics.GeoDistance; every=100, atols...)
    open(file) do io
        i = 0
        while !eof(io)
            if i % every == 1
                test_line(io, method; atols...)
            else
                readline(io)
            end
            i += 1
        end
    end
end

function test_all_files(method, files; every=100, atols...)
    @testset "$(basename(file))" for file in files
        test_file(file, method; atols...)
    end
end

# All the data files
data_files = filter(x -> endswith(x, ".txt"), readdir(joinpath(@__DIR__, "data"); join=true))

@testset "Dataset" begin
    @testset "Vincenty $T" for T in (BigFloat, Float32, Float64)
        method = Vincenty(; a=T(6378137), f=1/parse(T, "298.257223563"))
        atol_azi = atol_ang_dist = 0.1
        atol_dist = 0.01
        files = filter(x -> !occursin("antipodal.txt", x), data_files)
        # files = filter(x -> occursin("between_vertices", x), files)
        test_all_files(method, files; atol_azi, atol_ang_dist, atol_dist)
    end
end
