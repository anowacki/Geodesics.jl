using Geodesics
using Test

@testset "All tests" begin
    include("types.jl")
    include("vincenty.jl")
    include("haversine.jl")
end
