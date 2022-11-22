using Geodesics
using Test
using Unitful: @u_str

float_types = (Float16, Float32, Float64, BigFloat)

@testset "Types" begin
    @testset "Ellipsoid" begin
        @testset "Constructors" begin
            @testset "A $A" for A in (float_types..., typeof(1.0u"m"))
                @testset "T $T" for T in float_types
                    let a = rand(A), f = rand(T), ell = Ellipsoid(a, f)
                        if A == T
                            a_, f_ = a, f
                            A_, T_ = A, T
                        else
                            a_nounit, f_ = promote(a/oneunit(a), f)
                            a_ = a_nounit*oneunit(a)
                            A_, T_ = typeof.((a_, f_))
                        end
                        @test ell.a == a_
                        @test ell.f == f_
                        @test ell.a isa A_
                        @test ell.f isa T_
                        @test ell isa Ellipsoid{A_,T_}
                    end
                end
            end

            @testset "Default" begin
                let ell = Ellipsoid()
                    @test ell isa Ellipsoid{Float64,Float64}
                    @test ell.a == Geodesics.EARTH_R_MAJOR_WGS84
                    @test ell.f == Geodesics.EARTH_F_WGS84
                end
            end
        end
    end
end
