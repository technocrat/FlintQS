using Test, Random
using FlintQS
const F = FlintQS

@testset "modarith" begin
    # jacobi vs Euler's criterion for prime moduli
    for p in (3, 5, 7, 11, 13, 101, 997)
        for a in 1:p-1
            e = powermod(a, (p - 1) ÷ 2, p)
            @test F.jacobi(a, p) == (e == 1 ? 1 : -1)
        end
        @test F.jacobi(0, p) == 0
    end
    @test F.jacobi(2, 15) == 1       # composite modulus
    @test F.jacobi(5, 9) == 1
    @test F.jacobi(big"123456789012345678901234567890", 1009) ==
          F.jacobi(Int(mod(big"123456789012345678901234567890", 1009)), 1009)
    @test_throws ArgumentError F.jacobi(3, 8)

    # sqrt_mod_prime: p ≡ 3 mod 4, p ≡ 1 mod 4 (incl. p ≡ 1 mod 8, 17, 97 with large 2-power)
    for p in (3, 7, 11, 13, 17, 41, 97, 193, 1009, 786433, 1000003)
        for a in unique(rand(Xoshiro(p), 0:p-1, 50))
            if a == 0 || F.jacobi(a, p) == 1
                r = F.sqrt_mod_prime(a, p)
                @test mod(r * r, p) == a
            end
        end
    end
    @test F.sqrt_mod_prime(0, 13) == 0
end
