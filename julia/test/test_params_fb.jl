using Test, Random
using FlintQS
const F = FlintQS

@testset "params" begin
    p = F.params_for(40)
    @test (p.nprimes, p.M, p.firstprime, p.errorbits, p.threshold) == (1500, 32000, 8, 16, 66)
    @test F.params_for(60).nprimes == 3000
    @test F.params_for(91).nprimes == 80000
    @test F.params_for(91).threshold == 102
    @test F.params_for(92).nprimes == 64000
    @test_throws ArgumentError F.params_for(39)
    # table lengths are consistent (52 rows, digits 40..91)
    for d in 40:91
        q = F.params_for(d)
        @test q.largeprime > 1000 && q.threshold > q.errorbits
    end
end

@testset "factor base" begin
    n = big"1000003" * big"1000033" * big"999983"
    k = F.knuth_schroeppel(n)
    @test k in (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    kn = k * n
    fb = F.build_factorbase(kn, k, 200)
    @test length(fb.primes) == 200 && fb.primes[1] == 2
    @test issorted(fb.primes) && allunique(fb.primes)
    for (i, p) in enumerate(fb.primes)
        p == 2 && continue
        @test mod(fb.sqrts[i]^2, p) == mod(kn, p)
        @test F.is_probable_prime(BigInt(p))
    end
    @test fb.sizes[2] == round(Int, log2(fb.primes[2]) - 0.15)
    @test F.knuth_schroeppel(big"1") in (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
end
