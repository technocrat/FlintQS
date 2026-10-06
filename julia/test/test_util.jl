using Test, Random
using FlintQS
const F = FlintQS

@testset "util" begin
    @test F.primes_upto(30) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    @test F.primes_upto(1) == Int[]
    @test length(F.SMALL_PRIMES) == 168

    @test F.is_probable_prime(big"2") && F.is_probable_prime(big"97")
    @test !F.is_probable_prime(big"1") && !F.is_probable_prime(big"0")
    @test F.is_probable_prime(big"170141183460469231731687303715884105727")   # 2^127-1
    @test !F.is_probable_prime(big"561")           # Carmichael
    @test !F.is_probable_prime(big"3215031751")    # strong pseudoprime to bases 2,3,5,7
    @test !F.is_probable_prime(big"1000000007" * big"1000000009")

    @test F.integer_root(big"1000000", 2) == 1000
    @test F.integer_root(big"999999", 2) == 999
    @test F.integer_root(big"3"^40, 5) == big"3"^8
    @test F.perfect_power(big"1024") == (2, 10)
    @test F.perfect_power(big"7"^3) == (7, 3)
    @test F.perfect_power(big"1000003"^2) == (1000003, 2)
    @test F.perfect_power(big"1000003" * big"1000033") === nothing
    @test F.perfect_power(big"3") === nothing

    n = big"1000003" * big"1000033"
    d = F.pollard_rho(n; rng = Xoshiro(1))
    @test d !== nothing && 1 < d < n && n % d == 0
    @test F.pollard_rho(big"1000003"; maxiter = 10_000, rng = Xoshiro(1)) === nothing

    @test F.next_prime(big"100") == 101
    @test F.next_prime(big"101") == 101
end
