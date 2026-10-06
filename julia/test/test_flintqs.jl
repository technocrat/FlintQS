using Test, Random
using FlintQS
const F = FlintQS

function semiprime(rng, d1, d2)
    p = F.next_prime(BigInt(10)^(d1 - 1) + rand(rng, 1:10^6))
    q = F.next_prime(BigInt(10)^(d2 - 1) + rand(rng, 1:10^6))
    p == q ? semiprime(rng, d1, d2) : (p * q, p, q)
end

@testset "siqs_split" begin
    rng = Xoshiro(2024)
    for (d1, d2) in ((20, 21), (25, 26))        # 41- and 51-digit semiprimes
        n, p, q = semiprime(rng, d1, d2)
        d = F.siqs_split(n; rng = rng)
        @test 1 < d < n && n % d == 0
        @test d in (p, q)
    end
end

@testset "flintqs" begin
    rng = Xoshiro(99)
    # trivial inputs
    @test F.flintqs(1) == BigInt[]
    @test F.flintqs(97) == [97]
    @test_throws ArgumentError F.flintqs(0)
    @test_throws ArgumentError F.flintqs(-5)
    # small-factor stripping and ordinary factorization
    @test F.flintqs(2^5 * 3^3 * 7 * 1000003; rng = rng) == BigInt.([2, 2, 2, 2, 2, 3, 3, 3, 7, 1000003])
    # <= 18 digits (rho)
    n, p, q = semiprime(rng, 8, 9)
    @test F.flintqs(n; rng = rng) == sort([p, q])
    # Review Focus 2: prime powers
    P = F.next_prime(big"10"^24)
    @test F.flintqs(P^2; rng = rng) == [P, P]
    @test F.flintqs(P^3; rng = rng) == [P, P, P]
    # Review Focus 3: many small factors times a large prime
    @test F.flintqs(big"2"^5 * big"3"^3 * P; rng = rng) == vcat(fill(BigInt(2), 5), fill(BigInt(3), 3), [P])
    # Review Focus 4: 39-digit (lift path) and 40-digit (direct) semiprimes
    for (d1, d2) in ((19, 20), (20, 20))
        n, p, q = semiprime(rng, d1, d2)
        @test F.flintqs(n; rng = rng) == sort([p, q])
    end
    # Review Focus 5: very unbalanced semiprime (7-digit × 45-digit)
    n, p, q = semiprime(rng, 7, 45)
    @test F.flintqs(n; rng = rng) == sort([p, q])
    # product of three primes, 50 digits
    p1 = F.next_prime(big"10"^16 + 12345); p2 = F.next_prime(big"10"^16 + 99999); p3 = F.next_prime(big"10"^17 + 777)
    @test F.flintqs(p1 * p2 * p3; rng = rng) == sort([p1, p2, p3])
end

@testset "siqs_split accept predicate (lift path)" begin
    rng = Xoshiro(314)
    m, p, q = semiprime(rng, 19, 20)               # 39-digit m
    P = F.next_prime(big"1000")                    # 4-digit lift prime: in the factor base, sqrt 0
    N = m * P
    accept = h -> (c = gcd(h, m); 1 < c < m)
    for _ in 1:5                                   # P-only divisors must never be returned
        d = F.siqs_split(N; rng = rng, accept = accept)
        @test N % d == 0
        @test gcd(d, m) in (p, q)
    end
end
