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
