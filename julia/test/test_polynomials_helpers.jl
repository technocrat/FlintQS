using Test, Random
using FlintQS
const F = FlintQS

function make_state(digits; seed = 1)
    rng = Xoshiro(seed)
    p = F.next_prime(BigInt(10)^(digits ÷ 2) + rand(rng, 1:10^6))
    q = F.next_prime(BigInt(10)^(digits - digits ÷ 2) + rand(rng, 1:10^6))
    n = p * q
    k = F.knuth_schroeppel(n)
    kn = k * n
    P = F.params_for(ndigits(n))
    fb = F.build_factorbase(kn, k, P.nprimes)
    return F.SIQSState(kn, fb, P), rng
end
