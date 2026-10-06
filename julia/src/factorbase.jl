# SPDX-License-Identifier: GPL-2.0-or-later

const MULTIPLIERS = (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)

struct FactorBase
    primes::Vector{Int}
    sizes::Vector{UInt8}
    sqrts::Vector{Int}
    multiplier::Int
end

"Knuth-Schroeppel multiplier: favour k for which k*n is a residue modulo many small primes."
function knuth_schroeppel(n::BigInt)
    best, bestscore = 1, -Inf
    for k in MULTIPLIERS
        kn = k * n
        score = -0.5 * log(k)
        r8 = Int(mod(kn, 8))
        score += r8 == 1 ? 2 * log(2) : r8 == 5 ? log(2) : (r8 == 3 || r8 == 7) ? 0.5 * log(2) : 0.0
        for p in SMALL_PRIMES
            p == 2 && continue
            if k % p == 0
                score += log(p) / p
            elseif jacobi(kn, p) == 1
                score += 2 * log(p) / p
            end
        end
        if score > bestscore
            best, bestscore = k, score
        end
    end
    return best
end

"First `nprimes` primes p (including 2) for which `kn` is a square mod p, with roots and bit sizes."
function build_factorbase(kn::BigInt, k::Int, nprimes::Int)
    limit = max(1000, ceil(Int, 2.5 * nprimes * (log(nprimes) + 2)))
    while true
        primes = Int[]
        sqrts = Int[]
        for p in primes_upto(limit)
            if p == 2
                push!(primes, 2)
                push!(sqrts, Int(mod(kn, 2)))
            else
                a = Int(mod(kn, p))
                (a == 0 || jacobi(a, p) == 1) || continue
                push!(primes, p)
                push!(sqrts, a == 0 ? 0 : sqrt_mod_prime(a, p))
            end
            length(primes) == nprimes && break
        end
        if length(primes) == nprimes
            sizes = UInt8[UInt8(floor(Int, log2(p) - 0.15 + 0.5)) for p in primes]
            return FactorBase(primes, sizes, sqrts, k)
        end
        limit *= 2
    end
end
