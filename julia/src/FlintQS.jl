# SPDX-License-Identifier: GPL-2.0-or-later
"""
    FlintQS

Julia port of William Hart's FlintQS: a self-initializing quadratic sieve with single
large-prime merging and Block Lanczos linear algebra.

Exports [`flintqs`](@ref) (full prime factorization) and [`siqs_split`](@ref) (one divisor
of a ≥ 40-digit composite).
"""
module FlintQS

using Random
include("util.jl")
include("modarith.jl")
include("params.jl")
include("factorbase.jl")
include("polynomials.jl")
include("sieve.jl")
include("relations.jl")
include("lanczos.jl")

export flintqs, siqs_split

"Sieve with fresh A values until the store holds at least `target` relations."
function collect_relations!(store::RelationStore, st::SIQSState, target::Int, rng, verbose::Bool)
    sv = zeros(UInt8, 2 * st.params.M)
    curves = 0
    while nrelations(store) < target
        choose_a!(st, rng)
        init_a!(st)
        for i in 0:(2^(st.s - 1) - 1)
            i > 0 && next_poly!(st, i)
            sieve!(sv, st)
            scan!(store, st, sv)
            curves += 1
        end
        # fulls $(length(store.fulls)), combined $(length(store.combined)), partials $(store.npartials)
        verbose && print("$(nrelations(store)) of $target relations\r")
    end
    return store
end

"Square-root step for one dependency; returns a nontrivial divisor of `n` or `nothing`."
function try_dependency(rels::Vector{Relation}, dep::BitVector, fb::FactorBase, kn::BigInt, n::BigInt, accept)
    counts = zeros(Int, length(fb.primes) + 1)
    X = BigInt(1)
    for j in findall(dep)
        X = mod(X * rels[j].x, kn)
        for r in rels[j].fac
            counts[r] += 1
        end
    end
    any(isodd, counts) && return nothing
    Y = BigInt(1)
    for row in 2:length(counts)
        c = counts[row]
        c == 0 && continue
        Y = mod(Y * powermod(BigInt(fb.primes[row-1]), c ÷ 2, kn), kn)
    end
    for g in (gcd(X - Y, n), gcd(X + Y, n)), h in (g, n ÷ g)
        1 < h < n && accept(h) && return h
    end
    return nothing
end

function dependencies_to_factor(store::RelationStore, fb::FactorBase, kn::BigInt, n::BigInt, rng, accept)
    rels = vcat(store.fulls, store.combined)
    nrows = length(fb.primes) + 1
    cols = [parity_rows(r.fac) for r in rels]
    keep, nonempty = prune_singletons(cols, nrows)
    length(keep) - nonempty < 32 && return nothing
    kept_rels = rels[keep]
    kept_cols = cols[keep]
    for _ in 1:5
        for dep in block_lanczos(kept_cols, nrows; rng = rng)
            d = try_dependency(kept_rels, dep, fb, kn, n, accept)
            d !== nothing && return d
        end
    end
    return nothing
end

"""
    siqs_split(n; rng, verbose, accept) -> BigInt

Return a nontrivial divisor of the composite `n` using the self-initializing quadratic sieve.
`n` must have at least 40 decimal digits, no prime factor below 1000 and not be a perfect power.
`accept(d)` filters which divisors may be returned (default: any); the search continues through
the remaining dependencies, and further sieving, until one is accepted.
"""
function siqs_split(n::BigInt; rng = Random.default_rng(), verbose::Bool = false, accept = _ -> true)
    if verbose

    end
    verbose && println("Choosing Parameters...")
    P = params_for(ndigits(n))
    verbose && println("Choosing Multiplier...")
    k = knuth_schroeppel(n)
    kn = k * n
    verbose && println("Building Factor Base...")
    fb = build_factorbase(kn, k, P.nprimes)
    st = SIQSState(kn, fb, P)
    store = RelationStore(kn)
    target = P.nprimes + 64
    while true
        verbose && println("Sieving...")
        collect_relations!(store, st, target, rng, verbose)
        for g in store.found, h in (g, n ÷ g)
            1 < h < n && accept(h) && return h
        end
        verbose && println("\nExtracting Factors...")
        d = dependencies_to_factor(store, fb, kn, n, rng, accept)
        d !== nothing && return d
        verbose && println("Sieving Additional Relations...")
        target += 64
    end
end

"Nontrivial divisor of a composite, non-power `m` with no prime factor below 1000."
function split_composite(m::BigInt; rng, verbose::Bool = false)
    d = ndigits(m)
    if d <= 18
        g = pollard_rho(m; rng = rng)
        g === nothing && error("Pollard rho failed on a ≤18-digit number $m")
        return g
    elseif d >= 30
        return siqs_split(m; rng = rng, verbose = verbose)
    end
    # 19–29 digits: multiply by a prime (≥ 4 digits, so no factor < 1000) so the work number has
    # ≥ 41 digits. Only divisors that split `m` are accepted, so the lift prime itself is skipped
    # without discarding the sieve run.
    pdig = max(31 - d, 4)
    P = next_prime(BigInt(10)^(pdig - 1) + rand(rng, 1:10^min(pdig - 2, 15)))
    g = siqs_split(m * P; rng = rng, verbose = verbose, accept = h -> (c = gcd(h, m); 1 < c < m))
    return gcd(g, m)
end

"""
    flintqs(n; rng, verbose) -> Vector{BigInt}

Prime factorization of `n ≥ 1` as a sorted tuple of length 2,
Strips factors below 1000 and splits the rest with Pollard rho
(≤ 18 digits) or the self-initializing quadratic sieve.
Catches input that is prime or an integer power.
"""
function flintqs(n::Integer; rng = Random.default_rng(), verbose::Bool = false)
    n >= 1 || throw(ArgumentError("flintqs requires n ≥ 1, got $n"))
    m = BigInt(n)
    # catch input values that will cause quadratic sieve to fail, prime number and integer power
    is_probable_prime(m) && return (BigInt(1), BigInt(m))
    pp = perfect_power(m)
    pp !== nothing && return (pp[1], m ÷ pp[1])
    for p in SMALL_PRIMES
        while m % p == 0
            return (BigInt(p), m ÷ p)
        end
    end
    d = split_composite(m; rng = rng, verbose = verbose)
    sort((d, m ÷ d))
end

end # module

"Test a single number"
function test1()
    # some test numbers, uncomment one
    #n = "6"
    #n = "3344026937 " # 10 digits
    #n = "81937514831734261501" # 20 digits
    #n = "490516739493812924575428019499" # 30 digits
    #n = "3009945328299388957931300315646308313419" # 40 digits
    #n = "32466092539065621597006224985875602441118457873449" # 50 digits
    n = "750014700584359514878998176863189936656347267709304336104953" # 60 digits
    #n = "3821588930350792167954816009780238547235939701293631242896458394309323" # 70 digits
    #n = "444165670830100364807708934724975053225119501482759395905805810135248970251" # 75 digits
    #n = "32359223112252035483026442913476041867810505843753089370355888993328247493431081" # 80 digits
    #n = "210282352803628314297444623664410288223068234485244872675057323275627391582812832505918087" # 90 digits
    #n = "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139" # rsa-100
    #n = "35794234179725868774991807832568455403003778024228226193532908190484670252364677411513516111204504060317568667" # rsa-110
    nn = parse(BigInt, length(ARGS) > 0 ? ARGS[1] : n)
    start_time = time()
    factors = FlintQS.flintqs(nn, verbose=true)
    println("elapsed time = ", time() - start_time)
    println(factors)
end

include("util.jl")

"Test random semiprimes of length 2-42"
function test2()
    for i in 2:21
        println("i = $i")
        for _ in 1:200
            while true
                k = randomprime(i-1)
                l = randomprime(i-1)
                k != l && break
            end
            m = k * l
            f0, f1 = FlintQS.flintqs(m, verbose=false)
            if m != f0 * f1
                println("test2() failed for $k, $m")
                exit(0)
            end
            k = randomprime(i)
            l = randomprime(i-1)
            m = k * l
            f0, f1 = FlintQS.flintqs(m, verbose=false)
            if m != f0 * f1
                println("test2() failed for $k, $m")
                exit(0)
            end
            while true
                k = randomprime(i)
                l = randomprime(i)
                k != l && break
            end
            m = k * l
            f0, f1 = FlintQS.flintqs(m, verbose=false)
            if m != f0 * f1
                println("test2() failed for $k, $m")
                exit(0)
            end
        end
    end
    println("test2 passed")
end

"Test prime and prime power"
function test3()
    f0, f1 = FlintQS.flintqs(BigInt(11), verbose=false)
    if (f0, f1) != (BigInt(1), BigInt(11))
        println("flintqs(prime) failed")
        exit(0)
    end
    f0, f1 = FlintQS.flintqs(BigInt(8), verbose=false)
    if (f0, f1) != (BigInt(2), BigInt(4))
        println("flintqs(prime_power) failed")
        exit(0)
    end
    println("test3 passed")
end

#test1()
#test2()
test3()




