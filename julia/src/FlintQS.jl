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
    # ≥ 31 digits. Only divisors that split `m` are accepted, so the lift prime itself is skipped
    # without discarding the sieve run.
    pdig = max(31 - d, 4)
    P = next_prime(BigInt(10)^(pdig - 1) + rand(rng, 1:10^min(pdig - 2, 15)))
    g = siqs_split(m * P; rng = rng, verbose = verbose, accept = h -> (c = gcd(h, m); 1 < c < m))
    return gcd(g, m)
end

"""
    flintqs(n; rng, verbose) -> Vector{BigInt}

Prime factorization of `n > 1` as a sorted vector with multiplicity (`[]` for `n == 1`).
Strips factors below 1000, reduces perfect powers, and splits the rest with Pollard rho
(≤ 18 digits) or the self-initializing quadratic sieve.
"""
function flintqs(n::Integer; rng = Random.default_rng(), verbose::Bool = false)
    n > 0 || throw(ArgumentError("flintqs requires n > 0, got $n"))
    n == 1 && return BigInt[]
    m = BigInt(n)
    out = BigInt[]
    for p in SMALL_PRIMES
        while m % p == 0
            push!(out, BigInt(p))
            m ÷= p
        end
    end
    stack = BigInt[m]
    while !isempty(stack)
        c = pop!(stack)
        c == 1 && continue
        if is_probable_prime(c)
            push!(out, c)
        elseif (pp = perfect_power(c)) !== nothing
            b, k = pp
            for _ in 1:k
                push!(stack, b)
            end
        else
            d = split_composite(c; rng = rng, verbose = verbose)
            push!(stack, d, c ÷ d)
        end
    end
    return sort!(out)
end

end # module




