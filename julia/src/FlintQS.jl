# SPDX-License-Identifier: GPL-2.0-or-later
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
        verbose && println("curves=$curves full=$(length(store.fulls)) combined=$(length(store.combined)) ",
                           "partials=$(store.npartials) target=$target")
    end
    return store
end

"Square-root step for one dependency; returns a nontrivial divisor of `n` or `nothing`."
function try_dependency(rels::Vector{Relation}, dep::BitVector, fb::FactorBase, kn::BigInt, n::BigInt)
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
    for g in (gcd(X - Y, n), gcd(X + Y, n))
        1 < g < n && return g
    end
    return nothing
end

function dependencies_to_factor(store::RelationStore, fb::FactorBase, kn::BigInt, n::BigInt, rng)
    rels = vcat(store.fulls, store.combined)
    nrows = length(fb.primes) + 1
    cols = [parity_rows(r.fac) for r in rels]
    keep, nonempty = prune_singletons(cols, nrows)
    length(keep) - nonempty < 32 && return nothing
    kept_rels = rels[keep]
    kept_cols = cols[keep]
    for _ in 1:5
        for dep in block_lanczos(kept_cols, nrows; rng = rng)
            d = try_dependency(kept_rels, dep, fb, kn, n)
            d !== nothing && return d
        end
    end
    return nothing
end

"""
    siqs_split(n; rng, verbose) -> BigInt

Return a nontrivial divisor of the composite `n` using the self-initializing quadratic sieve.
`n` must have at least 40 decimal digits, no prime factor below 1000 and not be a perfect power.
"""
function siqs_split(n::BigInt; rng = Random.default_rng(), verbose::Bool = false)
    P = params_for(ndigits(n))
    k = knuth_schroeppel(n)
    kn = k * n
    fb = build_factorbase(kn, k, P.nprimes)
    st = SIQSState(kn, fb, P)
    store = RelationStore(kn)
    target = P.nprimes + 64
    try
        while true
            collect_relations!(store, st, target, rng, verbose)
            d = dependencies_to_factor(store, fb, kn, n, rng)
            d !== nothing && return d
            target += 64
        end
    catch e
        e isa FoundFactor || rethrow()
        g = gcd(e.d, n)
        1 < g < n && return g
        rethrow()
    end
end

end # module
