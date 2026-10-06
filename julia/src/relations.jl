# SPDX-License-Identifier: GPL-2.0-or-later

struct Relation
    x::BigInt
    fac::Vector{Int32}
end

mutable struct RelationStore
    kn::BigInt
    fulls::Vector{Relation}
    combined::Vector{Relation}
    hub::Dict{Int,Relation}
    seen::Set{BigInt}
    npartials::Int
    found::Vector{BigInt}       # divisors of kn revealed by large primes (gcd(q, kn) > 1)
end
RelationStore(kn::BigInt) = RelationStore(kn, Relation[], Relation[], Dict{Int,Relation}(), Set{BigInt}(), 0, BigInt[])

nrelations(s::RelationStore) = length(s.fulls) + length(s.combined)

"""
Divide `p` out of `res` in place as often as possible and return the exponent. `qb`, `rb`, `pb`
are caller-owned scratch BigInts so a failed division allocates nothing.
"""
function strip_prime!(res::BigInt, p::Int, qb::BigInt, rb::BigInt, pb::BigInt)
    Base.GMP.MPZ.set_ui!(pb, p)
    e = 0
    while true
        Base.GMP.MPZ.tdiv_qr!(qb, rb, res, pb)
        iszero(rb) || return e
        Base.GMP.MPZ.set!(res, qb)
        e += 1
    end
end

"""
Fully trial-divide Q(x) at sieve index `i0` (0-based). Returns `nothing`, or `(Relation, q)` where
`q == 1` for a full relation and `q > 1` is the single large prime of a partial relation.
"""
function try_candidate(st::SIQSState, sv::Vector{UInt8}, i0::Int)
    P = st.params
    fb = st.fb
    x = i0 - P.M
    t = st.A * x + st.B                      # A x + B
    Q = (t + st.B) * x + st.C                # A x^2 + 2 B x + C;  t^2 ≡ A Q (mod kn)
    res = abs(Q)
    res == 0 && return nothing
    fac = Int32[]
    Q < 0 && push!(fac, Int32(1))
    target = ndigits(res, base = 2) - P.errorbits
    qb, rb, pb = BigInt(), BigInt(), BigInt()

    extra = 0
    for ip in 1:P.firstprime
        e = strip_prime!(res, fb.primes[ip], qb, rb, pb)
        if e > 0
            extra += fb.sizes[ip]
            for _ in 1:e
                push!(fac, Int32(ip + 1))
            end
        end
    end
    Int(sv[i0+1]) + extra < target && return nothing

    for ip in P.firstprime+1:length(fb.primes)
        p = fb.primes[ip]
        xm = mod(x, p)
        (xm == st.r1[ip] || xm == st.r2[ip]) || continue
        e = strip_prime!(res, p, qb, rb, pb)
        for _ in 1:e
            push!(fac, Int32(ip + 1))
        end
    end
    for ip in st.aidx
        push!(fac, Int32(ip + 1))
    end
    rel = Relation(mod(t, st.kn), fac)
    res == 1 && return rel, 1
    (res <= P.largeprime && res < BigInt(fb.primes[end])^2) || return nothing
    return rel, Int(res)
end

function combine(r1::Relation, r2::Relation, q::Int, kn::BigInt)
    x = mod(r1.x * r2.x * invmod(BigInt(q), kn), kn)
    return Relation(x, vcat(r1.fac, r2.fac))
end

function add_candidate!(store::RelationStore, rel::Relation, q::Int)
    rel.x in store.seen && return
    push!(store.seen, rel.x)
    if q == 1
        push!(store.fulls, rel)
    else
        store.npartials += 1
        if haskey(store.hub, q)
            g = gcd(BigInt(q), store.kn)
            if g != 1                        # q | kn: not invertible, but g is a divisor of kn
                push!(store.found, g)
                return nothing
            end
            push!(store.combined, combine(rel, store.hub[q], q, store.kn))
        else
            store.hub[q] = rel
        end
    end
    return nothing
end

function scan!(store::RelationStore, st::SIQSState, sv::Vector{UInt8})
    thr = st.params.threshold
    @inbounds for j in eachindex(sv)
        sv[j] >= thr || continue
        out = try_candidate(st, sv, j - 1)
        out === nothing && continue
        add_candidate!(store, out[1], out[2])
    end
    return store
end

"Sorted matrix rows that occur an odd number of times in `fac`."
function parity_rows(fac::AbstractVector{<:Integer})
    s = sort(fac)
    out = Int32[]
    i = 1
    while i <= length(s)
        j = i
        while j < length(s) && s[j+1] == s[i]
            j += 1
        end
        isodd(j - i + 1) && push!(out, Int32(s[i]))
        i = j + 1
    end
    return out
end

"Drop columns that contain a row occurring in no other column (repeatedly)."
function prune_singletons(cols::Vector{Vector{Int32}}, nrows::Int)
    keep = trues(length(cols))
    cnt = zeros(Int, nrows)
    for c in cols, r in c
        cnt[r] += 1
    end
    changed = true
    while changed
        changed = false
        for j in eachindex(cols)
            keep[j] || continue
            if any(r -> cnt[r] == 1, cols[j])
                keep[j] = false
                changed = true
                for r in cols[j]
                    cnt[r] -= 1
                end
            end
        end
    end
    return findall(keep), count(>(0), cnt)
end
