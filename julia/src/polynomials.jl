# SPDX-License-Identifier: GPL-2.0-or-later

mutable struct SIQSState
    kn::BigInt
    fb::FactorBase
    params::Params
    s::Int
    lo::Int                 # 1-based index of first candidate A-prime
    span::Int               # number of candidate A-primes
    target::BigInt          # ideal A ≈ sqrt(2 kn) / M
    aidx::Vector{Int}
    A::BigInt
    B::BigInt
    C::BigInt
    bterms::Vector{BigInt}
    ainv2b::Matrix{Int}     # s × nprimes
    r1::Vector{Int}
    r2::Vector{Int}
    isa::BitVector
end

function SIQSState(kn::BigInt, fb::FactorBase, params::Params)
    np = length(fb.primes)
    s = ndigits(kn, base = 2) ÷ 28 + 1
    target = isqrt(2 * kn) ÷ params.M
    # Port of the C++ window computation (done with 0-based indices there).
    t = integer_root(target, s)
    fact = searchsortedfirst(fb.primes, Int(t)) - 1
    span = np ÷ s ÷ s ÷ 2
    lo = fact - span ÷ 2
    while (fact * fact) ÷ lo - lo < span
        lo -= 1
    end
    lo = max(lo, params.firstprime + 1)
    lo + span + s <= np || error("factor base too small for $s A-primes (lo=$lo, span=$span)")
    return SIQSState(kn, fb, params, s, lo + 1, span, target, zeros(Int, s), BigInt(0),
                     BigInt(0), BigInt(0), BigInt[], zeros(Int, s, np), zeros(Int, np),
                     zeros(Int, np), falses(np))
end

"Pick `s` distinct factor-base indices whose product is within a factor 2 of `st.target`."
function choose_a!(st::SIQSState, rng)
    fb = st.fb
    np = length(fb.primes)
    hi = st.lo + st.span - 1
    while true
        idx = Int[]
        A = BigInt(1)
        while length(idx) < st.s - 1
            i = rand(rng, st.lo:hi)
            (i in idx || fb.sqrts[i] == 0) && continue    # sqrt 0 means p | kn: B would vanish mod p
            push!(idx, i)
            A *= fb.primes[i]
        end
        want = Int(clamp(cld(st.target, A), 2, typemax(Int) ÷ 2))
        j = searchsortedfirst(fb.primes, want)
        for cand in (j, j - 1, j + 1, j - 2, j + 2)
            (st.params.firstprime < cand <= np) || continue
            (cand in idx || fb.sqrts[cand] == 0) && continue
            Afull = A * fb.primes[cand]
            if 0.5 < Afull / st.target < 2.0
                st.aidx = sort!(vcat(idx, cand))
                return st
            end
        end
    end
end

function a_root!(st::SIQSState, ip::Int)
    p = st.fb.primes[ip]
    Cp = Int(mod(st.C, p))
    Bp = Int(mod(st.B, p))
    r = mod(-Cp * invmod(mod(2 * Bp, p), p), p)
    st.r1[ip] = r
    st.r2[ip] = r
    return nothing
end

"Set up A, the B-terms, the first B and all roots for the A chosen in `st.aidx`."
function init_a!(st::SIQSState)
    fb, kn = st.fb, st.kn
    np = length(fb.primes)
    A = prod(BigInt(fb.primes[i]) for i in st.aidx)
    bterms = BigInt[]
    for i in st.aidx
        q = fb.primes[i]
        Aq = A ÷ q
        g = Int(mod(Aq, q))
        γ = mod(fb.sqrts[i] * invmod(g, q), q)
        γ = min(γ, q - γ)
        push!(bterms, Aq * γ)
    end
    st.A = A
    st.bterms = bterms
    st.B = sum(bterms)
    st.C = (st.B^2 - kn) ÷ A
    fill!(st.isa, false)
    for i in st.aidx
        st.isa[i] = true
    end
    for ip in 1:np
        p = fb.primes[ip]
        if st.isa[ip]
            a_root!(st, ip)
            continue
        end
        Ainv = p == 2 ? 1 : invmod(Int(mod(A, p)), p)
        Bp = Int(mod(st.B, p))
        for j in 1:st.s
            st.ainv2b[j, ip] = mod(2 * Int(mod(bterms[j], p)) * Ainv, p)
        end
        t = fb.sqrts[ip]
        st.r1[ip] = mod((t - Bp) * Ainv, p)
        st.r2[ip] = mod((-t - Bp) * Ainv, p)
    end
    return st
end

"Gray-code step `i` (1 ≤ i < 2^(s-1)): B ← B ± 2·bterm_j, roots shift by ∓ainv2b."
function next_poly!(st::SIQSState, i::Int)
    j = trailing_zeros(i) + 1
    sgn = if ((i >> (j - 1)) & 2) != 0
        st.B += 2 * st.bterms[j]
        -1
    else
        st.B -= 2 * st.bterms[j]
        1
    end
    st.C = (st.B^2 - st.kn) ÷ st.A
    fb = st.fb
    for ip in eachindex(fb.primes)
        if st.isa[ip]
            a_root!(st, ip)
        else
            p = fb.primes[ip]
            d = sgn * st.ainv2b[j, ip]
            st.r1[ip] = mod(st.r1[ip] + d, p)
            st.r2[ip] = mod(st.r2[ip] + d, p)
        end
    end
    return st
end
