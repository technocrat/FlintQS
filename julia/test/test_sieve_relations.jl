using Test, Random
using FlintQS
const F = FlintQS
include("test_polynomials_helpers.jl")

@testset "sieve" begin
    st, rng = make_state(40)
    F.choose_a!(st, rng); F.init_a!(st)
    M = st.params.M
    sv = zeros(UInt8, 2M)
    F.sieve!(sv, st)
    for i0 in rand(rng, 0:2M-1, 20)
        x = i0 - M
        expect = 0
        for ip in st.params.firstprime+1:length(st.fb.primes)
            p = st.fb.primes[ip]
            xm = mod(x, p)
            if xm == st.r1[ip] || xm == st.r2[ip]
                expect += st.fb.sizes[ip]
            end
        end
        @test sv[i0+1] == expect
    end
end

@testset "relations" begin
    st, rng = make_state(40)
    store = F.RelationStore(st.kn)
    sv = zeros(UInt8, 2 * st.params.M)
    for _ in 1:3
        F.choose_a!(st, rng); F.init_a!(st)
        for i in 0:(2^(st.s - 1) - 1)
            i > 0 && F.next_poly!(st, i)
            F.sieve!(sv, st)
            F.scan!(store, st, sv)
        end
    end
    println("relations: full=", length(store.fulls), " combined=", length(store.combined),
            " partials=", store.npartials)
    @test F.nrelations(store) > 0
    rowval(r) = r == 1 ? BigInt(-1) : BigInt(st.fb.primes[r-1])
    for rel in vcat(store.fulls, store.combined)
        lhs = mod(rel.x^2, st.kn)
        rhs = mod(prod(rowval(r) for r in rel.fac; init = BigInt(1)), st.kn)
        @test lhs == rhs
    end
    @test F.parity_rows(Int32[3, 5, 3, 7, 7, 7]) == Int32[5, 7]
    cols = [Int32[1, 2], Int32[2, 3], Int32[3, 4], Int32[1, 5]]
    # cascade: row 5 appears once -> col 4 dropped -> row 1 now once -> col 1 dropped -> all pruned
    @test F.prune_singletons(cols, 5) == (Int[], 0)
    @test F.prune_singletons([Int32[1,2], Int32[1,2], Int32[3]], 3) == ([1, 2], 2)
end

@testset "large-prime merging" begin
    st, rng = make_state(40)
    P = st.params
    # relaxed threshold/error bits so cofactor candidates (partials) get through
    st.params = F.Params(P.digits, P.nprimes, P.M, P.largeprime, P.firstprime, P.errorbits + 8, P.threshold - 8)
    store = F.RelationStore(st.kn)
    sv = zeros(UInt8, 2 * st.params.M)
    for _ in 1:3
        F.choose_a!(st, rng); F.init_a!(st)
        for i in 0:(2^(st.s - 1) - 1)
            i > 0 && F.next_poly!(st, i)
            F.sieve!(sv, st)
            F.scan!(store, st, sv)
        end
    end
    println("relaxed: full=", length(store.fulls), " combined=", length(store.combined),
            " partials=", store.npartials)
    @test store.npartials > 0
    @test length(store.combined) > 0
    rowval(r) = r == 1 ? BigInt(-1) : BigInt(st.fb.primes[r-1])
    for rel in store.combined
        @test mod(rel.x^2, st.kn) == mod(prod(rowval(r) for r in rel.fac; init = BigInt(1)), st.kn)
    end
end
