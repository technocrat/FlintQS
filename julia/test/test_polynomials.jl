using Test, Random
using FlintQS
const F = FlintQS
include("test_polynomials_helpers.jl")

@testset "polynomials" begin
    st, rng = make_state(40)
    @test st.s == ndigits(st.kn, base = 2) ÷ 28 + 1
    F.choose_a!(st, rng)
    @test length(unique(st.aidx)) == st.s
    @test all(>(st.params.firstprime), st.aidx)
    F.init_a!(st)
    A = prod(BigInt(st.fb.primes[i]) for i in st.aidx)
    @test st.A == A
    @test 0.5 < st.A / st.target < 2.0

    function check_poly(st)
        @test mod(st.B^2 - st.kn, st.A) == 0
        @test st.C == (st.B^2 - st.kn) ÷ st.A
        for ip in rand(Xoshiro(7), 1:length(st.fb.primes), 60)
            p = st.fb.primes[ip]
            Q(x) = st.A * x^2 + 2 * st.B * x + st.C
            if st.isa[ip]
                @test st.r1[ip] == st.r2[ip]
                @test mod(Q(BigInt(st.r1[ip])), p) == 0
            elseif p > 2
                @test mod(Q(BigInt(st.r1[ip])), p) == 0
                @test mod(Q(BigInt(st.r2[ip])), p) == 0
            end
        end
    end
    check_poly(st)
    seen = Set([st.B])
    for i in 1:(2^(st.s - 1) - 1)
        F.next_poly!(st, i)
        check_poly(st)
        push!(seen, st.B)
    end
    @test length(seen) == 2^(st.s - 1)      # every Gray-code step gives a new B
end
