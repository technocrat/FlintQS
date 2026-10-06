using Test, Random
using FlintQS
const F = FlintQS

function random_cols(rng, N, nrows, w)
    [sort!(unique(Int32.(rand(rng, 1:nrows, w)))) for _ in 1:N]
end

function mulB_ref(cols, nrows, v::BitVector)
    out = falses(nrows)
    for j in eachindex(cols)
        v[j] || continue
        for r in cols[j]
            out[r] = !out[r]
        end
    end
    out
end

@testset "lanczos 64x64 helpers" begin
    rng = Xoshiro(3)
    X = rand(rng, UInt64, 64); Y = rand(rng, UInt64, 64)
    bits(M) = [((M[r] >> c) & 1) == 1 for r in 1:length(M), c in 0:63]
    @test bits(F.mulmat(X, Y)) == ((Int.(bits(X)) * Int.(bits(Y))) .% 2 .== 1)
    @test F.mulmat(X, F.identity64()) == X
    V = rand(rng, UInt64, 200); W = rand(rng, UInt64, 200)
    @test bits(F.tmul(V, W)) == ((Int.(bits(V))' * Int.(bits(W))) .% 2 .== 1)
    M = F.identity64()
    for _ in 1:200
        i, j = rand(rng, 1:64, 2)
        i != j && (M[i] ⊻= M[j])
    end
    @test F.mulmat(M, F.inv64(M)) == F.identity64()
    G = rand(rng, UInt64, 64)
    Gt = vcat(G[1:40], zeros(UInt64, 24))
    T = F.tmul(Gt, Gt)                      # symmetric, rank ≤ 40
    mask, winv = F.select_cols(T, typemax(UInt64))
    @test count_ones(mask) == 40
    WT = F.mulmat(winv, T)
    for c in 0:63
        (mask >> c) & 1 == 1 && @test (WT[c+1] & mask) == (UInt64(1) << c)
    end
end

@testset "block_lanczos" begin
    rng = Xoshiro(11)
    for (N, nrows, w) in ((164, 100, 5), (400, 300, 6), (1100, 1000, 8))
        cols = random_cols(rng, N, nrows, w)
        deps = BitVector[]
        for attempt in 1:6
            deps = F.block_lanczos(cols, nrows; rng = rng)
            isempty(deps) || break
        end
        @test !isempty(deps)
        for v in deps
            @test length(v) == N && any(v)
            @test !any(mulB_ref(cols, nrows, v))
        end
    end
end

@testset "block_lanczos reliability" begin
    # Breakdowns are legitimate but must be rare, otherwise the driver wastes sieving time.
    rng = Xoshiro(21)
    N, nrows, w = 3200, 3100, 20
    ok = 0
    for t in 1:10
        cols = random_cols(rng, N, nrows, w)
        deps = F.block_lanczos(cols, nrows; rng = rng)
        if !isempty(deps)
            ok += 1
            @test all(v -> !any(mulB_ref(cols, nrows, v)), deps)
        end
    end
    @test ok >= 8
end
