using Test
@testset "FlintQS" begin
    for f in ("test_util.jl", "test_modarith.jl", "test_params_fb.jl", "test_polynomials.jl", "test_sieve_relations.jl", "test_lanczos.jl")
        include(f)
    end
end
