using Test
@testset "FlintQS" begin
    for f in ("test_util.jl", "test_modarith.jl", "test_params_fb.jl")
        include(f)
    end
end
