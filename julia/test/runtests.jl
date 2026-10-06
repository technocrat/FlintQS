using Test
@testset "FlintQS" begin
    for f in ("test_util.jl", "test_modarith.jl")
        include(f)
    end
end
