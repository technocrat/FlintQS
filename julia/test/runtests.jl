using Test
@testset "FlintQS" begin
    for f in ("test_util.jl",)
        include(f)
    end
end
