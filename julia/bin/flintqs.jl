#!/usr/bin/env julia
# Usage: julia --project=julia julia/bin/flintqs.jl <n> [-v]
using FlintQS
isempty(ARGS) && (println(stderr, "usage: flintqs.jl <n> [-v]"); exit(1))
n = parse(BigInt, ARGS[1])
t = @elapsed fs = flintqs(n; verbose = "-v" in ARGS)
println(join(fs, " * "), "   (", round(t, digits = 3), " s)")
