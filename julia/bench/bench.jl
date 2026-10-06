# Usage: julia --project=julia julia/bench/bench.jl [digits...]   (default: 40 45 50 55 60)
using FlintQS, Random

function semiprime(rng, digits)
    h = digits ÷ 2
    p = FlintQS.next_prime(BigInt(10)^(h - 1) + rand(rng, 1:10^6))
    q = FlintQS.next_prime(BigInt(10)^(digits - h - 1) + rand(rng, 1:10^6))
    return p * q
end

peak_mb() = Sys.maxrss() / 2^20

digits = isempty(ARGS) ? [40, 45, 50, 55, 60] : parse.(Int, ARGS)
flintqs(semiprime(Xoshiro(0), 40); rng = Xoshiro(0))   # warm-up (compilation)
println("| digits | seconds | peak RSS (MB) |\n|:--:|:--:|:--:|")
for d in digits
    n = semiprime(Xoshiro(d), d)
    t = @elapsed fs = flintqs(n; rng = Xoshiro(d))
    @assert prod(fs) == n
    println("| $d | $(round(t, digits = 2)) | $(round(peak_mb(), digits = 0)) |")
end
