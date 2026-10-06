# Benchmarks

Measured on Apple M1 Max with julia version 1.12.7. Each row is one `flintqs` call on a random balanced semiprime of the given size, after a warm-up call; every result is checked by multiplying the factors back. Peak RSS is the process-wide running maximum (it includes the Julia runtime, about 400 MB), so each row is at least the rows above it.

| digits | seconds | peak RSS (MB) |
|:--:|:--:|:--:|
| 40 | 0.06 | 382.0 |
| 45 | 0.17 | 411.0 |
| 50 | 0.6 | 524.0 |
| 55 | 0.95 | 600.0 |
| 60 | 2.73 | 732.0 |

Reproduce with `julia --project=julia julia/bench/bench.jl`.
