# Benchmarks

Measured on Apple M1 Max with julia version 1.12.7. Each row is one `flintqs` call on a random balanced semiprime of the given size, after a warm-up call; every result is checked by multiplying the factors back. Peak RSS includes the Julia runtime (about 400 MB).

| digits | seconds | peak RSS (MB) |
|:--:|:--:|:--:|
| 40 | 0.07 | 417.0 |
| 45 | 0.21 | 435.0 |
| 50 | 0.57 | 552.0 |
| 55 | 1.13 | 746.0 |
| 60 | 2.98 | 917.0 |

Reproduce with `julia --project=julia julia/bench/bench.jl`.
