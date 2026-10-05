# FlintQS Julia Port — Design

## Goal
A Julia port of William Hart's FlintQS (multiple-polynomial quadratic sieve with
large-prime relations and Block Lanczos), living in `julia/` inside this fork.
The C++ in `src/` is left untouched so upstream merges stay clean.

## Constraints
- **Clean room:** no code from the author's separate Julia SIQS repo. All helpers are
  written fresh from the FlintQS C++ and standard algorithms.
- **Self-contained:** depends only on Julia stdlib (`Random`, `Test`, `Printf`) and GMP via `BigInt`.
- **In-memory relations:** the C++ file spooling is not ported. This also avoids the file/buffer
  handling behind CVE-2023-29465.
- **Out of scope:** parallelism, Pari integration, factor-tree output.

## Success criteria
- `flintqs(n::BigInt)` returns a nontrivial factorization of composite `n`, verified by
  multiplying the factors back.
- Correct on random semiprimes of 20, 30, 40 and 50 digits (60 is a stretch target).
- Benchmark script reports wall time and peak memory per digit size.

## Layout (`julia/`)
| File | Ports | Purpose |
|---|---|---|
| `src/util.jl` | (new) | Miller-Rabin primality, perfect-power test, integer sqrt, small primes sieve, Pollard rho for small inputs |
| `src/modarith.jl` | `ModuloArith`, `TonelliShanks` | modular inverse, Jacobi symbol, square roots mod p |
| `src/params.jl` | tuning tables in `QS.cpp` | factor-base size, sieve size, large-prime cutoff by digit count |
| `src/factorbase.jl` | `knuthSchroeppel`, `computeFactorBase`, `computeSizes` | multiplier and factor base |
| `src/polynomials.jl` | polynomial selection in `mainRoutine` | choose A from factor-base primes, derive B, C |
| `src/sieve.jl` | `sieveInterval`, `sieve2` | sieve over `Vector{UInt8}` |
| `src/relations.jl` | `evaluateSieve`, `lprels` | trial-divide candidates, full relations, large-prime merging via cycle graph |
| `src/lanczos.jl` | `lanczos`, `F2matrix` | Block Lanczos (64-bit blocks) over sparse GF(2) matrix |
| `src/FlintQS.jl` | `main` | module and `flintqs(n)` driver |
| `test/` | | per-file unit tests plus end-to-end factoring tests |
| `bench/bench.jl` | | timing and memory benchmark |

## Data flow
`n` → prechecks (prime, perfect power, small factors, Pollard rho for small `n`) →
multiplier and factor base → loop {pick polynomial, sieve, evaluate candidates, record
relations}. The loop ends when relations exceed factor-base size plus a margin. Then
build the sparse matrix, run Block Lanczos for null vectors, and try gcds. If every
null vector fails, sieve further and retry (a gap listed in upstream `QStodo`).

## Error handling
Invalid input (n < 2, prime) returns early with a clear result. Internal failures
throw Julia exceptions rather than aborting. Lanczos failure falls back to more sieving.

## Testing
- Unit: modular inverse, square roots, Jacobi, Miller-Rabin vs. known primes and
  pseudoprimes, Lanczos on small matrices with known null spaces, relations satisfy
  x² ≡ Q(x) (mod n).
- End-to-end: random semiprimes at 20/30/40/50 digits, factors multiply back to `n`.

## Build order
util → modarith → params/factorbase → polynomials → sieve → relations → lanczos → driver → bench.
Each stage is tested before the next begins.
