# FlintQS Julia Port Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Julia port of FlintQS's self-initializing quadratic sieve (SIQS) with single-large-prime merging and Block Lanczos, in `julia/`, with tests, a benchmark, and MaterialDocs documentation.

**Architecture:** One Julia package `FlintQS` (`julia/`), eight small source files that mirror the C++ units (utilities, modular arithmetic, parameter tables, factor base, polynomial generation, sieve, relation handling, Block Lanczos) plus a driver. Relations live in memory. The C++ in `src/` is untouched.

**Tech Stack:** Julia 1.12 (stdlib only: `Random`, `Test`, `Printf`), GMP via `BigInt`; Documenter.jl + MaterialDocs.jl for docs.

**Spec:** `docs/superpowers/specs/2026-10-05-flintqs-julia-design.md`

## Global Constraints

- Package root is `julia/`; package name `FlintQS`, uuid `f2b7871c-ac93-4404-9113-4d8c0e711c3b`.
- Runtime dependencies: Julia stdlib only. No code from the author's separate Julia SIQS repo (clean room w.r.t. that repo).
- This is a port of GPL-2.0-or-later code (see `COPYING`): every new source file carries the header `# SPDX-License-Identifier: GPL-2.0-or-later` and `julia/` is distributed under the same terms.
- FlintQS tables cover 40-91 decimal digits (`MINDIG 40`). Inputs of 19-39 digits are handled by multiplying by a prime so the work number has >= 40 digits ("lift"); inputs of <= 18 digits use Pollard rho.
- Factor-base primes must be < 2^30 so `Int` products cannot overflow.
- Relations are in memory only; no temp files.
- Sign of Q(x) is tracked as matrix row 1 (the C++ ignores sign; this is a deliberate deviation).
- Test command (from repo root): `julia --project=julia julia/test/<file>.jl`; full suite: `julia --project=julia -e 'using Pkg; Pkg.test()'`.
- Every commit message ends with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.

## Review Focus

Failure modes the spec implies but no happy-path test exercises (each is pinned in Task 8):
1. `n` prime, 0, 1, or negative: return `[n]` for primes, `BigInt[]` for 1, throw `ArgumentError` for `n < 1`.
2. `n` a perfect power (e.g. `p^2`, `p^3` with p ~ 25 digits): must return p repeated, not loop in SIQS (SIQS cannot split prime powers).
3. `n` with many small factors times one large prime (e.g. `2^5 * 3^3 * P`): small factors stripped first, remainder prime.
4. Digit-count boundaries: 39 digits (lift path) and 40 digits (direct) both factor correctly.
5. Very unbalanced semiprime (7-digit prime times 45-digit prime): result still correct (large factor found by sieve, small one is just another gcd).

## File Structure

```
julia/
  Project.toml
  .gitignore                 # Manifest.toml
  src/FlintQS.jl             # module, includes, exports, flintqs(), siqs_split()
  src/util.jl                # primes_upto, is_probable_prime, integer_root, perfect_power, pollard_rho, next_prime
  src/modarith.jl            # jacobi, sqrt_mod_prime
  src/params.jl              # Params + tables + params_for
  src/factorbase.jl          # FactorBase, knuth_schroeppel, build_factorbase
  src/polynomials.jl         # SIQSState, a_window, choose_a!, init_a!, next_poly!
  src/sieve.jl               # sieve!
  src/relations.jl           # Relation, RelationStore, try_candidate, scan!, combine, parity_rows, prune_singletons
  src/lanczos.jl             # 64x64 GF(2) helpers, block_lanczos
  test/runtests.jl           # includes every test_*.jl
  test/test_util.jl  test_modarith.jl  test_params_fb.jl  test_polynomials.jl
  test/test_sieve_relations.jl  test_lanczos.jl  test_flintqs.jl
  bin/flintqs.jl             # CLI
  bench/bench.jl
  docs/Project.toml  docs/make.jl  docs/src/{index,algorithm,api,benchmarks}.md
```

Every test file begins with:
```julia
using Test, Random
using FlintQS
const F = FlintQS
```
`src/FlintQS.jl` grows one `include` per task (listed in each task).

---

### Task 1: Package skeleton and utilities

**Files:**
- Create: `julia/Project.toml`, `julia/.gitignore`, `julia/src/FlintQS.jl`, `julia/src/util.jl`, `julia/test/runtests.jl`, `julia/test/test_util.jl`

**Interfaces:**
- Produces: `primes_upto(limit)::Vector{Int}`, `SMALL_PRIMES::Vector{Int}` (primes <= 1000), `is_probable_prime(n::BigInt)::Bool`, `integer_root(n::BigInt, k::Int)::BigInt`, `perfect_power(n::BigInt)::Union{Nothing,Tuple{BigInt,Int}}`, `pollard_rho(n::BigInt; maxiter, rng)::Union{Nothing,BigInt}`, `next_prime(n::BigInt)::BigInt`.

- [ ] **Step 1: Create skeleton files**

`julia/Project.toml`:
```toml
name = "FlintQS"
uuid = "f2b7871c-ac93-4404-9113-4d8c0e711c3b"
version = "0.1.0"

[deps]
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Printf = "de0858da-6303-5e9b-ae3b-b3bb9a6d47d5"

[compat]
julia = "1.10"

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]
```
`julia/.gitignore`:
```
Manifest.toml
docs/build/
```
`julia/src/FlintQS.jl`:
```julia
# SPDX-License-Identifier: GPL-2.0-or-later
module FlintQS

using Random

include("util.jl")

end # module
```
`julia/test/runtests.jl`:
```julia
using Test
@testset "FlintQS" begin
    for f in ("test_util.jl",)
        include(f)
    end
end
```

- [ ] **Step 2: Write the failing test** — `julia/test/test_util.jl`

```julia
using Test, Random
using FlintQS
const F = FlintQS

@testset "util" begin
    @test F.primes_upto(30) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    @test F.primes_upto(1) == Int[]
    @test length(F.SMALL_PRIMES) == 168

    @test F.is_probable_prime(big"2") && F.is_probable_prime(big"97")
    @test !F.is_probable_prime(big"1") && !F.is_probable_prime(big"0")
    @test F.is_probable_prime(big"170141183460469231731687303715884105727")   # 2^127-1
    @test !F.is_probable_prime(big"561")           # Carmichael
    @test !F.is_probable_prime(big"3215031751")    # strong pseudoprime to bases 2,3,5,7
    @test !F.is_probable_prime(big"1000000007" * big"1000000009")

    @test F.integer_root(big"1000000", 2) == 1000
    @test F.integer_root(big"999999", 2) == 999
    @test F.integer_root(big"3"^40, 5) == big"3"^8
    @test F.perfect_power(big"1024") == (2, 10)
    @test F.perfect_power(big"7"^3) == (7, 3)
    @test F.perfect_power(big"1000003"^2) == (1000003, 2)
    @test F.perfect_power(big"1000003" * big"1000033") === nothing
    @test F.perfect_power(big"3") === nothing

    n = big"1000003" * big"1000033"
    d = F.pollard_rho(n; rng = Xoshiro(1))
    @test d !== nothing && 1 < d < n && n % d == 0
    @test F.pollard_rho(big"1000003"; maxiter = 10_000, rng = Xoshiro(1)) === nothing

    @test F.next_prime(big"100") == 101
    @test F.next_prime(big"101") == 101
end
```

- [ ] **Step 3: Run test to verify it fails**

Run: `julia --project=julia julia/test/test_util.jl`
Expected: FAIL (`UndefVarError: primes_upto` or similar).

- [ ] **Step 4: Implement** — `julia/src/util.jl`

```julia
# SPDX-License-Identifier: GPL-2.0-or-later

"Primes <= `limit` (sieve of Eratosthenes)."
function primes_upto(limit::Integer)
    limit < 2 && return Int[]
    isp = trues(limit)
    isp[1] = false
    for i in 2:isqrt(limit)
        isp[i] || continue
        for j in i*i:i:limit
            isp[j] = false
        end
    end
    return findall(isp)
end

const SMALL_PRIMES = primes_upto(1000)

"Miller-Rabin with the first 24 primes as bases (deterministic below 3.3e24, probabilistic above)."
function is_probable_prime(n::BigInt)
    n < 2 && return false
    for p in SMALL_PRIMES
        n == p && return true
        n % p == 0 && return false
    end
    d = n - 1
    r = 0
    while iseven(d)
        d >>= 1
        r += 1
    end
    for a in SMALL_PRIMES[1:24]
        x = powermod(BigInt(a), d, n)
        (x == 1 || x == n - 1) && continue
        composite = true
        for _ in 1:r-1
            x = mod(x * x, n)
            if x == n - 1
                composite = false
                break
            end
        end
        composite && return false
    end
    return true
end

"Smallest prime >= n."
function next_prime(n::BigInt)
    n <= 2 && return big"2"
    m = iseven(n) ? n + 1 : n
    while !is_probable_prime(m)
        m += 2
    end
    return m
end

"floor(n^(1/k)) by Newton iteration from above."
function integer_root(n::BigInt, k::Int)
    (n < 2 || k == 1) && return n
    x = BigInt(1) << cld(ndigits(n, base = 2), k)
    while true
        y = ((k - 1) * x + n ÷ x^(k - 1)) ÷ k
        y >= x && return x
        x = y
    end
end

"If `n = b^k` (k >= 2) return `(b, k)` with the largest such k, else `nothing`."
function perfect_power(n::BigInt)
    n < 4 && return nothing
    for k in ndigits(n, base = 2):-1:2
        b = integer_root(n, k)
        b > 1 && b^k == n && return (b, k)
    end
    return nothing
end

"Brent's variant of Pollard rho. Returns a nontrivial factor or `nothing` after ~`maxiter` steps."
function pollard_rho(n::BigInt; maxiter::Int = 5_000_000, rng = Random.default_rng())
    iseven(n) && return BigInt(2)
    iters = 0
    while iters < maxiter
        y = BigInt(rand(rng, 1:1_000_000))
        c = BigInt(rand(rng, 1:1_000_000))
        m = 128
        g = BigInt(1)
        r = 1
        q = BigInt(1)
        x = y
        ys = y
        while g == 1 && iters < maxiter
            x = y
            for _ in 1:r
                y = mod(y * y + c, n)
            end
            k = 0
            while k < r && g == 1
                ys = y
                for _ in 1:min(m, r - k)
                    y = mod(y * y + c, n)
                    q = mod(q * abs(x - y), n)
                end
                g = gcd(q, n)
                k += m
            end
            iters += r
            r *= 2
        end
        if g == n
            g = BigInt(1)
            while g == 1
                ys = mod(ys * ys + c, n)
                g = gcd(abs(x - ys), n)
            end
        end
        (g != 1 && g != n) && return g
    end
    return nothing
end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `julia --project=julia julia/test/test_util.jl`
Expected: all tests PASS. (The prime `1000003` check with `maxiter=10_000` must return `nothing`; if rho loops forever on a prime, fix the inner `while g == n` fallback to also respect `maxiter`.)

- [ ] **Step 6: Commit**

```bash
git add julia && git commit -m "Add Julia package skeleton and number-theory utilities

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Modular arithmetic

**Files:**
- Create: `julia/src/modarith.jl`, `julia/test/test_modarith.jl`
- Modify: `julia/src/FlintQS.jl` (add `include("modarith.jl")`), `julia/test/runtests.jl` (add `"test_modarith.jl"`)

**Interfaces:**
- Produces: `jacobi(a::Int, n::Int)::Int`, `jacobi(a::BigInt, n::Int)::Int` (n odd positive), `sqrt_mod_prime(a::Int, p::Int)::Int` (a a QR or 0 mod p; p prime < 2^30; returns r with r^2 ≡ a).

- [ ] **Step 1: Write the failing test** — `julia/test/test_modarith.jl`

```julia
using Test, Random
using FlintQS
const F = FlintQS

@testset "modarith" begin
    # jacobi vs Euler's criterion for prime moduli
    for p in (3, 5, 7, 11, 13, 101, 997)
        for a in 1:p-1
            e = powermod(a, (p - 1) ÷ 2, p)
            @test F.jacobi(a, p) == (e == 1 ? 1 : -1)
        end
        @test F.jacobi(0, p) == 0
    end
    @test F.jacobi(2, 15) == 1       # composite modulus
    @test F.jacobi(5, 9) == 1
    @test F.jacobi(big"123456789012345678901234567890", 1009) ==
          F.jacobi(Int(mod(big"123456789012345678901234567890", 1009)), 1009)
    @test_throws ArgumentError F.jacobi(3, 8)

    # sqrt_mod_prime: p ≡ 3 mod 4, p ≡ 1 mod 4 (incl. p ≡ 1 mod 8, 17, 97 with large 2-power)
    for p in (3, 7, 11, 13, 17, 41, 97, 193, 1009, 786433, 1000003)
        for a in unique(rand(Xoshiro(p), 0:p-1, 50))
            if a == 0 || F.jacobi(a, p) == 1
                r = F.sqrt_mod_prime(a, p)
                @test mod(r * r, p) == a
            end
        end
    end
    @test F.sqrt_mod_prime(0, 13) == 0
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=julia julia/test/test_modarith.jl` — Expected: FAIL (`jacobi` undefined).

- [ ] **Step 3: Implement** — `julia/src/modarith.jl`

```julia
# SPDX-License-Identifier: GPL-2.0-or-later

"Jacobi symbol (a/n) for odd positive n."
function jacobi(a::Int, n::Int)
    (n > 0 && isodd(n)) || throw(ArgumentError("jacobi: n must be odd and positive"))
    a = mod(a, n)
    result = 1
    while a != 0
        while iseven(a)
            a >>= 1
            r = n & 7
            (r == 3 || r == 5) && (result = -result)
        end
        a, n = n, a
        (a & 3 == 3 && n & 3 == 3) && (result = -result)
        a = mod(a, n)
    end
    return n == 1 ? result : 0
end

jacobi(a::BigInt, n::Int) = jacobi(Int(mod(a, n)), n)

"Square root of `a` modulo the odd prime `p` (Tonelli-Shanks). `a` must be 0 or a residue; p < 2^30."
function sqrt_mod_prime(a::Int, p::Int)
    a = mod(a, p)
    (a == 0 || p == 2) && return a
    p & 3 == 3 && return powermod(a, (p + 1) >> 2, p)
    q = p - 1
    s = 0
    while iseven(q)
        q >>= 1
        s += 1
    end
    z = 2
    while jacobi(z, p) != -1
        z += 1
    end
    m = s
    c = powermod(z, q, p)
    t = powermod(a, q, p)
    r = powermod(a, (q + 1) >> 1, p)
    while t != 1
        i = 0
        tt = t
        while tt != 1
            tt = mod(tt * tt, p)
            i += 1
        end
        b = c
        for _ in 1:(m - i - 1)
            b = mod(b * b, p)
        end
        m = i
        c = mod(b * b, p)
        t = mod(t * c, p)
        r = mod(r * b, p)
    end
    return r
end
```

- [ ] **Step 4: Run to verify it passes** — `julia --project=julia julia/test/test_modarith.jl` — Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add julia && git commit -m "Add Jacobi symbol and Tonelli-Shanks square roots

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Parameters and factor base

**Files:**
- Create: `julia/src/params.jl`, `julia/src/factorbase.jl`, `julia/test/test_params_fb.jl`
- Modify: `FlintQS.jl` (include `params.jl`, `factorbase.jl` after `modarith.jl`), `runtests.jl`

**Interfaces:**
- Consumes: `jacobi`, `sqrt_mod_prime`, `primes_upto`, `SMALL_PRIMES`.
- Produces:
  - `struct Params; digits::Int; nprimes::Int; M::Int; largeprime::Int; firstprime::Int; errorbits::Int; threshold::Int; end` — sieve covers `x ∈ [-M, M)`; `firstprime` = number of leading factor-base primes NOT sieved (trial-divided instead).
  - `params_for(digits::Int)::Params` (throws `ArgumentError` for digits < 40).
  - `struct FactorBase; primes::Vector{Int}; sizes::Vector{UInt8}; sqrts::Vector{Int}; multiplier::Int; end` (`primes[1] == 2`, `sqrts[i]` = √(kn) mod primes[i], 0 if p | kn).
  - `knuth_schroeppel(n::BigInt)::Int`, `build_factorbase(kn::BigInt, k::Int, nprimes::Int)::FactorBase`.

- [ ] **Step 1: Write the failing test** — `julia/test/test_params_fb.jl`

```julia
using Test, Random
using FlintQS
const F = FlintQS

@testset "params" begin
    p = F.params_for(40)
    @test (p.nprimes, p.M, p.firstprime, p.errorbits, p.threshold) == (1500, 32000, 8, 16, 66)
    @test F.params_for(60).nprimes == 3000
    @test F.params_for(91).nprimes == 80000
    @test F.params_for(91).threshold == 102
    @test F.params_for(92).nprimes == 64000
    @test_throws ArgumentError F.params_for(39)
    # table lengths are consistent (52 rows, digits 40..91)
    for d in 40:91
        q = F.params_for(d)
        @test q.largeprime > 1000 && q.threshold > q.errorbits
    end
end

@testset "factor base" begin
    n = big"1000003" * big"1000033" * big"999983"
    k = F.knuth_schroeppel(n)
    @test k in (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    kn = k * n
    fb = F.build_factorbase(kn, k, 200)
    @test length(fb.primes) == 200 && fb.primes[1] == 2
    @test issorted(fb.primes) && allunique(fb.primes)
    for (i, p) in enumerate(fb.primes)
        p == 2 && continue
        @test mod(fb.sqrts[i]^2, p) == mod(kn, p)
        @test F.is_probable_prime(BigInt(p))
    end
    @test fb.sizes[2] == round(Int, log2(fb.primes[2]) - 0.15)
    # Knuth-Schroeppel prefers a multiplier making kn a QR mod many small primes
    @test F.knuth_schroeppel(big"1") in (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
end
```

- [ ] **Step 2: Run to verify it fails** — `julia --project=julia julia/test/test_params_fb.jl` — Expected: FAIL.

- [ ] **Step 3: Implement `julia/src/params.jl`** (tables copied from `src/QS.cpp`, digits 40..91; the C++ value `8600000` for 49 digits is an apparent typo and is `860000` here)

```julia
# SPDX-License-Identifier: GPL-2.0-or-later

# Tuning tables from FlintQS QS.cpp, index i = digits - 39 (digits 40..91).
const LARGE_PRIMES = [
    250000, 300000, 370000, 440000, 510000, 580000, 650000, 720000, 790000, 860000,
    930000, 1000000, 1700000, 2400000, 3100000, 3800000, 4500000, 5200000, 5900000, 6600000,
    7300000, 8000000, 8900000, 10000000, 11300000, 12800000, 14500000, 16300000, 18100000, 20000000,
    22000000, 24000000, 27000000, 32000000, 39000000,
    53000000, 65000000, 75000000, 87000000, 100000000,
    114000000, 130000000, 150000000, 172000000, 195000000,
    220000000, 250000000, 300000000, 350000000, 400000000,
    450000000, 500000000]

const PRIMES_NO = [
    1500, 1500, 1600, 1700, 1750, 1800, 1900, 2000, 2050, 2100,
    2150, 2200, 2250, 2300, 2400, 2500, 2600, 2700, 2800, 2900,
    3000, 3150, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000,
    9500, 10000, 11500, 13000, 15000,
    17000, 24000, 27000, 30000, 37000,
    45000, 47000, 53000, 57000, 58000,
    59000, 60000, 64000, 68000, 72000,
    76000, 80000]

const FIRST_PRIMES = [
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    9, 8, 9, 9, 9, 9, 10, 10, 10, 10,
    10, 10, 11, 11, 12, 12, 13, 14, 15, 17,
    19, 21, 22, 22, 23,
    24, 25, 25, 26, 26,
    27, 27, 27, 27, 28,
    28, 28, 28, 29, 29,
    29, 29]

const ERROR_BITS = [
    16, 17, 17, 18, 18, 19, 19, 19, 20, 20,
    21, 21, 21, 22, 22, 22, 23, 23, 23, 24,
    24, 24, 25, 25, 25, 25, 26, 26, 26, 26,
    27, 27, 28, 28, 29,
    29, 30, 30, 30, 31,
    31, 31, 31, 32, 32,
    32, 32, 32, 33, 33,
    33, 33]

const THRESHOLDS = [
    66, 67, 67, 68, 68, 68, 69, 69, 69, 69,
    70, 70, 70, 71, 71, 71, 72, 72, 73, 73,
    74, 74, 75, 75, 76, 76, 77, 77, 78, 79,
    80, 81, 82, 83, 84,
    85, 86, 87, 88, 89,
    91, 92, 93, 93, 94,
    95, 96, 97, 98, 100,
    101, 102]

# Half the sieve length (sieve covers x in [-M, M)).
const SIEVE_HALF = vcat(fill(32000, 32), [64000, 64000, 64000],
                        [96000, 96000, 96000, 128000, 128000],
                        fill(160000, 5), fill(192000, 7))

struct Params
    digits::Int
    nprimes::Int
    M::Int
    largeprime::Int
    firstprime::Int
    errorbits::Int
    threshold::Int
end

function params_for(digits::Integer)
    digits >= 40 || throw(ArgumentError("SIQS parameters start at 40 digits (got $digits)"))
    if digits <= 91
        i = digits - 39
        return Params(digits, PRIMES_NO[i], SIEVE_HALF[i], LARGE_PRIMES[i],
                      FIRST_PRIMES[i], ERROR_BITS[i], THRESHOLDS[i])
    end
    # "all bets are off" branch of the C++ main()
    return Params(digits, 64000, 192000, 64000 * 10 * digits, 30,
                  digits ÷ 4 + 2, 43 + (7 * digits) ÷ 10)
end
```

- [ ] **Step 4: Implement `julia/src/factorbase.jl`**

```julia
# SPDX-License-Identifier: GPL-2.0-or-later

const MULTIPLIERS = (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)

struct FactorBase
    primes::Vector{Int}
    sizes::Vector{UInt8}
    sqrts::Vector{Int}
    multiplier::Int
end

"Knuth-Schroeppel multiplier: favour k for which k*n is a residue modulo many small primes."
function knuth_schroeppel(n::BigInt)
    best, bestscore = 1, -Inf
    for k in MULTIPLIERS
        kn = k * n
        score = -0.5 * log(k)
        r8 = Int(mod(kn, 8))
        score += r8 == 1 ? 2 * log(2) : r8 == 5 ? log(2) : (r8 == 3 || r8 == 7) ? 0.5 * log(2) : 0.0
        for p in SMALL_PRIMES
            p == 2 && continue
            if k % p == 0
                score += log(p) / p
            elseif jacobi(kn, p) == 1
                score += 2 * log(p) / p
            end
        end
        if score > bestscore
            best, bestscore = k, score
        end
    end
    return best
end

"First `nprimes` primes p (including 2) for which `kn` is a square mod p, with roots and bit sizes."
function build_factorbase(kn::BigInt, k::Int, nprimes::Int)
    limit = max(1000, ceil(Int, 2.5 * nprimes * (log(nprimes) + 2)))
    while true
        primes = Int[]
        sqrts = Int[]
        for p in primes_upto(limit)
            if p == 2
                push!(primes, 2)
                push!(sqrts, Int(mod(kn, 2)))
            else
                a = Int(mod(kn, p))
                (a == 0 || jacobi(a, p) == 1) || continue
                push!(primes, p)
                push!(sqrts, a == 0 ? 0 : sqrt_mod_prime(a, p))
            end
            length(primes) == nprimes && break
        end
        if length(primes) == nprimes
            sizes = UInt8[UInt8(floor(Int, log2(p) - 0.15 + 0.5)) for p in primes]
            return FactorBase(primes, sizes, sqrts, k)
        end
        limit *= 2
    end
end
```

- [ ] **Step 5: Wire up and run** — add the two includes and `"test_params_fb.jl"`; run `julia --project=julia julia/test/test_params_fb.jl` — Expected: PASS. (`SIEVE_HALF` has 32+3+5+5+7 = 52 entries; the `params` testset loop over 40:91 catches a wrong length via `BoundsError`.)

- [ ] **Step 6: Commit**

```bash
git add julia && git commit -m "Add parameter tables and factor base construction

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Polynomial generation (self-initialization)

**Files:**
- Create: `julia/src/polynomials.jl`, `julia/test/test_polynomials.jl`
- Modify: `FlintQS.jl` (include after `factorbase.jl`), `runtests.jl`

**Interfaces:**
- Consumes: `Params`, `FactorBase`, `integer_root`.
- Produces:
  - `mutable struct SIQSState` with fields `kn::BigInt, fb::FactorBase, params::Params, s::Int, lo::Int, span::Int, target::BigInt, aidx::Vector{Int}, A::BigInt, B::BigInt, C::BigInt, bterms::Vector{BigInt}, ainv2b::Matrix{Int}, r1::Vector{Int}, r2::Vector{Int}, isa::BitVector` and constructor `SIQSState(kn, fb, params)`.
  - `choose_a!(st, rng)`: fills `st.aidx` (1-based indices into `fb.primes`, distinct, all > `params.firstprime`) with `A ≈ st.target` (within factor 2).
  - `init_a!(st)`: given `st.aidx`, sets `A, bterms, B, C, ainv2b, r1, r2, isa` for the first polynomial of this A.
  - `next_poly!(st, i::Int)` for `i` in `1:2^(s-1)-1`: Gray-code step to the next B (updates `B, C, r1, r2`).
  - Invariants: `B^2 ≡ kn (mod A)`, `C = (B^2 - kn) ÷ A` exactly, and for each FB prime index `ip` with `p ∤ A`, `Q(r) ≡ 0 (mod p)` for `r ∈ {r1[ip], r2[ip]}` where `Q(x) = A x^2 + 2 B x + C`; for `isa[ip]`, `r2[ip] == r1[ip]` is the single root.

- [ ] **Step 1: Write the failing test** — `julia/test/test_polynomials.jl`

```julia
using Test, Random
using FlintQS
const F = FlintQS

function make_state(digits; seed = 1)
    rng = Xoshiro(seed)
    p = F.next_prime(BigInt(10)^(digits ÷ 2) + rand(rng, 1:10^6))
    q = F.next_prime(BigInt(10)^(digits - digits ÷ 2) + rand(rng, 1:10^6))
    n = p * q
    k = F.knuth_schroeppel(n)
    kn = k * n
    P = F.params_for(ndigits(n))
    fb = F.build_factorbase(kn, k, P.nprimes)
    return F.SIQSState(kn, fb, P), rng
end

@testset "polynomials" begin
    st, rng = make_state(40)
    @test st.s == ndigits(st.kn, base = 2) ÷ 28 + 1
    F.choose_a!(st, rng)
    @test length(unique(st.aidx)) == st.s
    @test all(>(st.params.firstprime), st.aidx)
    F.init_a!(st)
    A = prod(BigInt(st.fb.primes[i]) for i in st.aidx)
    @test st.A == A
    @test 0.5 < st.A / st.target < 2.0

    function check_poly(st)
        @test mod(st.B^2 - st.kn, st.A) == 0
        @test st.C == (st.B^2 - st.kn) ÷ st.A
        for ip in rand(Xoshiro(7), 1:length(st.fb.primes), 60)
            p = st.fb.primes[ip]
            Q(x) = st.A * x^2 + 2 * st.B * x + st.C
            if st.isa[ip]
                @test st.r1[ip] == st.r2[ip]
                @test mod(Q(BigInt(st.r1[ip])), p) == 0
            elseif p > 2
                @test mod(Q(BigInt(st.r1[ip])), p) == 0
                @test mod(Q(BigInt(st.r2[ip])), p) == 0
            end
        end
    end
    check_poly(st)
    seen = Set([st.B])
    for i in 1:(2^(st.s - 1) - 1)
        F.next_poly!(st, i)
        check_poly(st)
        push!(seen, st.B)
    end
    @test length(seen) == 2^(st.s - 1)      # every Gray-code step gives a new B
end
```

- [ ] **Step 2: Run to verify it fails** — `julia --project=julia julia/test/test_polynomials.jl` — Expected: FAIL (`SIQSState` undefined).

- [ ] **Step 3: Implement** — `julia/src/polynomials.jl`

```julia
# SPDX-License-Identifier: GPL-2.0-or-later

mutable struct SIQSState
    kn::BigInt
    fb::FactorBase
    params::Params
    s::Int
    lo::Int                 # 1-based index of first candidate A-prime
    span::Int               # number of candidate A-primes
    target::BigInt          # ideal A ≈ sqrt(2 kn) / M
    aidx::Vector{Int}
    A::BigInt
    B::BigInt
    C::BigInt
    bterms::Vector{BigInt}
    ainv2b::Matrix{Int}     # s × nprimes
    r1::Vector{Int}
    r2::Vector{Int}
    isa::BitVector
end

function SIQSState(kn::BigInt, fb::FactorBase, params::Params)
    np = length(fb.primes)
    s = ndigits(kn, base = 2) ÷ 28 + 1
    target = isqrt(2 * kn) ÷ params.M
    # Port of the C++ window computation (done with 0-based indices there).
    t = integer_root(target, s)
    fact = searchsortedfirst(fb.primes, Int(t)) - 1
    span = np ÷ s ÷ s ÷ 2
    lo = fact - span ÷ 2
    while (fact * fact) ÷ lo - lo < span
        lo -= 1
    end
    lo = max(lo, params.firstprime + 1)
    lo + span + s <= np || error("factor base too small for $s A-primes (lo=$lo, span=$span)")
    return SIQSState(kn, fb, params, s, lo + 1, span, target, zeros(Int, s), BigInt(0),
                     BigInt(0), BigInt(0), BigInt[], zeros(Int, s, np), zeros(Int, np),
                     zeros(Int, np), falses(np))
end

"Pick `s` distinct factor-base indices whose product is within a factor 2 of `st.target`."
function choose_a!(st::SIQSState, rng)
    fb = st.fb
    np = length(fb.primes)
    hi = st.lo + st.span - 1
    while true
        idx = Int[]
        A = BigInt(1)
        while length(idx) < st.s - 1
            i = rand(rng, st.lo:hi)
            i in idx && continue
            push!(idx, i)
            A *= fb.primes[i]
        end
        want = Int(clamp(cld(st.target, A), 2, typemax(Int) ÷ 2))
        j = searchsortedfirst(fb.primes, want)
        for cand in (j, j - 1, j + 1, j - 2, j + 2)
            (st.params.firstprime < cand <= np) || continue
            cand in idx && continue
            Afull = A * fb.primes[cand]
            if 0.5 < Afull / st.target < 2.0
                st.aidx = sort!(vcat(idx, cand))
                return st
            end
        end
    end
end

function a_root!(st::SIQSState, ip::Int)
    p = st.fb.primes[ip]
    Cp = Int(mod(st.C, p))
    Bp = Int(mod(st.B, p))
    r = mod(-Cp * invmod(mod(2 * Bp, p), p), p)
    st.r1[ip] = r
    st.r2[ip] = r
    return nothing
end

"Set up A, the B-terms, the first B and all roots for the A chosen in `st.aidx`."
function init_a!(st::SIQSState)
    fb, kn = st.fb, st.kn
    np = length(fb.primes)
    A = prod(BigInt(fb.primes[i]) for i in st.aidx)
    bterms = BigInt[]
    for i in st.aidx
        q = fb.primes[i]
        Aq = A ÷ q
        g = Int(mod(Aq, q))
        γ = mod(fb.sqrts[i] * invmod(g, q), q)
        γ = min(γ, q - γ)
        push!(bterms, Aq * γ)
    end
    st.A = A
    st.bterms = bterms
    st.B = sum(bterms)
    st.C = (st.B^2 - kn) ÷ A
    fill!(st.isa, false)
    for i in st.aidx
        st.isa[i] = true
    end
    for ip in 1:np
        p = fb.primes[ip]
        if st.isa[ip]
            a_root!(st, ip)
            continue
        end
        Ainv = p == 2 ? 1 : invmod(Int(mod(A, p)), p)
        Bp = Int(mod(st.B, p))
        for j in 1:st.s
            st.ainv2b[j, ip] = mod(2 * Int(mod(bterms[j], p)) * Ainv, p)
        end
        t = fb.sqrts[ip]
        st.r1[ip] = mod((t - Bp) * Ainv, p)
        st.r2[ip] = mod((-t - Bp) * Ainv, p)
    end
    return st
end

"Gray-code step `i` (1 ≤ i < 2^(s-1)): B ← B ± 2·bterm_j, roots shift by ∓ainv2b."
function next_poly!(st::SIQSState, i::Int)
    j = trailing_zeros(i) + 1
    sgn = if ((i >> (j - 1)) & 2) != 0
        st.B += 2 * st.bterms[j]
        -1
    else
        st.B -= 2 * st.bterms[j]
        1
    end
    st.C = (st.B^2 - st.kn) ÷ st.A
    fb = st.fb
    for ip in eachindex(fb.primes)
        if st.isa[ip]
            a_root!(st, ip)
        else
            p = fb.primes[ip]
            d = sgn * st.ainv2b[j, ip]
            st.r1[ip] = mod(st.r1[ip] + d, p)
            st.r2[ip] = mod(st.r2[ip] + d, p)
        end
    end
    return st
end
```

- [ ] **Step 4: Wire up and run** — add `include("polynomials.jl")` and `"test_polynomials.jl"`; run `julia --project=julia julia/test/test_polynomials.jl` — Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add julia && git commit -m "Add SIQS polynomial generation with Gray-code B switching

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Sieve and relation extraction

**Files:**
- Create: `julia/src/sieve.jl`, `julia/src/relations.jl`, `julia/test/test_sieve_relations.jl`
- Modify: `FlintQS.jl` (include `sieve.jl`, `relations.jl`), `runtests.jl`

**Interfaces:**
- Consumes: `SIQSState` and its invariants (Task 4), `Params`.
- Produces:
  - `sieve!(sv::Vector{UInt8}, st)`: zero `sv` (length `2M`) then add `fb.sizes[ip]` at every index `i` (0-based) with `(i - M) ≡ r (mod p)` for `r ∈ {r1,r2}`, for `ip > params.firstprime`.
  - `struct Relation; x::BigInt; fac::Vector{Int32}; end` — meaning `x^2 ≡ ∏ row_value(fac) (mod kn)`, where `fac` lists matrix rows with multiplicity (row 1 = sign −1, row `i+1` = `fb.primes[i]`).
  - `mutable struct RelationStore` with `kn::BigInt`, `fulls::Vector{Relation}`, `combined::Vector{Relation}`, `hub::Dict{Int,Relation}`, `seen::Set{BigInt}`, `npartials::Int`, and constructor `RelationStore(kn)`; `nrelations(store)::Int` = `length(fulls)+length(combined)`.
  - `struct FoundFactor <: Exception; d::BigInt; end` (thrown when a large prime shares a factor with `kn`).
  - `scan!(store, st, sv)`: for every index with `sv[i+1] >= threshold`, call `try_candidate` and file the result.
  - `parity_rows(fac)::Vector{Int32}`: sorted rows with odd multiplicity.
  - `prune_singletons(cols, nrows)::Tuple{Vector{Int},Int}`: indices of kept columns, number of non-empty rows among kept columns.

- [ ] **Step 1: Write the failing test** — `julia/test/test_sieve_relations.jl`

```julia
using Test, Random
using FlintQS
const F = FlintQS
include("test_polynomials_helpers.jl")   # make_state (moved from test_polynomials.jl, see Step 2)

@testset "sieve" begin
    st, rng = make_state(40)
    F.choose_a!(st, rng); F.init_a!(st)
    M = st.params.M
    sv = zeros(UInt8, 2M)
    F.sieve!(sv, st)
    # Spot-check: sieve value at an index equals the summed sizes of sieved primes whose roots hit it.
    for i0 in rand(rng, 0:2M-1, 20)
        x = i0 - M
        expect = 0
        for ip in st.params.firstprime+1:length(st.fb.primes)
            p = st.fb.primes[ip]
            xm = mod(x, p)
            if xm == st.r1[ip] || xm == st.r2[ip]
                expect += st.fb.sizes[ip]
            end
        end
        @test sv[i0+1] == expect
    end
end

@testset "relations" begin
    st, rng = make_state(40)
    store = F.RelationStore(st.kn)
    sv = zeros(UInt8, 2 * st.params.M)
    for _ in 1:3
        F.choose_a!(st, rng); F.init_a!(st)
        for i in 0:(2^(st.s - 1) - 1)
            i > 0 && F.next_poly!(st, i)
            F.sieve!(sv, st)
            F.scan!(store, st, sv)
        end
    end
    @test F.nrelations(store) > 0
    # Every stored relation satisfies x^2 ≡ ∏ (row values) (mod kn)
    rowval(r) = r == 1 ? BigInt(-1) : BigInt(st.fb.primes[r-1])
    for rel in vcat(store.fulls, store.combined)
        lhs = mod(rel.x^2, st.kn)
        rhs = mod(prod(rowval(r) for r in rel.fac; init = BigInt(1)), st.kn)
        @test lhs == rhs
    end
    # parity rows: odd multiplicities only, sorted
    @test F.parity_rows(Int32[3, 5, 3, 7, 7, 7]) == Int32[5, 7]
    cols = [Int32[1, 2], Int32[2, 3], Int32[3, 4], Int32[1, 5]]
    keep, nonempty = F.prune_singletons(cols, 5)
    # cascade: row 5 appears once -> col 4 dropped -> row 1 now once -> col 1 dropped -> rows 2,3,4 ... all pruned
    @test F.prune_singletons(cols, 5) == (Int[], 0)
    @test F.prune_singletons([Int32[1,2], Int32[1,2], Int32[3]], 3) == ([1, 2], 2)
end
```

- [ ] **Step 2: Move `make_state`** from `test_polynomials.jl` into `julia/test/test_polynomials_helpers.jl` (the function body unchanged, no `@testset`), and in `test_polynomials.jl` replace the definition with `include("test_polynomials_helpers.jl")`. The helpers file begins with the same three `using` lines.

- [ ] **Step 3: Run to verify it fails** — `julia --project=julia julia/test/test_sieve_relations.jl` — Expected: FAIL (`sieve!` undefined).

- [ ] **Step 4: Implement `julia/src/sieve.jl`**

```julia
# SPDX-License-Identifier: GPL-2.0-or-later

function sieve!(sv::Vector{UInt8}, st::SIQSState)
    fill!(sv, 0)
    M = st.params.M
    len = length(sv)
    fb = st.fb
    @inbounds for ip in st.params.firstprime+1:length(fb.primes)
        p = fb.primes[ip]
        sz = fb.sizes[ip]
        a = st.r1[ip]
        b = st.r2[ip]
        for j in mod(a + M, p)+1:p:len
            sv[j] += sz
        end
        if b != a
            for j in mod(b + M, p)+1:p:len
                sv[j] += sz
            end
        end
    end
    return sv
end
```

- [ ] **Step 5: Implement `julia/src/relations.jl`**

```julia
# SPDX-License-Identifier: GPL-2.0-or-later

struct Relation
    x::BigInt
    fac::Vector{Int32}
end

struct FoundFactor <: Exception
    d::BigInt
end

mutable struct RelationStore
    kn::BigInt
    fulls::Vector{Relation}
    combined::Vector{Relation}
    hub::Dict{Int,Relation}
    seen::Set{BigInt}
    npartials::Int
end
RelationStore(kn::BigInt) = RelationStore(kn, Relation[], Relation[], Dict{Int,Relation}(), Set{BigInt}(), 0)

nrelations(s::RelationStore) = length(s.fulls) + length(s.combined)

"Divide `p` out of `res` as often as possible; returns (new res, exponent)."
function strip_prime(res::BigInt, p::Int)
    e = 0
    while true
        q, r = divrem(res, p)
        r == 0 || return res, e
        res = q
        e += 1
    end
end

"""
Fully trial-divide Q(x) at sieve index `i0` (0-based). Returns `nothing`, or `(Relation, q)` where
`q == 1` for a full relation and `q > 1` is the single large prime of a partial relation.
"""
function try_candidate(st::SIQSState, sv::Vector{UInt8}, i0::Int)
    P = st.params
    fb = st.fb
    x = i0 - P.M
    t = st.A * x + st.B                      # A x + B
    Q = (t + st.B) * x + st.C                # A x^2 + 2 B x + C;  t^2 ≡ A Q (mod kn)
    res = abs(Q)
    res == 0 && return nothing
    fac = Int32[]
    Q < 0 && push!(fac, Int32(1))
    target = ndigits(res, base = 2) - P.errorbits

    extra = 0
    for ip in 1:P.firstprime
        res, e = strip_prime(res, fb.primes[ip])
        if e > 0
            extra += fb.sizes[ip]
            for _ in 1:e
                push!(fac, Int32(ip + 1))
            end
        end
    end
    Int(sv[i0+1]) + extra < target && return nothing

    for ip in P.firstprime+1:length(fb.primes)
        p = fb.primes[ip]
        xm = mod(x, p)
        (xm == st.r1[ip] || xm == st.r2[ip]) || continue
        res, e = strip_prime(res, p)
        for _ in 1:e
            push!(fac, Int32(ip + 1))
        end
    end
    for ip in st.aidx
        push!(fac, Int32(ip + 1))
    end
    rel = Relation(mod(t, st.kn), fac)
    res == 1 && return rel, 1
    (res <= P.largeprime && res < BigInt(fb.primes[end])^2) || return nothing
    return rel, Int(res)
end

function combine(r1::Relation, r2::Relation, q::Int, kn::BigInt)
    g = gcd(BigInt(q), kn)
    g == 1 || throw(FoundFactor(g))
    x = mod(r1.x * r2.x * invmod(BigInt(q), kn), kn)
    return Relation(x, vcat(r1.fac, r2.fac))
end

function add_candidate!(store::RelationStore, rel::Relation, q::Int)
    rel.x in store.seen && return
    push!(store.seen, rel.x)
    if q == 1
        push!(store.fulls, rel)
    else
        store.npartials += 1
        if haskey(store.hub, q)
            push!(store.combined, combine(rel, store.hub[q], q, store.kn))
        else
            store.hub[q] = rel
        end
    end
    return nothing
end

function scan!(store::RelationStore, st::SIQSState, sv::Vector{UInt8})
    thr = st.params.threshold
    @inbounds for j in eachindex(sv)
        sv[j] >= thr || continue
        out = try_candidate(st, sv, j - 1)
        out === nothing && continue
        add_candidate!(store, out[1], out[2])
    end
    return store
end

"Sorted matrix rows that occur an odd number of times in `fac`."
function parity_rows(fac::AbstractVector{<:Integer})
    s = sort(fac)
    out = Int32[]
    i = 1
    while i <= length(s)
        j = i
        while j < length(s) && s[j+1] == s[i]
            j += 1
        end
        isodd(j - i + 1) && push!(out, Int32(s[i]))
        i = j + 1
    end
    return out
end

"Drop columns that contain a row occurring in no other column (repeatedly)."
function prune_singletons(cols::Vector{Vector{Int32}}, nrows::Int)
    keep = trues(length(cols))
    cnt = zeros(Int, nrows)
    for c in cols, r in c
        cnt[r] += 1
    end
    changed = true
    while changed
        changed = false
        for j in eachindex(cols)
            keep[j] || continue
            if any(r -> cnt[r] == 1, cols[j])
                keep[j] = false
                changed = true
                for r in cols[j]
                    cnt[r] -= 1
                end
            end
        end
    end
    return findall(keep), count(>(0), cnt)
end
```

- [ ] **Step 6: Wire up and run** — add the includes and test file; run `julia --project=julia julia/test/test_sieve_relations.jl`.
Expected: PASS. **Calibration check:** print `F.nrelations(store)` and `store.npartials`; for the 40-digit run over 3 A-values it should be > 0. If it is 0, the sieve threshold is mis-calibrated for this port's sieve scale: lower `THRESHOLDS` for the digit class by 4 at a time until relations appear, and record the adjustment in a comment above the table.

- [ ] **Step 7: Commit**

```bash
git add julia && git commit -m "Add sieve, relation extraction and large-prime merging

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Block Lanczos over GF(2)

**Files:**
- Create: `julia/src/lanczos.jl`, `julia/test/test_lanczos.jl`
- Modify: `FlintQS.jl` (include `lanczos.jl`), `runtests.jl`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `block_lanczos(cols::Vector{Vector{Int32}}, nrows::Int; rng)::Vector{BitVector}` — every returned vector `v` (length `length(cols)`) is nonzero and satisfies `B v = 0` over GF(2), where column `j` of `B` has ones in the rows listed in `cols[j]`. Returns `BitVector[]` when the iteration breaks down (caller retries with a fresh `rng` state). Also internal helpers `mulB!`, `mulBt!`, `tmul`, `mulmat`, `identity64`, `select_cols`, `inv64` used by the tests.
- Representation: an N×64 GF(2) block is a `Vector{UInt64}` of length N (entry i = row i, bit c = column c). A 64×64 matrix is a `Vector{UInt64}` of length 64 (entry r+1 = row r).

- [ ] **Step 1: Write the failing test** — `julia/test/test_lanczos.jl`

```julia
using Test, Random
using FlintQS
const F = FlintQS

dense(cols, nrows) = [Int(r) in c for r in 1:nrows, c in cols]   # nrows × N Bool

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
    bits(M) = [((M[r] >> c) & 1) == 1 for r in 1:64, c in 0:63]
    @test bits(F.mulmat(X, Y)) == (Int.(bits(X)) * Int.(bits(Y))) .% 2 .== 1
    @test F.mulmat(X, F.identity64()) == X
    V = rand(rng, UInt64, 200); W = rand(rng, UInt64, 200)
    bv(V) = [((V[i] >> c) & 1) == 1 for i in 1:200, c in 0:63]
    @test bits(F.tmul(V, W)) == (Int.(bv(V))' * Int.(bv(W))) .% 2 .== 1
    # inv64 of a random invertible matrix
    M = F.identity64()
    for _ in 1:200; i, j = rand(rng, 1:64, 2); i != j && (M[i] ⊻= M[j]); end
    @test F.mulmat(M, F.inv64(M)) == F.identity64()
    # select_cols: symmetric matrix of rank 40 -> mask has 40 bits and the block is invertible
    G = rand(rng, UInt64, 64)
    T = F.tmul(vcat(G[1:40], zeros(UInt64, 24)), vcat(G[1:40], zeros(UInt64, 24)))   # symmetric, rank ≤ 40
    mask, winv = F.select_cols(T, typemax(UInt64))
    @test count_ones(mask) == 40
    sub = F.mulmat(winv, T)       # winv * T restricted: equals identity on the selected columns
    for c in 0:63
        (mask >> c) & 1 == 1 && @test (F.mulmat(winv, T)[c+1] & mask) == (UInt64(1) << c)
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
```

- [ ] **Step 2: Run to verify it fails** — `julia --project=julia julia/test/test_lanczos.jl` — Expected: FAIL (`mulmat` undefined).

- [ ] **Step 3: Implement** — `julia/src/lanczos.jl`

```julia
# SPDX-License-Identifier: GPL-2.0-or-later
# Montgomery's block Lanczos over GF(2), 64 vectors at a time, for A = Bᵀ B.

identity64() = UInt64[UInt64(1) << r for r in 0:63]

"out = B*v where column j of B has ones in rows cols[j]; v has one UInt64 per column."
function mulB!(out::Vector{UInt64}, cols, v::Vector{UInt64})
    fill!(out, 0)
    @inbounds for j in eachindex(cols)
        x = v[j]
        x == 0 && continue
        for r in cols[j]
            out[r] ⊻= x
        end
    end
    return out
end

"out = Bᵀ*w."
function mulBt!(out::Vector{UInt64}, cols, w::Vector{UInt64})
    @inbounds for j in eachindex(cols)
        acc = UInt64(0)
        for r in cols[j]
            acc ⊻= w[r]
        end
        out[j] = acc
    end
    return out
end

"Rows-times-matrix: `X` (any number of 64-bit rows) times a 64×64 matrix `M`."
function mulmat(X::Vector{UInt64}, M::Vector{UInt64})
    out = zeros(UInt64, length(X))
    @inbounds for i in eachindex(X)
        x = X[i]
        acc = UInt64(0)
        while x != 0
            acc ⊻= M[trailing_zeros(x)+1]
            x &= x - 1
        end
        out[i] = acc
    end
    return out
end

"Xᵀ Y for two N×64 blocks (a 64×64 matrix)."
function tmul(X::Vector{UInt64}, Y::Vector{UInt64})
    M = zeros(UInt64, 64)
    @inbounds for i in eachindex(X)
        x = X[i]
        y = Y[i]
        while x != 0
            M[trailing_zeros(x)+1] ⊻= y
            x &= x - 1
        end
    end
    return M
end

"Inverse of a nonsingular 64×64 GF(2) matrix (Gauss-Jordan), or `nothing`."
function inv64(M::Vector{UInt64})
    a = copy(M)
    b = identity64()
    for c in 0:63
        p = 0
        for r in c:63
            if (a[r+1] >> c) & 1 == 1
                p = r
                break
            end
        end
        p == 0 && (a[c+1] >> c) & 1 == 0 && return nothing
        a[p+1], a[c+1] = a[c+1], a[p+1]
        b[p+1], b[c+1] = b[c+1], b[p+1]
        for r in 0:63
            r != c && (a[r+1] >> c) & 1 == 1 && (a[r+1] ⊻= a[c+1]; b[r+1] ⊻= b[c+1])
        end
    end
    return b
end

"""
Choose a maximal set S of columns (bitmask) with `T[S,S]` nonsingular, trying the columns in
`prefer` first, and return `(mask, winv)` where `winv = S (SᵀTS)⁻¹ Sᵀ`.
`T` must be symmetric.
"""
function select_cols(T::Vector{UInt64}, prefer::UInt64)
    order = vcat([c for c in 0:63 if (prefer >> c) & 1 == 1],
                 [c for c in 0:63 if (prefer >> c) & 1 == 0])
    R = copy(T)                       # Schur complement on the not-yet-selected columns
    mask = UInt64(0)
    progress = true
    while progress
        progress = false
        for c in order                # 1×1 pivots
            (mask >> c) & 1 == 1 && continue
            if (R[c+1] >> c) & 1 == 1
                piv = R[c+1]
                for r in 0:63
                    (R[r+1] >> c) & 1 == 1 && (R[r+1] ⊻= piv)
                end
                mask |= UInt64(1) << c
                progress = true
            end
        end
        progress && continue
        for c in order                # 2×2 pivot on an alternating remainder
            (mask >> c) & 1 == 1 && continue
            rc = R[c+1]
            rc == 0 && continue
            d = trailing_zeros(rc)
            rd = R[d+1]
            for r in 0:63
                x = R[r+1]
                nx = x
                (x >> c) & 1 == 1 && (nx ⊻= rd)
                (x >> d) & 1 == 1 && (nx ⊻= rc)
                R[r+1] = nx
            end
            mask |= (UInt64(1) << c) | (UInt64(1) << d)
            progress = true
            break
        end
    end
    sub = UInt64[(mask >> r) & 1 == 1 ? (T[r+1] & mask) : (UInt64(1) << r) for r in 0:63]
    inv = inv64(sub)
    inv === nothing && return mask, nothing
    winv = UInt64[(mask >> r) & 1 == 1 ? (inv[r+1] & mask) : UInt64(0) for r in 0:63]
    return mask, winv
end

"""
Null vectors of B (as `BitVector`s of length `length(cols)`), or `BitVector[]` if the Lanczos
iteration breaks down. B is `nrows × length(cols)`; column j has ones at rows `cols[j]`.
"""
function block_lanczos(cols::Vector{Vector{Int32}}, nrows::Int; rng = Random.default_rng())
    N = length(cols)
    tmp = zeros(UInt64, nrows)
    function mulA(v)
        mulB!(tmp, cols, v)
        return mulBt!(zeros(UInt64, N), cols, tmp)
    end
    I64 = identity64()
    Y = rand(rng, UInt64, N)
    V = mulA(Y)
    V0 = copy(V)
    X = zeros(UInt64, N)
    Vm1 = Vm2 = nothing
    prev = nothing                    # (winv, vtav, vtav2, mask) of step i-1
    winv_pp = nothing                 # winv of step i-2
    prev_mask = nothing
    for _ in 1:(N ÷ 32 + 20)
        AV = mulA(V)
        vtav = tmul(V, AV)
        all(iszero, vtav) && return nullvectors(cols, nrows, X .⊻ Y, V)
        vtav2 = tmul(AV, AV)
        prefer = prev_mask === nothing ? typemax(UInt64) : ~prev_mask
        mask, winv = select_cols(vtav, prefer)
        (mask == 0 || winv === nothing) && return BitVector[]
        prev_mask !== nothing && (prefer & ~mask) != 0 && return BitVector[]
        X .⊻= mulmat(V, mulmat(winv, tmul(V, V0)))
        D = I64 .⊻ mulmat(winv, (vtav2 .& mask) .⊻ vtav)
        Vn = (AV .& mask) .⊻ mulmat(V, D)
        if prev !== nothing
            E = mulmat(prev.winv, vtav .& mask)
            Vn .⊻= mulmat(Vm1, E)
        end
        if winv_pp !== nothing
            inner = (prev.vtav2 .& prev.mask) .⊻ prev.vtav
            F = mulmat(winv_pp, mulmat(I64 .⊻ mulmat(prev.vtav, prev.winv), inner) .& mask)
            Vn .⊻= mulmat(Vm2, F)
        end
        winv_pp = prev === nothing ? nothing : prev.winv
        prev = (winv = winv, vtav = vtav, vtav2 = vtav2, mask = mask)
        prev_mask = mask
        Vm2, Vm1, V = Vm1, V, Vn
    end
    return BitVector[]
end

"Combine the 128 columns of [Z V] into vectors in ker B by GF(2) elimination on B·[Z V]."
function nullvectors(cols, nrows, Z::Vector{UInt64}, V::Vector{UInt64})
    N = length(cols)
    blocks = (Z, V)
    BW = [mulB!(zeros(UInt64, nrows), cols, b) for b in blocks]
    basis = Dict{Int,Tuple{BitVector,UInt128}}()
    result = BitVector[]
    for j in 0:127
        blk, bit = (j >> 6) + 1, j & 63
        v = BitVector([(BW[blk][r] >> bit) & 1 == 1 for r in 1:nrows])
        comb = UInt128(1) << j
        while true
            p = findfirst(v)
            if p === nothing
                vec = falses(N)
                for jj in 0:127
                    (comb >> jj) & 1 == 1 || continue
                    b2, bt2 = (jj >> 6) + 1, jj & 63
                    for i in 1:N
                        (blocks[b2][i] >> bt2) & 1 == 1 && (vec[i] = !vec[i])
                    end
                end
                any(vec) && push!(result, vec)
                break
            end
            if haskey(basis, p)
                bv, bc = basis[p]
                v .⊻= bv
                comb ⊻= bc
            else
                basis[p] = (v, comb)
                break
            end
        end
    end
    return result
end
```

- [ ] **Step 4: Wire up and run** — add `include("lanczos.jl")` and `"test_lanczos.jl"`; run `julia --project=julia julia/test/test_lanczos.jl`.
Expected: PASS. If `select_cols` or the recurrence tests fail, debug in this order: (a) `mulmat`/`tmul` tests, (b) `select_cols` rank test, (c) `block_lanczos` on the smallest matrix (164×100). The `E` and `F` formulas follow Montgomery (1995): `E = Wᵢ₋₁⁻¹ (VᵢᵀAVᵢ) SᵢSᵢᵀ`, `F = Wᵢ₋₂⁻¹ (I + Vᵢ₋₁ᵀAVᵢ₋₁ Wᵢ₋₁⁻¹)(Vᵢ₋₁ᵀA²Vᵢ₋₁ Sᵢ₋₁Sᵢ₋₁ᵀ + Vᵢ₋₁ᵀAVᵢ₋₁) SᵢSᵢᵀ` (all over GF(2), so minus = plus).

- [ ] **Step 5: Commit**

```bash
git add julia && git commit -m "Add block Lanczos null-space solver over GF(2)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 7: SIQS driver (`siqs_split`)

**Files:**
- Modify: `julia/src/FlintQS.jl`
- Create: `julia/test/test_flintqs.jl` (the part below; extended in Task 8)

**Interfaces:**
- Consumes: everything above.
- Produces: `siqs_split(n::BigInt; rng, verbose)::BigInt` — a nontrivial divisor of `n` (n composite, odd-or-even, not a perfect power, >= 40 digits, no factors < 1000). Internal: `collect_relations!`, `dependencies_to_factor`, `try_dependency`.

- [ ] **Step 1: Write the failing test** — start `julia/test/test_flintqs.jl`

```julia
using Test, Random
using FlintQS
const F = FlintQS

function semiprime(rng, d1, d2)
    p = F.next_prime(BigInt(10)^(d1 - 1) + rand(rng, 1:10^(d1 - 2)))
    q = F.next_prime(BigInt(10)^(d2 - 1) + rand(rng, 1:10^(d2 - 2)))
    p == q ? semiprime(rng, d1, d2) : (p * q, p, q)
end

@testset "siqs_split" begin
    rng = Xoshiro(2024)
    for (d1, d2) in ((20, 21), (25, 26))        # 41- and 51-digit semiprimes
        n, p, q = semiprime(rng, d1, d2)
        d = F.siqs_split(n; rng = rng)
        @test 1 < d < n && n % d == 0
        @test d in (p, q)
    end
end
```

- [ ] **Step 2: Run to verify it fails** — `julia --project=julia julia/test/test_flintqs.jl` — Expected: FAIL (`siqs_split` undefined).

- [ ] **Step 3: Implement the driver** — append to `julia/src/FlintQS.jl` before `end # module`

```julia
include("modarith.jl")
include("params.jl")
include("factorbase.jl")
include("polynomials.jl")
include("sieve.jl")
include("relations.jl")
include("lanczos.jl")

export flintqs, siqs_split

"Sieve with fresh A values until the store holds at least `target` relations."
function collect_relations!(store::RelationStore, st::SIQSState, target::Int, rng, verbose::Bool)
    sv = zeros(UInt8, 2 * st.params.M)
    curves = 0
    while nrelations(store) < target
        choose_a!(st, rng)
        init_a!(st)
        for i in 0:(2^(st.s - 1) - 1)
            i > 0 && next_poly!(st, i)
            sieve!(sv, st)
            scan!(store, st, sv)
            curves += 1
        end
        verbose && println("curves=$curves full=$(length(store.fulls)) combined=$(length(store.combined)) ",
                           "partials=$(store.npartials) target=$target")
    end
    return store
end

"Square-root step for one dependency; returns a nontrivial divisor of `n` or `nothing`."
function try_dependency(rels::Vector{Relation}, dep::BitVector, fb::FactorBase, kn::BigInt, n::BigInt)
    counts = zeros(Int, length(fb.primes) + 1)
    X = BigInt(1)
    for j in findall(dep)
        X = mod(X * rels[j].x, kn)
        for r in rels[j].fac
            counts[r] += 1
        end
    end
    any(isodd, counts) && return nothing
    Y = BigInt(1)
    for row in 2:length(counts)
        c = counts[row]
        c == 0 && continue
        Y = mod(Y * powermod(BigInt(fb.primes[row-1]), c ÷ 2, kn), kn)
    end
    for g in (gcd(X - Y, n), gcd(X + Y, n))
        1 < g < n && return g
    end
    return nothing
end

function dependencies_to_factor(store::RelationStore, fb::FactorBase, kn::BigInt, n::BigInt, rng)
    rels = vcat(store.fulls, store.combined)
    nrows = length(fb.primes) + 1
    cols = [parity_rows(r.fac) for r in rels]
    keep, nonempty = prune_singletons(cols, nrows)
    length(keep) - nonempty < 32 && return nothing
    kept_rels = rels[keep]
    kept_cols = cols[keep]
    for _ in 1:5
        for dep in block_lanczos(kept_cols, nrows; rng = rng)
            d = try_dependency(kept_rels, dep, fb, kn, n)
            d !== nothing && return d
        end
    end
    return nothing
end

"""
    siqs_split(n; rng, verbose) -> BigInt

Return a nontrivial divisor of the composite `n` using the self-initializing quadratic sieve.
`n` must have at least 40 decimal digits, no prime factor below 1000 and not be a perfect power.
"""
function siqs_split(n::BigInt; rng = Random.default_rng(), verbose::Bool = false)
    P = params_for(ndigits(n))
    k = knuth_schroeppel(n)
    kn = k * n
    fb = build_factorbase(kn, k, P.nprimes)
    st = SIQSState(kn, fb, P)
    store = RelationStore(kn)
    target = P.nprimes + 64
    try
        while true
            collect_relations!(store, st, target, rng, verbose)
            d = dependencies_to_factor(store, fb, kn, n, rng)
            d !== nothing && return d
            target += 64
        end
    catch e
        e isa FoundFactor || rethrow()
        g = gcd(e.d, n)
        1 < g < n && return g
        rethrow()
    end
end
```
(Remove the earlier `include("util.jl")`-only body duplication: keep `include("util.jl")` first, then the includes above in order.)

- [ ] **Step 4: Run to verify it passes** — `julia --project=julia julia/test/test_flintqs.jl`
Expected: PASS in well under a minute for 41 and 51 digits. If it never terminates, run with `verbose = true` and inspect the relation counts: (a) zero relations → threshold calibration (Task 5 Step 6); (b) relations but no factor → `try_dependency` returning `nothing` repeatedly means Lanczos vectors are not in ker B (re-run Task 6 tests) or the sign row is mishandled.

- [ ] **Step 5: Commit**

```bash
git add julia && git commit -m "Add SIQS driver: relation collection, linear algebra, square-root step

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Public `flintqs` API, small-input handling, CLI

**Files:**
- Modify: `julia/src/FlintQS.jl`, `julia/test/test_flintqs.jl`, `julia/test/runtests.jl`
- Create: `julia/bin/flintqs.jl`

**Interfaces:**
- Produces: `flintqs(n::Integer; rng, verbose)::Vector{BigInt}` — sorted prime factors with multiplicity; `[]` for 1; `[n]` for primes; `ArgumentError` for `n < 1`.
- Internal: `split_composite(m::BigInt; rng, verbose)::BigInt` — nontrivial divisor of composite non-power `m` with no factor < 1000, any size (<= 18 digits: rho; 19-39 digits: lift; >= 40: `siqs_split`).

- [ ] **Step 1: Add failing tests** — append to `julia/test/test_flintqs.jl`

```julia
@testset "flintqs" begin
    rng = Xoshiro(99)
    # trivial inputs
    @test F.flintqs(1) == BigInt[]
    @test F.flintqs(97) == [97]
    @test_throws ArgumentError F.flintqs(0)
    @test_throws ArgumentError F.flintqs(-5)
    # small-factor stripping and ordinary factorization
    @test F.flintqs(2^5 * 3^3 * 7 * 1000003; rng = rng) == BigInt.([2, 2, 2, 2, 2, 3, 3, 3, 7, 1000003])
    # <= 18 digits (rho)
    n, p, q = semiprime(rng, 8, 9)
    @test F.flintqs(n; rng = rng) == sort([p, q])
    # Review Focus 2: prime powers
    P = F.next_prime(big"10"^24)
    @test F.flintqs(P^2; rng = rng) == [P, P]
    @test F.flintqs(P^3; rng = rng) == [P, P, P]
    # Review Focus 3: many small factors times a large prime
    @test F.flintqs(big"2"^5 * big"3"^3 * P; rng = rng) == vcat(fill(BigInt(2), 5), fill(BigInt(3), 3), [P])
    # Review Focus 4: 39-digit (lift path) and 40-digit (direct) semiprimes
    for (d1, d2) in ((19, 20), (20, 20))
        n, p, q = semiprime(rng, d1, d2)
        @test F.flintqs(n; rng = rng) == sort([p, q])
    end
    # Review Focus 5: very unbalanced semiprime (7-digit × 45-digit)
    n, p, q = semiprime(rng, 7, 45)
    @test F.flintqs(n; rng = rng) == sort([p, q])
    # product of three primes, 50 digits
    p1 = F.next_prime(big"10"^16 + 12345); p2 = F.next_prime(big"10"^16 + 99999); p3 = F.next_prime(big"10"^17 + 777)
    @test F.flintqs(p1 * p2 * p3; rng = rng) == sort([p1, p2, p3])
end
```

- [ ] **Step 2: Run to verify it fails** — `julia --project=julia julia/test/test_flintqs.jl` — Expected: FAIL (`flintqs` undefined).

- [ ] **Step 3: Implement** — append to `julia/src/FlintQS.jl` (inside the module, after `siqs_split`)

```julia
"Nontrivial divisor of a composite, non-power `m` with no prime factor below 1000."
function split_composite(m::BigInt; rng, verbose::Bool = false)
    d = ndigits(m)
    if d <= 18
        g = pollard_rho(m; rng = rng)
        g === nothing && error("Pollard rho failed on a ≤18-digit number $m")
        return g
    elseif d >= 40
        return siqs_split(m; rng = rng, verbose = verbose)
    end
    # 19–39 digits: multiply by a prime so the work number has ≥ 40 digits, then keep gcds with m.
    while true
        P = next_prime(BigInt(10)^(41 - d) + rand(rng, 1:10^(40 - d)))
        N = m * P
        g = siqs_split(N; rng = rng, verbose = verbose)
        for c in (gcd(g, m), gcd(N ÷ g, m))
            1 < c < m && return c
        end
    end
end

"""
    flintqs(n; rng, verbose) -> Vector{BigInt}

Prime factorization of `n ≥ 1` as a sorted vector with multiplicity (`[]` for `n == 1`).
Strips factors below 1000, reduces perfect powers, and splits the rest with Pollard rho
(≤ 18 digits) or the self-initializing quadratic sieve.
"""
function flintqs(n::Integer; rng = Random.default_rng(), verbose::Bool = false)
    n >= 1 || throw(ArgumentError("flintqs requires n ≥ 1, got $n"))
    m = BigInt(n)
    out = BigInt[]
    for p in SMALL_PRIMES
        while m % p == 0
            push!(out, BigInt(p))
            m ÷= p
        end
    end
    stack = BigInt[m]
    while !isempty(stack)
        c = pop!(stack)
        c == 1 && continue
        if is_probable_prime(c)
            push!(out, c)
        elseif (pp = perfect_power(c)) !== nothing
            b, k = pp
            for _ in 1:k
                push!(stack, b)
            end
        else
            d = split_composite(c; rng = rng, verbose = verbose)
            push!(stack, d, c ÷ d)
        end
    end
    return sort!(out)
end
```

`julia/bin/flintqs.jl`:
```julia
#!/usr/bin/env julia
# Usage: julia --project=julia julia/bin/flintqs.jl <n> [-v]
using FlintQS
isempty(ARGS) && (println(stderr, "usage: flintqs.jl <n> [-v]"); exit(1))
n = parse(BigInt, ARGS[1])
t = @elapsed fs = flintqs(n; verbose = "-v" in ARGS)
println(join(fs, " * "), "   (", round(t, digits = 3), " s)")
```

- [ ] **Step 4: Run to verify it passes** — `julia --project=julia julia/test/test_flintqs.jl`
Expected: PASS. Prime-power cases: `P^2` has 49 digits so it must be caught by `perfect_power` before SIQS (SIQS would never split it).

- [ ] **Step 5: Run the whole suite and the CLI**

Run: `julia --project=julia -e 'using Pkg; Pkg.test()'` — Expected: all testsets PASS.
Run: `julia --project=julia julia/bin/flintqs.jl $(julia --project=julia -e 'using FlintQS; print(FlintQS.next_prime(big"10"^25)*FlintQS.next_prime(big"10"^26))')` — Expected: two primes printed.

- [ ] **Step 6: Commit**

```bash
git add julia && git commit -m "Add flintqs() with small-input handling, prime powers, and CLI

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Benchmark script

**Files:**
- Create: `julia/bench/bench.jl`

**Interfaces:** Consumes `flintqs`. Produces a Markdown table on stdout (digits, seconds, peak RSS MB).

- [ ] **Step 1: Write the script**

```julia
# Usage: julia --project=julia julia/bench/bench.jl [digits...]   (default: 40 45 50 55 60)
using FlintQS, Random

function semiprime(rng, digits)
    h = digits ÷ 2
    p = FlintQS.next_prime(BigInt(10)^(h - 1) + rand(rng, 1:10^(h - 2)))
    q = FlintQS.next_prime(BigInt(10)^(digits - h - 1) + rand(rng, 1:10^(digits - h - 2)))
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
```

- [ ] **Step 2: Run it** — `julia --project=julia julia/bench/bench.jl 40 45 50` — Expected: a table with three rows and no assertion error. Save the output of the default run (`40 45 50 55 60`) for Task 10's `benchmarks.md`.

- [ ] **Step 3: Commit**

```bash
git add julia/bench && git commit -m "Add benchmark script

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 10: Documentation with MaterialDocs.jl

**Files:**
- Create: `julia/docs/Project.toml`, `julia/docs/make.jl`, `julia/docs/src/index.md`, `julia/docs/src/algorithm.md`, `julia/docs/src/api.md`, `julia/docs/src/benchmarks.md`
- Modify: `julia/src/FlintQS.jl` (module docstring), `julia/src/util.jl` (docstrings already present for exported-API-relevant helpers)

**Interfaces:** Consumes the exported API (`flintqs`, `siqs_split`) and its docstrings.

- [ ] **Step 1: Docs environment**

Run:
```bash
cd julia/docs && julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.add(["Documenter", "MaterialDocs"])'
```
Expected: resolves and installs without error (MaterialDocs is in the General registry; source `github.com/mthelm85/MaterialDocs.jl`). `docs/Project.toml` and `Manifest.toml` are created (Manifest is gitignored).

- [ ] **Step 2: Module docstring** — in `julia/src/FlintQS.jl` put this immediately above `module FlintQS`:

```julia
"""
    FlintQS

Julia port of William Hart's FlintQS: a self-initializing quadratic sieve with single
large-prime merging and Block Lanczos linear algebra.

Exports [`flintqs`](@ref) (full prime factorization) and [`siqs_split`](@ref) (one divisor
of a ≥ 40-digit composite).
"""
```

- [ ] **Step 3: `julia/docs/make.jl`**

```julia
using Documenter, MaterialDocs, FlintQS

DocMeta.setdocmeta!(FlintQS, :DocTestSetup, :(using FlintQS); recursive = true)

makedocs(;
    modules = [FlintQS],
    authors = "technocrat",
    sitename = "FlintQS.jl",
    format = Material3(theme = :ocean_depth, dark_mode = :toggle),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "How it works" => "algorithm.md",
        "API Reference" => "api.md",
        "Benchmarks" => "benchmarks.md",
    ],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(; repo = "github.com/technocrat/FlintQS", devbranch = "julia-port")
end
```

- [ ] **Step 4: Pages**

`julia/docs/src/index.md`:
````markdown
# FlintQS.jl

A Julia port of [FlintQS](https://github.com/sagemath/FlintQS), William Hart's quadratic sieve
(now obsolete upstream; superseded by `qsieve_factor` in FLINT).

```julia
using FlintQS
flintqs(big"1000000000000000000000000000057" * big"1000000000000000000000000000063")
```

## Installation

```julia
using Pkg
Pkg.develop(path = "julia")      # from a checkout of this repository
```

## Scope

- Factors any positive integer: small factors by trial division, perfect powers by integer
  roots, ≤ 18-digit composites by Pollard rho, larger composites by SIQS.
- SIQS parameters are the FlintQS tables for 40–91 digits. Inputs of 19–39 digits are
  multiplied by a prime to reach 40 digits.
- Relations are held in memory; there are no temporary files.

## License

GPL-2.0-or-later, as the upstream code it is ported from.
````

`julia/docs/src/algorithm.md`:
````markdown
# How it works

| Stage | Julia | FlintQS C++ |
|:--|:--|:--|
| Multiplier, factor base | `knuth_schroeppel`, `build_factorbase` | `knuthSchroeppel`, `computeFactorBase` |
| Parameters | `params_for` | tuning tables in `QS.cpp` |
| Polynomials | `choose_a!`, `init_a!`, `next_poly!` | polynomial selection in `mainRoutine` |
| Sieve | `sieve!` | `sieveInterval`, `sieve2` |
| Relations | `try_candidate`, `scan!` | `evaluateSieve` |
| Large primes | `combine`, `add_candidate!` | `lprels.cpp` |
| Linear algebra | `prune_singletons`, `block_lanczos` | `lanczos.cpp` |

## Polynomials

For ``A`` a product of ``s`` factor-base primes and ``B^2 \equiv kn \pmod A``, set
``Q(x) = Ax^2 + 2Bx + C`` with ``C = (B^2 - kn)/A``. Then
``(Ax + B)^2 \equiv A\,Q(x) \pmod{kn}``. Switching ``B`` by Gray code gives ``2^{s-1}``
polynomials per ``A`` at the cost of a few machine-integer operations per prime.

## Relations

A sieve index passes when its accumulated log size reaches the table threshold. Candidates are
fully trial-divided over the factor base. A cofactor of 1 is a full relation; a prime cofactor
below the large-prime bound is a *partial*; two partials with the same large prime ``q`` merge
into one relation after dividing ``x`` by ``q``.

The sign of ``Q(x)`` is matrix row 1. (The C++ ignores sign; this port tracks it.)

## Linear algebra

Columns containing a row that occurs nowhere else are removed, then Montgomery's block Lanczos
(64 vectors at a time, over ``A = B^T B``) produces vectors in the kernel of the exponent-parity
matrix. Each is turned into a congruence of squares ``X^2 \equiv Y^2`` and ``\gcd(X - Y, n)``.
````

`julia/docs/src/api.md`:
````markdown
# API Reference

```@docs
FlintQS
flintqs
siqs_split
```
````

`julia/docs/src/benchmarks.md`: paste the Markdown table produced by `julia --project=julia julia/bench/bench.jl` (default sizes), preceded by a sentence stating the machine (`Sys.cpu_info()[1].model`, Julia version) and that times include one `flintqs` call per size after a warm-up. No placeholder text may remain.

- [ ] **Step 5: Build the docs**

Run: `julia --project=julia/docs julia/docs/make.jl`
Expected: build succeeds with no warnings about missing docstrings (`checkdocs = :exports`); `julia/docs/build/index.html` exists. Open it and confirm the Material 3 theme renders and the search box works.

- [ ] **Step 6: Commit**

```bash
git add julia/docs julia/src && git commit -m "Add MaterialDocs documentation

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:** util (T1), modarith (T2), params/factorbase (T3), polynomials (T4), sieve + relations (T5), lanczos (T6), driver (T7-8), tests per file (T1-8), benchmark (T9), docs (T10). Clean-room and standalone constraints in Global Constraints. Spec items updated in this plan versus the spec text: FlintQS is a *self-initializing* QS (spec said multiple-polynomial); inputs below 40 digits use rho/lift because the C++ tables start at 40 digits; sign is tracked (spec silent); license is GPL-2.0-or-later (spec silent). The spec is amended to match in the same commit as this plan.

**Placeholders:** none; the only open item is the threshold calibration in Task 5 Step 6, which has a concrete procedure.

**Type consistency:** `Relation(x, fac::Vector{Int32})`, rows `1 = sign`, `i+1 = fb.primes[i]` used identically in `try_candidate`, `combine`, `parity_rows`, `try_dependency` (`counts` indexed by row, prime index `row-1`). `SIQSState.aidx` are 1-based FB indices everywhere. `Params.firstprime` = count of unsieved leading primes in `sieve!`, `try_candidate`, and `choose_a!` (A-primes `> firstprime`).

**Review Focus coverage:** all five items have tests in Task 8 (`flintqs` testset).
