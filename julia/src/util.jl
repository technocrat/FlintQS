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

