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
