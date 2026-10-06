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
        p = -1
        for r in c:63
            if (a[r+1] >> c) & 1 == 1
                p = r
                break
            end
        end
        p < 0 && return nothing
        a[p+1], a[c+1] = a[c+1], a[p+1]
        b[p+1], b[c+1] = b[c+1], b[p+1]
        for r in 0:63
            if r != c && (a[r+1] >> c) & 1 == 1
                a[r+1] ⊻= a[c+1]
                b[r+1] ⊻= b[c+1]
            end
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
        for c in order
            (mask >> c) & 1 == 1 && continue
            if (R[c+1] >> c) & 1 == 1                 # 1×1 pivot
                piv = R[c+1]
                for r in 0:63
                    (R[r+1] >> c) & 1 == 1 && (R[r+1] ⊻= piv)
                end
                mask |= UInt64(1) << c
                progress = true
            elseif (prefer >> c) & 1 == 1 && R[c+1] != 0
                # Preferred column with zero Schur diagonal: pair it with any coupled column d.
                # P = [[0,1],[1,e]] is nonsingular with inverse [[e,1],[1,0]] over GF(2).
                d = trailing_zeros(R[c+1])
                rc, rd = R[c+1], R[d+1]
                e = (rd >> d) & 1 == 1
                for r in 0:63
                    x = R[r+1]
                    add = UInt64(0)
                    (x >> c) & 1 == 1 && (add ⊻= (e ? rc : UInt64(0)) ⊻ rd)
                    (x >> d) & 1 == 1 && (add ⊻= rc)
                    R[r+1] = x ⊻ add
                end
                mask |= (UInt64(1) << c) | (UInt64(1) << d)
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
        # Krylov space exhausted (A is singular): salvage whatever the current state yields.
        # `nullvectors` only returns vectors it has verified to satisfy B v = 0.
        (mask == 0 || winv === nothing) && return nullvectors(cols, nrows, X .⊻ Y, V)
        prev_mask !== nothing && (prefer & ~mask) != 0 && return nullvectors(cols, nrows, X .⊻ Y, V)
        X .⊻= mulmat(V, mulmat(winv, tmul(V, V0)))
        D = I64 .⊻ mulmat(winv, (vtav2 .& mask) .⊻ vtav)
        Vn = (AV .& mask) .⊻ mulmat(V, D)
        if prev !== nothing
            E = mulmat(prev.winv, vtav .& mask)
            Vn .⊻= mulmat(Vm1, E)
        end
        if winv_pp !== nothing
            inner = (prev.vtav2 .& prev.mask) .⊻ prev.vtav
            Fm = mulmat(winv_pp, mulmat(I64 .⊻ mulmat(prev.vtav, prev.winv), inner) .& mask)
            Vn .⊻= mulmat(Vm2, Fm)
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
