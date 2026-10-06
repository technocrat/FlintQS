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
polynomials per ``A`` at the cost of a few machine-integer operations per prime. Primes
dividing ``kn`` (square root 0) are never used in ``A``.

## Relations

A sieve index passes when its accumulated log size reaches the table threshold. Candidates are
fully trial-divided over the factor base. A cofactor of 1 is a full relation; a prime cofactor
below the large-prime bound is a *partial*; two partials with the same large prime ``q`` merge
into one relation after dividing ``x`` by ``q``.

The sign of ``Q(x)`` is matrix row 1. (The C++ ignores sign; this port tracks it.)

## Linear algebra

Columns containing a row that occurs nowhere else are removed, then Montgomery's block Lanczos
(64 vectors at a time, over ``A = B^T B``) produces vectors in the kernel of the exponent-parity
matrix. If the Krylov space is exhausted early, the vectors found so far are combined and
verified rather than discarded. Each kernel vector is turned into a congruence of squares
``X^2 \equiv Y^2`` and ``\gcd(X - Y, n)`` gives a factor.
