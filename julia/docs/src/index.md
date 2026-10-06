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
