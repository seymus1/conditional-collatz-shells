# Conditional Collatz Shells

Supporting code for the paper:

**Conditional Delayed Descent in Non-Certified Collatz Shells via Dyadic Gaps and Residue-Cover Trees**

by **Seyma Yaman Kayadibi**.

## Overview

This project studies delayed first descent for the accelerated odd Collatz map

T(n) = (3n + 1) / 2^v2(3n + 1)

using canonical shells, dyadic gaps, valuation cylinders, and residue-cover trees.

The central conditional route is:

Global Dyadic Gap + Non-Certified Tree Closure  
=> Universal First Descent  
=> Collatz Convergence.

## Code

Main scripts:

scripts/canonical_noncertified_shell_verification.py

scripts/dyadic_gap_verification.py

The first script performs exact-integer finite-range verification for canonical non-certified shell elements.

The second script verifies the dyadic gap condition and checks the upper continued-fraction obstruction levels.

Paper parameters:

N = 10,000,000  
L = 10,000  
M_MAX = 100,000

To run:

python3 scripts/canonical_noncertified_shell_verification.py

python3 scripts/dyadic_gap_verification.py

Both scripts use only the Python standard library.

## Status

This repository provides finite-range computational support for the conditional framework developed in the paper.

It does not constitute an unconditional proof of the Collatz conjecture.

## Author

Seyma Yaman Kayadibi  
Melbourne, Australia

## License

MIT License
