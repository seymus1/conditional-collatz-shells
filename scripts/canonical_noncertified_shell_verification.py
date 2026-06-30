#!/usr/bin/env python3
"""
Canonical-shell exact verification for delayed-descent coverage
of non-certified accelerated odd Collatz spines.

Paper 3 target:
    Delayed-Descent Coverage of Non-Certified Accelerated Odd Collatz Spines

This script counts each odd integer only once, at its canonical shell level

    M(n) = v2(n + 1).

It verifies, up to a finite bound N, whether canonical non-certified shell
elements descend below their initial value after the first possible exit step.

This is finite exact-integer evidence only.
It does not prove the Collatz conjecture.
"""

from statistics import median
import csv
import time

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

N = 10_000_000
STEP_LIMIT = 10_000
PROGRESS_EVERY = 250_000

WRITE_CSV = True

# ------------------------------------------------------------
# Basic arithmetic
# ------------------------------------------------------------

def v2(x: int) -> int:
    """
    Return v_2(x), the exponent of the highest power of 2 dividing x.
    """
    if x <= 0:
        raise ValueError("v2(x) requires x > 0.")
    return (x & -x).bit_length() - 1

def T_odd(n: int) -> int:
    """
    Accelerated odd Collatz map:
        T(n) = (3n + 1) / 2^{v2(3n + 1)}.
    """
    x = 3 * n + 1
    return x >> v2(x)

def first_descent_time(n: int, limit: int = STEP_LIMIT):
    """
    Return tau(n) = min{r >= 1 : T^r(n) < n},
    if found within the iteration limit.

    Return None if unresolved within the limit.
    """
    x = n
    for r in range(1, limit + 1):
        x = T_odd(x)
        if x < n:
            return r
    return None

def kappa(m: int) -> int:
    """
    Compute kappa_m as the smallest integer k >= 0 such that

        2^k (2^m - 1) >= 3^m.

    This avoids floating-point logarithms.
    """
    k = 0
    lhs_base = (1 << m) - 1
    rhs = pow(3, m)

    while (1 << k) * lhs_base < rhs:
        k += 1

    return k

def q_star(m: int, k: int) -> int:
    """
    Compute q_m^* = (3^m)^(-1) mod 2^k.
    """
    modulus = 1 << k
    return pow(pow(3, m), -1, modulus)

def canonical_shell_level(n: int) -> int:
    """
    Return M(n) = v2(n + 1).
    """
    return v2(n + 1)

def canonical_parameter(n: int, m: int) -> int:
    """
    If n in S_m^circ, then n = 2^m q - 1 with q odd.
    Return q.
    """
    return (n + 1) >> m

# ------------------------------------------------------------
# Cumulative valuation certificate
# ------------------------------------------------------------

def cumulative_data_until(n: int, r: int):
    """
    Compute A_r(n), B_r(n), and T^r(n) exactly.

    Formula:
        T^r(n) = (3^r n + B_r) / 2^{A_r}

    where
        A_r = sum_{i=0}^{r-1} v2(3T^i(n)+1)

    and
        B_r = sum_{j=0}^{r-1} 3^{r-1-j} 2^{A_j}.
    """
    x = n
    valuations = []

    for _ in range(r):
        y = 3 * x + 1
        a = v2(y)
        valuations.append(a)
        x = y >> a

    A_values = [0]
    total = 0
    for a in valuations:
        total += a
        A_values.append(total)

    A_r = A_values[-1]

    B_r = 0
    for j in range(r):
        B_r += pow(3, r - 1 - j) * (1 << A_values[j])

    lhs = pow(3, r) * n + B_r
    rhs = 1 << A_r

    if lhs % rhs != 0:
        raise RuntimeError("Cumulative formula divisibility check failed.")

    T_r = lhs // rhs

    return {
        "A_r": A_r,
        "B_r": B_r,
        "T_r": T_r,
        "valuations": valuations,
    }

def descent_certificate_holds(n: int, r: int) -> bool:
    """
    Check the cumulative descent inequality:

        B_r < n(2^{A_r} - 3^r).

    This is equivalent to T^r(n) < n.
    """
    data = cumulative_data_until(n, r)
    A_r = data["A_r"]
    B_r = data["B_r"]

    return B_r < n * ((1 << A_r) - pow(3, r))

# ------------------------------------------------------------
# Main computation
# ------------------------------------------------------------

def main() -> None:
    start = time.time()

    # Determine admissible shell levels.
    max_m = 0
    while (1 << (max_m + 1)) - 1 <= N:
        max_m += 1

    print(f"N = {N:,}")
    print(f"STEP_LIMIT = {STEP_LIMIT:,}")
    print(f"Admissible nontrivial shell levels: 2 <= m <= {max_m}")
    print()

    # Precompute kappa_m and q_m^*.
    kappas = {}
    qstars = {}

    for m in range(2, max_m + 1):
        k = kappa(m)
        kappas[m] = k
        qstars[m] = q_star(m, k)

    # Canonical shell counts.
    shell_rows = []

    total_shell_cases = 0
    total_certified = 0
    total_noncertified = 0

    for m in range(2, max_m + 1):
        qmax = (N + 1) // (1 << m)

        # Canonical shell requires q odd.
        shell_count = (qmax + 1) // 2

        k = kappas[m]
        modulus = 1 << k
        qs = qstars[m]

        # q_s is odd; q_s + modulus*t remains odd.
        if qs <= qmax:
            certified_count = (qmax - qs) // modulus + 1
        else:
            certified_count = 0

        noncertified_count = shell_count - certified_count

        shell_rows.append(
            {
                "m": m,
                "kappa": k,
                "modulus": modulus,
                "shell_count": shell_count,
                "certified": certified_count,
                "noncertified": noncertified_count,
            }
        )

        total_shell_cases += shell_count
        total_certified += certified_count
        total_noncertified += noncertified_count

    # Per-m dynamic statistics for non-certified shell cases.
    per_m_stats = {
        m: {
            "noncertified": 0,
            "resolved": 0,
            "unresolved": 0,
            "early": 0,
            "exact_at_m": 0,
            "delayed": 0,
            "tau_values": [],
            "E_values": [],
        }
        for m in range(2, max_m + 1)
    }

    # Global dynamic statistics.
    resolved = 0
    unresolved = 0
    early = 0
    exact_at_m = 0
    delayed = 0

    tau_values = []
    E_values = []

    # Cumulative certificate sanity check.
    certificate_failures = 0

    processed_noncert = 0

    print("Starting exact iteration over canonical non-certified shell cases...")
    print()

    # n ≡ 3 mod 4 are exactly canonical shells with m >= 2.
    for n in range(3, N + 1, 4):
        m = canonical_shell_level(n)
        q = canonical_parameter(n, m)

        if q % 2 != 1:
            raise RuntimeError(f"Canonical parameter q is not odd for n={n}.")

        k = kappas[m]
        modulus = 1 << k
        qs = qstars[m]

        certified = (q % modulus == qs)

        if certified:
            continue

        per_m_stats[m]["noncertified"] += 1
        processed_noncert += 1

        tau_n = first_descent_time(n)

        if tau_n is None:
            unresolved += 1
            per_m_stats[m]["unresolved"] += 1
        else:
            resolved += 1
            tau_values.append(tau_n)
            per_m_stats[m]["resolved"] += 1
            per_m_stats[m]["tau_values"].append(tau_n)

            E = tau_n - m
            E_values.append(E)
            per_m_stats[m]["E_values"].append(E)

            if tau_n < m:
                early += 1
                per_m_stats[m]["early"] += 1
            elif tau_n == m:
                exact_at_m += 1
                per_m_stats[m]["exact_at_m"] += 1
            else:
                delayed += 1
                per_m_stats[m]["delayed"] += 1

            # Verify cumulative descent certificate at the found first descent time.
            if not descent_certificate_holds(n, tau_n):
                certificate_failures += 1

        if processed_noncert % PROGRESS_EVERY == 0:
            elapsed = time.time() - start
            print(
                f"Processed non-certified cases: {processed_noncert:,} "
                f"| elapsed: {elapsed:.1f}s"
            )

    # ------------------------------------------------------------
    # Print canonical shell count table.
    # ------------------------------------------------------------

    print()
    print("Canonical-shell count table up to N =", f"{N:,}")
    print("-" * 86)
    print(
        f"{'m':>3} {'kappa':>6} {'2^kappa':>8} "
        f"{'|S_m^circ|':>14} {'|D_m^circ|':>14} {'|N_m^circ|':>14}"
    )
    print("-" * 86)

    for row in shell_rows:
        print(
            f"{row['m']:>3} "
            f"{row['kappa']:>6} "
            f"{row['modulus']:>8} "
            f"{row['shell_count']:>14,} "
            f"{row['certified']:>14,} "
            f"{row['noncertified']:>14,}"
        )

    print("-" * 86)
    print(
        f"{'Total':>3} {'--':>6} {'--':>8} "
        f"{total_shell_cases:>14,} "
        f"{total_certified:>14,} "
        f"{total_noncertified:>14,}"
    )

    # ------------------------------------------------------------
    # Print global delayed-descent summary.
    # ------------------------------------------------------------

    print()
    print("Canonical non-certified delayed-descent summary")
    print("-" * 86)
    print(f"Search bound N: {N:,}")
    print(f"STEP_LIMIT: {STEP_LIMIT:,}")
    print(f"Canonical shell cases with m >= 2: {total_shell_cases:,}")
    print(f"Certified canonical-shell cases: {total_certified:,}")
    print(f"Canonical non-certified shell cases: {total_noncertified:,}")
    print(f"Processed non-certified cases: {processed_noncert:,}")
    print(f"Resolved first-descent times: {resolved:,}")
    print(f"Unresolved cases: {unresolved:,}")
    print(f"Early descents tau < m: {early:,}")
    print(f"Exact descents tau = m: {exact_at_m:,}")
    print(f"Delayed descents tau > m: {delayed:,}")
    print(f"Cumulative certificate failures at tau: {certificate_failures:,}")

    if tau_values:
        print(f"Mean tau: {sum(tau_values) / len(tau_values):.10f}")
        print(f"Median tau: {median(tau_values)}")
        print(f"Maximum tau: {max(tau_values)}")

    if E_values:
        print(f"Mean E = tau - m: {sum(E_values) / len(E_values):.10f}")
        print(f"Median E: {median(E_values)}")
        print(f"Maximum E: {max(E_values)}")

    # ------------------------------------------------------------
    # Per-m summary.
    # ------------------------------------------------------------

    print()
    print("Per-shell delayed-descent summary")
    print("-" * 110)
    print(
        f"{'m':>3} {'noncert':>10} {'resolved':>10} {'unresolved':>10} "
        f"{'exact_m':>10} {'delayed':>10} {'mean_tau':>12} {'mean_E':>12} {'max_E':>8}"
    )
    print("-" * 110)

    per_m_csv_rows = []

    for m in range(2, max_m + 1):
        stats = per_m_stats[m]
        taus = stats["tau_values"]
        Es = stats["E_values"]

        mean_tau = sum(taus) / len(taus) if taus else None
        mean_E = sum(Es) / len(Es) if Es else None
        max_E = max(Es) if Es else None

        print(
            f"{m:>3} "
            f"{stats['noncertified']:>10,} "
            f"{stats['resolved']:>10,} "
            f"{stats['unresolved']:>10,} "
            f"{stats['exact_at_m']:>10,} "
            f"{stats['delayed']:>10,} "
            f"{mean_tau if mean_tau is not None else 0:>12.6f} "
            f"{mean_E if mean_E is not None else 0:>12.6f} "
            f"{max_E if max_E is not None else 0:>8}"
        )

        per_m_csv_rows.append(
            {
                "m": m,
                "noncertified": stats["noncertified"],
                "resolved": stats["resolved"],
                "unresolved": stats["unresolved"],
                "early": stats["early"],
                "exact_at_m": stats["exact_at_m"],
                "delayed": stats["delayed"],
                "mean_tau": mean_tau,
                "median_tau": median(taus) if taus else None,
                "max_tau": max(taus) if taus else None,
                "mean_E": mean_E,
                "median_E": median(Es) if Es else None,
                "max_E": max_E,
            }
        )

    # ------------------------------------------------------------
    # Write CSV outputs.
    # ------------------------------------------------------------

    if WRITE_CSV:
        with open("paper3_canonical_shell_counts.csv", "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "m",
                    "kappa",
                    "modulus",
                    "shell_count",
                    "certified",
                    "noncertified",
                ],
            )
            writer.writeheader()
            writer.writerows(shell_rows)

        with open("paper3_per_m_delayed_summary.csv", "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "m",
                    "noncertified",
                    "resolved",
                    "unresolved",
                    "early",
                    "exact_at_m",
                    "delayed",
                    "mean_tau",
                    "median_tau",
                    "max_tau",
                    "mean_E",
                    "median_E",
                    "max_E",
                ],
            )
            writer.writeheader()
            writer.writerows(per_m_csv_rows)

        with open("paper3_global_summary.txt", "w") as f:
            f.write(f"N = {N}\n")
            f.write(f"STEP_LIMIT = {STEP_LIMIT}\n")
            f.write(f"Canonical shell cases m >= 2 = {total_shell_cases}\n")
            f.write(f"Certified canonical-shell cases = {total_certified}\n")
            f.write(f"Canonical non-certified shell cases = {total_noncertified}\n")
            f.write(f"Resolved = {resolved}\n")
            f.write(f"Unresolved = {unresolved}\n")
            f.write(f"Early tau < m = {early}\n")
            f.write(f"Exact tau = m = {exact_at_m}\n")
            f.write(f"Delayed tau > m = {delayed}\n")
            f.write(f"Cumulative certificate failures at tau = {certificate_failures}\n")

            if tau_values:
                f.write(f"Mean tau = {sum(tau_values) / len(tau_values):.10f}\n")
                f.write(f"Median tau = {median(tau_values)}\n")
                f.write(f"Maximum tau = {max(tau_values)}\n")

            if E_values:
                f.write(f"Mean E = {sum(E_values) / len(E_values):.10f}\n")
                f.write(f"Median E = {median(E_values)}\n")
                f.write(f"Maximum E = {max(E_values)}\n")

        print()
        print("CSV/text outputs written:")
        print("  paper3_canonical_shell_counts.csv")
        print("  paper3_per_m_delayed_summary.csv")
        print("  paper3_global_summary.txt")

    elapsed = time.time() - start

    print()
    print(f"Finished in {elapsed:.2f} seconds.")
    print()
    print("Important:")
    print("This is finite-range exact-integer evidence only.")
    print("It supports the non-certified coverage problem but does not prove Collatz.")

if __name__ == "__main__":
    main()

