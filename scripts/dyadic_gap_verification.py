# paper3_dyadic_gap_verification_100k.py
# Exact dyadic gap verification for Paper 3
# Range: 2 <= m <= 100000
# Precision: 800 decimal digits

from decimal import Decimal, getcontext
import csv
import time


# ============================================================
# Parameters
# ============================================================

M_MAX = 100_000
DECIMAL_PRECISION = 800
CF_TERMS = 120

getcontext().prec = DECIMAL_PRECISION

LOG2_3 = Decimal(3).ln() / Decimal(2).ln()


# ============================================================
# Basic functions
# ============================================================

def ceil_decimal(x: Decimal) -> int:
    """
    Exact ceiling for a Decimal value.
    """
    n = int(x)

    if Decimal(n) == x:
        return n

    return n + 1


def K_m(m: int) -> int:
    """
    Compute K_m = ceil(m log_2 3).
    """
    return ceil_decimal(Decimal(m) * LOG2_3)


def dyadic_gap_data(m: int) -> dict:
    """
    Compute dyadic gap data for a given m.

    Dyadic gap condition:
        2^{K_m} - 3^m >= 2^{K_m-m}

    where:
        K_m = ceil(m log_2 3).
    """
    K = K_m(m)

    two_K = 1 << K
    three_m = pow(3, m)

    gap = two_K - three_m
    lower_bound = 1 << (K - m)
    margin = gap - lower_bound

    return {
        "m": m,
        "K": K,
        "gap_holds": gap >= lower_bound,
        "gap_minus_lower_sign": 1 if margin > 0 else 0 if margin == 0 else -1,
        "gap_bit_length": gap.bit_length() if gap > 0 else None,
        "lower_bound_bit_length": lower_bound.bit_length(),
        "margin_bit_length": margin.bit_length() if margin > 0 else None,
        "K_valid": two_K >= three_m and (1 << (K - 1)) < three_m,
    }


# ============================================================
# Continued fraction functions
# ============================================================

def continued_fraction_decimal(x: Decimal, n_terms: int) -> list[int]:
    """
    Compute the continued fraction expansion of Decimal x.
    """
    terms = []
    y = x

    for _ in range(n_terms):
        a = int(y)
        terms.append(a)

        frac = y - Decimal(a)

        if frac == 0:
            break

        y = Decimal(1) / frac

    return terms


def convergents_from_cf(cf: list[int]) -> list[tuple[int, int]]:
    """
    Return convergents p/q from a continued fraction list.
    """
    convergents = []

    p_minus2, p_minus1 = 0, 1
    q_minus2, q_minus1 = 1, 0

    for a in cf:
        p = a * p_minus1 + p_minus2
        q = a * q_minus1 + q_minus2

        convergents.append((p, q))

        p_minus2, p_minus1 = p_minus1, p
        q_minus2, q_minus1 = q_minus1, q

    return convergents


def upper_convergents_log2_3(q_limit: int) -> list[dict]:
    """
    Return upper convergents p/q to log_2 3 with q <= q_limit.
    """
    cf = continued_fraction_decimal(LOG2_3, CF_TERMS)
    convergents = convergents_from_cf(cf)

    rows = []

    for j, (p, q) in enumerate(convergents):
        if q > q_limit:
            break

        value = Decimal(p) / Decimal(q)

        if value > LOG2_3 and q >= 2:
            next_partial = cf[j + 1] if j + 1 < len(cf) else None
            epsilon = Decimal(p) - Decimal(q) * LOG2_3

            rows.append({
                "j": j,
                "p": p,
                "q": q,
                "next_partial_quotient": next_partial,
                "epsilon_p_minus_q_alpha": str(epsilon),
            })

    return rows


# ============================================================
# Direct dyadic gap verification
# ============================================================

def run_direct_gap_test(m_max: int, csv_path: str) -> dict:
    """
    Verify the dyadic gap condition for all 2 <= m <= m_max.

    Saves detailed results to a CSV file.
    """
    failures = []
    invalid_K = []
    checked = 0

    start = time.time()

    fieldnames = [
        "m",
        "K",
        "gap_holds",
        "gap_minus_lower_sign",
        "gap_bit_length",
        "lower_bound_bit_length",
        "margin_bit_length",
        "K_valid",
    ]

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for m in range(2, m_max + 1):
            row = dyadic_gap_data(m)
            writer.writerow(row)

            checked += 1

            if not row["gap_holds"]:
                failures.append(row)

            if not row["K_valid"]:
                invalid_K.append(row)

            if m % 5000 == 0:
                elapsed = time.time() - start
                print(f"Checked up to m={m:,} | elapsed={elapsed:.2f}s")

    elapsed = time.time() - start

    return {
        "checked": checked,
        "failures": failures,
        "invalid_K": invalid_K,
        "elapsed_seconds": elapsed,
    }


# ============================================================
# Upper convergent verification
# ============================================================

def run_upper_convergent_tests(q_limit: int, csv_path: str) -> list[dict]:
    """
    Test the dyadic gap at all upper continued-fraction convergent
    denominators q <= q_limit.
    """
    rows = upper_convergents_log2_3(q_limit)
    output_rows = []

    for row in rows:
        q = row["q"]
        p = row["p"]

        gap_info = dyadic_gap_data(q)

        output = {
            "j": row["j"],
            "p": str(p),
            "q": str(q),
            "next_partial_quotient": str(row["next_partial_quotient"]),
            "epsilon_p_minus_q_alpha": row["epsilon_p_minus_q_alpha"],
            "gap_holds_at_q": gap_info["gap_holds"],
            "K_valid_at_q": gap_info["K_valid"],
            "K_q": gap_info["K"],
            "direct_margin_D_sign": gap_info["gap_minus_lower_sign"],
            "direct_margin_D_bit_length": gap_info["margin_bit_length"],
        }

        output_rows.append(output)

    fieldnames = [
        "j",
        "p",
        "q",
        "next_partial_quotient",
        "epsilon_p_minus_q_alpha",
        "gap_holds_at_q",
        "K_valid_at_q",
        "K_q",
        "direct_margin_D_sign",
        "direct_margin_D_bit_length",
    ]

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for output in output_rows:
            writer.writerow(output)

    return output_rows


# ============================================================
# Summary file
# ============================================================

def write_summary(
    summary_path: str,
    direct_result: dict,
    upper_rows: list[dict],
) -> None:
    """
    Write a text summary for Paper 3 computational reporting.
    """
    cf = continued_fraction_decimal(LOG2_3, CF_TERMS)

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Paper 3 dyadic gap / continued fraction verification\n")
        f.write(f"Precision: {DECIMAL_PRECISION} decimal digits\n")
        f.write(f"Direct dyadic gap test: 2 <= m <= {M_MAX}\n")
        f.write(f"Direct checked levels: {direct_result['checked']}\n")
        f.write(f"Direct failures: {len(direct_result['failures'])}\n")
        f.write(f"Invalid K count: {len(direct_result['invalid_K'])}\n")
        f.write(f"Continued fraction terms: {CF_TERMS}\n")
        f.write(f"Upper convergent q limit: q <= {M_MAX}\n")
        f.write("Initial continued fraction terms:\n")
        f.write(str(cf[:50]) + "\n")
        f.write(f"Upper convergent rows within q limit: {len(upper_rows)}\n\n")

        f.write("Upper convergent checks:\n")

        for row in upper_rows:
            f.write(str(row) + "\n")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    direct_csv = "paper3_direct_dyadic_gap_100k.csv"
    upper_csv = "paper3_upper_convergent_gap_tests_100k.csv"
    summary_txt = "paper3_dyadic_gap_summary_100k.txt"

    print("Paper 3 Dyadic Gap Verification")
    print("-------------------------------")
    print(f"M_MAX = {M_MAX:,}")
    print(f"Precision = {DECIMAL_PRECISION} decimal digits")
    print(f"log2(3) = {LOG2_3}")

    direct_result = run_direct_gap_test(M_MAX, direct_csv)

    print("\nSummary")
    print("-------")
    print(f"Checked m = 2 to {M_MAX:,}")
    print(f"Elapsed seconds: {direct_result['elapsed_seconds']:.2f}")

    if direct_result["failures"]:
        print("Status: FAILURE FOUND")
        for failure in direct_result["failures"][:5]:
            print(failure)
    else:
        print("Status: NO FAILURES FOUND")

    if direct_result["invalid_K"]:
        print("Invalid K values found:")
        for row in direct_result["invalid_K"][:5]:
            print(row)
    else:
        print("Invalid K count: 0")

    upper_rows = run_upper_convergent_tests(M_MAX, upper_csv)

    print("\nUpper convergent checks:")
    for row in upper_rows:
        print(row)

    write_summary(summary_txt, direct_result, upper_rows)

    print("\nFiles written:")
    print(f"- {direct_csv}")
    print(f"- {upper_csv}")
    print(f"- {summary_txt}")
