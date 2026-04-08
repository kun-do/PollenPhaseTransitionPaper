#include "wigner_compat.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

std::vector<long double>& log_factorials()
{
    static std::vector<long double> cache(1, 0.0L);
    return cache;
}

void ensure_log_factorials(int max_argument)
{
    if (max_argument < 0) {
        return;
    }

    std::vector<long double>& cache = log_factorials();
    if (max_argument < static_cast<int>(cache.size())) {
        return;
    }

    const int start = static_cast<int>(cache.size());
    cache.resize(max_argument + 1);
    for (int i = start; i <= max_argument; ++i) {
        cache[i] = std::lgammal(static_cast<long double>(i) + 1.0L);
    }
}

long double log_factorial(int n)
{
    ensure_log_factorials(n);
    return log_factorials()[n];
}

bool valid_quantum_pair(int two_j, int two_m)
{
    if (two_j < 0) {
        return false;
    }
    if (std::abs(two_m) > two_j) {
        return false;
    }
    return ((two_j - two_m) % 2) == 0;
}

int half_to_int(int doubled_value)
{
    return doubled_value / 2;
}

long double phase_from_integer(int exponent)
{
    return (exponent % 2) == 0 ? 1.0L : -1.0L;
}

} // namespace

void wig_table_init(int max_two_j, int /*wigner_type*/)
{
    const int cache_limit = std::max(0, (3 * max_two_j) / 2 + 4);
    ensure_log_factorials(cache_limit);
}

void wig_temp_init(int max_two_j)
{
    wig_table_init(max_two_j, 3);
}

void wig_thread_temp_init(int max_two_j)
{
    wig_temp_init(max_two_j);
}

double wig3jj(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3)
{
    if (two_m1 + two_m2 + two_m3 != 0) {
        return 0.0;
    }
    if (!valid_quantum_pair(two_j1, two_m1) ||
        !valid_quantum_pair(two_j2, two_m2) ||
        !valid_quantum_pair(two_j3, two_m3)) {
        return 0.0;
    }
    if (two_j3 < std::abs(two_j1 - two_j2) || two_j3 > (two_j1 + two_j2)) {
        return 0.0;
    }
    if (((two_j1 + two_j2 + two_j3) % 2) != 0) {
        return 0.0;
    }

    const int triangle_a = half_to_int(two_j1 + two_j2 - two_j3);
    const int triangle_b = half_to_int(two_j1 - two_j2 + two_j3);
    const int triangle_c = half_to_int(-two_j1 + two_j2 + two_j3);
    const int triangle_d = half_to_int(two_j1 + two_j2 + two_j3) + 1;

    if (triangle_a < 0 || triangle_b < 0 || triangle_c < 0) {
        return 0.0;
    }

    const int j1_plus_m1 = half_to_int(two_j1 + two_m1);
    const int j1_minus_m1 = half_to_int(two_j1 - two_m1);
    const int j2_plus_m2 = half_to_int(two_j2 + two_m2);
    const int j2_minus_m2 = half_to_int(two_j2 - two_m2);
    const int j3_plus_m3 = half_to_int(two_j3 + two_m3);
    const int j3_minus_m3 = half_to_int(two_j3 - two_m3);

    const int z_min = std::max(
        0,
        std::max(
            half_to_int(two_j2 - two_j3 - two_m1),
            half_to_int(two_j1 + two_m2 - two_j3)));
    const int z_max = std::min(
        triangle_a,
        std::min(j1_minus_m1, j2_plus_m2));

    if (z_min > z_max) {
        return 0.0;
    }

    const int cache_limit = std::max(
        triangle_d,
        std::max(
            std::max(j1_plus_m1, j1_minus_m1),
            std::max(
                std::max(j2_plus_m2, j2_minus_m2),
                std::max(j3_plus_m3, j3_minus_m3))));
    ensure_log_factorials(cache_limit);

    const long double log_prefactor =
        0.5L * (
            log_factorial(triangle_a) +
            log_factorial(triangle_b) +
            log_factorial(triangle_c) -
            log_factorial(triangle_d) +
            log_factorial(j1_plus_m1) +
            log_factorial(j1_minus_m1) +
            log_factorial(j2_plus_m2) +
            log_factorial(j2_minus_m2) +
            log_factorial(j3_plus_m3) +
            log_factorial(j3_minus_m3));

    long double max_term_log = -std::numeric_limits<long double>::infinity();
    std::vector<long double> term_logs;
    std::vector<long double> term_signs;
    term_logs.reserve(z_max - z_min + 1);
    term_signs.reserve(z_max - z_min + 1);

    for (int z = z_min; z <= z_max; ++z) {
        const int denominator_a = z;
        const int denominator_b = triangle_a - z;
        const int denominator_c = j1_minus_m1 - z;
        const int denominator_d = j2_plus_m2 - z;
        const int denominator_e = half_to_int(two_j3 - two_j2 + two_m1) + z;
        const int denominator_f = half_to_int(two_j3 - two_j1 - two_m2) + z;

        if (denominator_b < 0 || denominator_c < 0 || denominator_d < 0 ||
            denominator_e < 0 || denominator_f < 0) {
            continue;
        }

        const long double term_log =
            -(
                log_factorial(denominator_a) +
                log_factorial(denominator_b) +
                log_factorial(denominator_c) +
                log_factorial(denominator_d) +
                log_factorial(denominator_e) +
                log_factorial(denominator_f));

        term_logs.push_back(term_log);
        term_signs.push_back(phase_from_integer(z));
        max_term_log = std::max(max_term_log, term_log);
    }

    if (!std::isfinite(static_cast<double>(max_term_log))) {
        return 0.0;
    }

    long double scaled_sum = 0.0L;
    for (std::size_t i = 0; i < term_logs.size(); ++i) {
        scaled_sum += term_signs[i] * std::exp(term_logs[i] - max_term_log);
    }

    if (scaled_sum == 0.0L) {
        return 0.0;
    }

    const int phase_exponent = half_to_int(two_j1 - two_j2 - two_m3);
    const long double prefactor_sign = phase_from_integer(phase_exponent);
    const long double combined = prefactor_sign * std::exp(log_prefactor + max_term_log) * scaled_sum;

    if (std::abs(combined) < 1.0e-18L) {
        return 0.0;
    }

    return static_cast<double>(combined);
}

void wig_temp_free()
{
}

void wig_table_free()
{
    std::vector<long double>& cache = log_factorials();
    cache.assign(1, 0.0L);
}
