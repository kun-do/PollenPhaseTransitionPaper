#include "../PhaseDiagramCalculations/Common/wigner_compat.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

bool nearly_equal(double actual, double expected, double tolerance = 1.0e-12)
{
    return std::fabs(actual - expected) <= tolerance;
}

int fail(const char* label, double actual, double expected)
{
    std::cerr << label << " failed: actual=" << actual << " expected=" << expected << std::endl;
    return 1;
}

} // namespace

int main()
{
    wig_table_init(200, 3);
    wig_temp_init(200);

    if (!nearly_equal(wig3jj(0, 0, 0, 0, 0, 0), 1.0)) {
        return fail("origin", wig3jj(0, 0, 0, 0, 0, 0), 1.0);
    }

    const double minus_one_over_root_three = -1.0 / std::sqrt(3.0);
    if (!nearly_equal(wig3jj(2, 2, 0, 0, 0, 0), minus_one_over_root_three)) {
        return fail("closed_form_m0", wig3jj(2, 2, 0, 0, 0, 0), minus_one_over_root_three);
    }

    const double plus_one_over_root_three = 1.0 / std::sqrt(3.0);
    if (!nearly_equal(wig3jj(2, 2, 0, 2, -2, 0), plus_one_over_root_three)) {
        return fail("closed_form_m1", wig3jj(2, 2, 0, 2, -2, 0), plus_one_over_root_three);
    }

    if (wig3jj(2, 2, 0, 2, 0, 0) != 0.0) {
        return fail("m_sum_rule", wig3jj(2, 2, 0, 2, 0, 0), 0.0);
    }

    if (wig3jj(2, 2, 6, 0, 0, 0) != 0.0) {
        return fail("triangle_rule", wig3jj(2, 2, 6, 0, 0, 0), 0.0);
    }

    if (wig3jj(3, 2, 1, 0, 0, 0) != 0.0) {
        return fail("parity_rule", wig3jj(3, 2, 1, 0, 0, 0), 0.0);
    }

    const double value = wig3jj(2, 2, 0, 2, -2, 0);
    const double reflected = wig3jj(2, 2, 0, -2, 2, 0);
    if (!nearly_equal(value, reflected)) {
        return fail("reflection_symmetry", reflected, value);
    }

    const double odd_sum_value = wig3jj(2, 2, 2, 2, 0, -2);
    const double odd_sum_negated = wig3jj(2, 2, 2, -2, 0, 2);
    if (!nearly_equal(odd_sum_value, -odd_sum_negated)) {
        return fail("sign_symmetry", odd_sum_value, -odd_sum_negated);
    }

    for (int two_j = 0; two_j <= 20; ++two_j) {
        for (int two_m = -two_j; two_m <= two_j; two_m += 2) {
            const double stable = wig3jj(two_j, two_j, 0, two_m, -two_m, 0);
            if (!std::isfinite(stable)) {
                std::cerr << "non-finite value at two_j=" << two_j << " two_m=" << two_m << std::endl;
                return 1;
            }
        }
    }

    wig_temp_free();
    wig_table_free();
    return 0;
}
