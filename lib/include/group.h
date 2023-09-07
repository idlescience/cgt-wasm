#ifndef GROUP_LIBRARY_H
#define GROUP_LIBRARY_H

#include <vector>

/**
 * Implementation of the computation of the game theoretic Shapley value.
 */
namespace Group
{
    /** IMPLEMENTATION */

    static unsigned binomial_coef(unsigned n, unsigned k)
    {
        if (k > n)
        {
            return 0;
        }
        if (k * 2 > n)
        {
            k = n - k;
        }
        if (k == 0)
        {
            return 1;
        }

        int result = n;
        for (int i = 2; i <= k; ++i)
        {
            result *= (n - i + 1);
            result /= i;
        }
        return result;
    }

} // namespace Group

#endif
