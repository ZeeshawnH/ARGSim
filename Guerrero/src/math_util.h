/**
 * math_util.hpp
 * Contains bionmial coefficient method
 *
 * @author Zeeshawn Hasnain
 */


double binomial_coefficient(int n, int k) {
    if (k > n) return 0;
    double result = 1.0;
    
    for (int i = 1; i <= k; ++i) {
        result *= static_cast<double>(n - (k - i)) / i;
    }
    return result;
}
