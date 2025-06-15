#include "utility.h"

/**
 * Removing leading and trailing whitespaces of a string
 */
std::string trim_spaces(const std::string& line)
{
    const char* WhiteSpace = " \t\v\r\n";
    std::size_t start = line.find_first_not_of(WhiteSpace);
    // If the whole string is white spaces, return an empty string
    if (start == std::string::npos)
        return std::string();
    std::size_t end = line.find_last_not_of(WhiteSpace);
    return line.substr(start, end - start + 1);
}

/**
 * Define laplacian, use idx to preserve periodicity
*/
Grid laplacian(const Grid& f, double h) {
    int n = f.size();
    Grid result = zero_grid(n);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) default(none) shared(f, result, n, h)
#endif //_OPENMP
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = (f[idx(i + 1, n)][j] + f[idx(i - 1, n)][j]
                            + f[i][idx(j + 1, n)] + f[i][idx(j - 1, n)]
                            - 4.0 * f[i][j]) / (h * h);
        }
    }
    return result;
}

/**
 * Genrerate a random double number
 */
double rand_real(double lower_bound, double upper_bound)
{
    static std::default_random_engine re{std::random_device{}()};
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    return unif(re);
}