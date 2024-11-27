#include <iostream>
#include <fstream>
#include <random>
#include <format>

int main() {
    const double lower_bound = -10000, upper_bound = 10000;
    int n;
    double x, y;
    std::cin >> n;
    std::ofstream output(std::format("sample_{}.txt", n));
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;
    for (int i = 0; i < n; ++i) {
        x = unif(re), y = unif(re);
        output << x << ' ' << y << '\n';
    }
    output.close();

    return 0;
}