/**
 * instance.cpp
 *
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <random>
#include <vector>
#include <exception>

#include "xorshift.hpp"

using namespace std;

static void sampleGraphWithGroundTruth(
    vector<int> &x0,
    vector<Eigen::Triplet<unsigned>> &entries,
    unsigned n,
    double a,
    double b,
    unsigned seed)
{
  if (a < 0.0 || a > n)
    throw invalid_argument("a is out of range");
  if (b < 0.0 || b > n)
    throw invalid_argument("b is out of range");

  xorshift prng(seed);
  bernoulli_distribution faircoin(0.5);
  uniform_real_distribution<double> unif(0.0, 1.0);

  x0.clear();
  x0.reserve(n);

  // generate ground truth
  for (unsigned i = 0; i < n; i++)
    x0.push_back(faircoin(prng) ? +1.0 : -1.0);

  const double inp = a/static_cast<double>(n);
  const double outp = b/static_cast<double>(n);

  for (unsigned src = 0; src < n; src++) {
    for (unsigned dst = src; dst < n; dst++) {
      const double r = unif(prng);
      const int x0ij = x0[src] * x0[dst];
      if ((x0ij == +1 && r <= inp) || (x0ij == -1 && r <= outp)) {
        entries.emplace_back(src, dst, 1);
        if (src != dst)
          entries.emplace_back(dst, src, 1);
      }
    }
  }
}

int main(int argc, char **argv)
{

  const unsigned n = 128000;
  const double a = 20.0;
  const double b = 3.0;

  vector<int> x0;
  vector<Eigen::Triplet<unsigned>> entries;
  sampleGraphWithGroundTruth(x0, entries, n, a, b, 8323753);

  return 0;
}
