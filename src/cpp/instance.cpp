/**
 * instance.cpp
 *
 */

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <exception>

#include <vector>
#include <stdexcept>
#include <utility>

#include <getopt.h>
#include <unistd.h>

#include "timer.hpp"

typedef double SamplingFloatType;

// Boost.Random's implementation of mt19937 seems to be faster than the
// built-in one from in <random>. Unfortunately, Boost.Random generators do not
// seem to be directly compatible with C++11 distributions, so we have this
// mess.
#ifdef USE_BOOST_RANDOM
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
typedef ::boost::random::mt19937 PRNG;
static ::boost::random::bernoulli_distribution<SamplingFloatType> faircoin(0.5);
static ::boost::random::uniform_real_distribution<SamplingFloatType> unif(0.0, 1.0);
#else
#include <random>
typedef ::std::minstd_rand PRNG;
static ::std::bernoulli_distribution faircoin(0.5);
static ::std::uniform_real_distribution<SamplingFloatType> unif(0.0, 1.0);
#endif

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

typedef Eigen::SparseMatrix<unsigned, Eigen::ColMajor> AdjMatrix;

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

  PRNG prng(seed);

  x0.clear();
  x0.reserve(n);

  // generate ground truth
  for (unsigned i = 0; i < n; i++)
    x0.push_back(faircoin(prng) ? +1.0 : -1.0);

  const SamplingFloatType inp = a/static_cast<SamplingFloatType>(n);
  const SamplingFloatType outp = b/static_cast<SamplingFloatType>(n);

  for (unsigned src = 0; src < n; src++) {
    const int x0src = x0[src];
    for (unsigned dst = src; dst < n; dst++) {
      const SamplingFloatType r = unif(prng);
      const int x0ij = x0src * x0[dst];
      if (unlikely((x0ij == +1 && r <= inp) || (x0ij == -1 && r <= outp))) {
        entries.emplace_back(src, dst, 1);
        if (likely(src != dst))
          entries.emplace_back(dst, src, 1);
      }
    }
  }
}

static inline unsigned numNodes(const AdjMatrix &g)
{
  return g.rows();
}

static inline bool nonZeroDegree(const AdjMatrix &g, unsigned node)
{
  unsigned s = 0;
  for (AdjMatrix::InnerIterator it(g, node); it && s < 1; ++it, ++s)
    ;
  return s;
}

static inline unsigned degree(const AdjMatrix &g, unsigned node)
{
  unsigned s = 0;
  for (AdjMatrix::InnerIterator it(g, node); it; ++it, ++s)
    ;
  return s;
}

static inline double avgDegree(const AdjMatrix &g)
{
  unsigned s = 0;
  for (unsigned node = 0; node < numNodes(g); node++)
    s += degree(g, node);
  return static_cast<double>(s) / static_cast<double>(numNodes(g));
}

static unsigned numEdges(const AdjMatrix &g)
{
  unsigned s = 0;
  for (unsigned i = 0; i < g.outerSize(); i++)
    for (AdjMatrix::InnerIterator it(g, i); it; ++it)
      if (it.col() >= it.row())
        s++;
  return s;
}

static inline vector<int>
estimator(const Eigen::VectorXd &predictions)
{
  vector<int> ret(predictions.size());
  for (unsigned i = 0; i < ret.size(); i++)
    ret[i] = (predictions[i] >= 0.0) ? +1 : -1;
  return ret;
}

static pair<vector<int>, bool> solve(
    const AdjMatrix &g,
    unsigned r,
    unsigned seed,
    bool verbose,
    double eta = 0.0,
    unsigned epochs = 100000,
    double conv_tol = 1e-4)
{
  PRNG prng(seed);

  const unsigned n = numNodes(g);
  const unsigned nedges = numEdges(g);

  if (eta <= 0.0)
    eta = 2.0 * nedges / static_cast<double>(n * n);

  if (verbose)
    cout << "n " << n << ", r " << r << ", avgDegree " << avgDegree(g)
         << ", numDistinctEdges " << nedges << ", numTotalEdges " << g.nonZeros()
         << ", eta " << eta << ", epochs " << epochs << ", conv_tol " << conv_tol
         << endl;

  Eigen::MatrixXd spins(r, n); // Eigen dense matrices are optimized for column major
  for (unsigned node = 0; node < n; node++) {
    if (!nonZeroDegree(g, node))
      continue;
    for (unsigned i = 0; i < r; i++)
      spins(i, node) = faircoin(prng) ? +1.0 : -1.0;
    spins.col(node) /= spins.col(node).norm();
  }

  Eigen::VectorXd M = spins.rowwise().sum();
  Eigen::VectorXd cur, diff;

  vector<unsigned> pi;
  pi.reserve(n);
  for (unsigned i = 0; i < n; i++)
    pi.push_back(i);

  double maxdiff = 0;
  for (unsigned epoch = 0; epoch < epochs; epoch++) {
#ifdef USE_BOOST_RANDOM
    random_shuffle(pi.begin(), pi.end(), [&prng](int i) {
        ::boost::random::uniform_int_distribution<> unif(0, i-1);
        return unif(prng);
    });
#else
    shuffle(pi.begin(), pi.end(), prng);
#endif
    maxdiff = 0;
    for (auto node : pi) {
      cur = -(2.0 * eta) * M;
      for (AdjMatrix::InnerIterator it(g, node); it; ++it)
        cur += spins.col(it.index());
      const double curnorm = cur.norm();
      if (fabs(curnorm) <= 1e-10)
        continue;
      cur /= curnorm;
      diff = cur - spins.col(node);
      M += diff;
      maxdiff = std::max(maxdiff, diff.squaredNorm());
      spins.col(node) = cur;
    }
    if (verbose)
      cout << "epoch " << epoch << ", maxdiff " << maxdiff << endl;
    if (maxdiff <= conv_tol)
      break;
  }

  // compute covariance matrix
  const Eigen::MatrixXd cov = 1.0 / static_cast<double>(n) * (spins * spins.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(cov);

  // principle component
  const Eigen::VectorXd v1 = eigh.eigenvectors().col(r - 1);

  return std::make_pair(estimator(spins.transpose() * v1), maxdiff <= conv_tol);
}

static double errorRate(const vector<int> &truth, const vector<int> &predict)
{
  unsigned match0 = 0;
  unsigned match1 = 0;
  for (unsigned i = 0; i < truth.size(); i++) {
    if (truth[i] == predict[i])
      match0++;
    else
      match1++;
  }
  return std::min(match0, match1) / static_cast<double>(truth.size());
}

static double correlation(const vector<int> &truth, const vector<int> &predict)
{
  int s = 0;
  for (unsigned i = 0; i < truth.size(); i++)
    s += truth[i] * predict[i];
  return s / static_cast<double>(truth.size());
}

static unsigned numPositiveLabels(const vector<int> &truth)
{
  unsigned s = 0;
  for (auto label : truth)
    if (label == +1)
      s++;
  return s;
}

int main(int argc, char **argv)
{
	int verbose = 0;
  unsigned n = 0, r = 0, seed_gen = 0, seed_opt = 0;
  double a = 0, b = 0;
  double eta = 0.0;

  while (1) {
    static struct option long_options[] =
    {
      {"n"        , required_argument , 0        , 'n'} ,
      {"a"        , required_argument , 0        , 'a'} ,
      {"b"        , required_argument , 0        , 'b'} ,
      {"r"        , required_argument , 0        , 'r'} ,
      {"seed-gen" , required_argument , 0        , 'g'} ,
      {"seed-opt" , required_argument , 0        , 'o'} ,
      {"eta"      , required_argument , 0        , 't'} ,
      {"verbose"  , no_argument       , &verbose , 1}   ,
      {0          , 0                 , 0        , 0}
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "n:a:b:r:g:o:t:", long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
    case 0:
      if (long_options[option_index].flag != 0)
        break;
      abort();
      break;

    case 'n':
      n = strtoul(optarg, nullptr, 10);
      break;

    case 'a':
      a = strtod(optarg, nullptr);
      break;

    case 'b':
      b = strtod(optarg, nullptr);
      break;

    case 'r':
      r = strtoul(optarg, nullptr, 10);
      break;

    case 'g':
      seed_gen = strtoul(optarg, nullptr, 10);
      break;

    case 'o':
      seed_opt = strtoul(optarg, nullptr, 10);
      break;

    case 't':
      eta = strtod(optarg, nullptr);
      break;

    case '?':
      exit(1);

    default:
      abort();
    }
  }

  if (!n) {
    cerr << "n must be positive" << endl;
    exit(1);
  }

  if (!r) {
    cerr << "r must be positive" << endl;
    exit(1);
  }

  if (r >= n) {
    cerr << "warning: r >= n, setting r = n" << endl;
    r = n;
  }

  if (a < 0.0 || a > static_cast<double>(n)) {
    cerr << "a is invalid" << endl;
    exit(1);
  }

  if (b < 0.0 || b > static_cast<double>(n)) {
    cerr << "b is invalid" << endl;
    exit(1);
  }

  if (a <= b)
    cerr << "warning: a > b is needed for positive lambda" << endl;

  if (eta < 0.0) {
    cerr << "eta must be >= 0" << endl;
    exit(1);
  }

  vector<int> x0;
  vector<Eigen::Triplet<unsigned>> entries;

  if (verbose) {
    const double d = (a + b) / 2.0;
    const double lambda = (a - b) / sqrt(2.0 * (a + b));
    const double expectedTotalEdges = a + static_cast<double>(n-1) * d;
    cout << "graph summary: d " << d << ", lambda " << lambda
         << ", expectedTotalEdges " << expectedTotalEdges << endl;
  }

  timer t;
  sampleGraphWithGroundTruth(x0, entries, n, a, b, seed_gen);
  const double sampleMs = t.lap_ms();

  if (verbose) {
    cout << "graph sampling took " << sampleMs << "ms" << endl;
    cout << "ground truth: numPositiveLabels " << numPositiveLabels(x0) << endl;
  }

  AdjMatrix g(n, n);
  g.setFromTriplets(entries.begin(), entries.end());
  g.makeCompressed();

  t.lap();
  const auto pp = solve(g, r, seed_opt, verbose, eta);
  const auto predictions = pp.first;
  const auto converged = pp.second;
  const double solveMs = t.lap_ms();

  const double error = errorRate(x0, predictions);
  const double corr = correlation(x0, predictions);

  cout.precision(std::numeric_limits<double>::max_digits10);
  cout << "{"
       << "\"n\": " << n
       << ", \"r\": " << r
       << ", \"a\": " << a
       << ", \"b\": " << b
       << ", \"seed_gen\": " << seed_gen
       << ", \"seed_opt\": " << seed_opt
       << ", \"error\": " << error
       << ", \"correlation\": " << corr
       << ", \"converged\": " << (converged ? "true" : "false")
       << ", \"sample_ms\": " << sampleMs
       << ", \"solve_ms\": " << solveMs
       << "}" << endl;

  return 0;
}
