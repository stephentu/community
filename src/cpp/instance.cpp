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
#include <random>
#include <vector>
#include <exception>

#include <getopt.h>
#include <unistd.h>

#include "timer.hpp"

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

typedef Eigen::SparseMatrix<unsigned, Eigen::ColMajor> AdjMatrix;

typedef std::minstd_rand PRNG;

using namespace std;

static bernoulli_distribution faircoin(0.5);
static uniform_real_distribution<float> unif(0.0, 1.0);

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

  const float inp = a/static_cast<float>(n);
  const float outp = b/static_cast<float>(n);

  for (unsigned src = 0; src < n; src++) {
    const int x0src = x0[src];
    for (unsigned dst = src; dst < n; dst++) {
      const float r = unif(prng);
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

static vector<int> solve(
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

  if (eta <= 0.0)
    eta = 2.0 * numEdges(g) / static_cast<double>(n * n);

  if (verbose)
    cout << "n " << n << ", r " << r << ", avgDegree " << avgDegree(g)
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

  for (unsigned epoch = 0; epoch < epochs; epoch++) {
    double maxdiff = 0;
    shuffle(pi.begin(), pi.end(), prng);
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

  //cout << "evals " << eigh.eigenvalues() << endl;

  // principle component
  const Eigen::VectorXd v1 = eigh.eigenvectors().col(r - 1);

  return estimator(spins.transpose() * v1);
}

static double errorRate(const vector<int> &truth, const vector<int> &predict)
{
  unsigned match0 = 0;
  unsigned match1 = 1;
  for (unsigned i = 0; i < truth.size(); i++) {
    if (truth[i] == predict[i])
      match0++;
    else
      match1++;
  }
  return std::min(match0, match1) / static_cast<double>(truth.size());
}

int main(int argc, char **argv)
{

	int verbose = 0;
  unsigned n = 0, r = 0, seed_gen = 0, seed_opt = 0;
  double a = 0, b = 0;

  while (1) {
    static struct option long_options[] =
    {
      {"n"        , required_argument , 0        , 'n'} ,
      {"a"        , required_argument , 0        , 'a'} ,
      {"b"        , required_argument , 0        , 'b'} ,
      {"r"        , required_argument , 0        , 'r'} ,
      {"seed-gen" , required_argument , 0        , 'g'} ,
      {"seed-opt" , required_argument , 0        , 'o'} ,
      {"verbose"  , no_argument       , &verbose , 1}   ,
      {0, 0, 0, 0}
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "n:a:b:r:g:o:", long_options, &option_index);
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

  //const unsigned n = 128000;
  //const double a = 20.0;
  //const double b = 3.0;
  //const unsigned r = 100;

  //const unsigned n = 5;
  //const double a = 3;
  //const double b = 2;
  //const unsigned r = 2;

  vector<int> x0;
  vector<Eigen::Triplet<unsigned>> entries;
  {
    scoped_timer t("sample graph", verbose);
    sampleGraphWithGroundTruth(x0, entries, n, a, b, seed_gen);
  }

  if (verbose)
    cout << "graph summary: d " << (a + b)/2.0
         << ", lambda " << ((a - b) / sqrt(2.0 * (a + b)))
         << endl;

  AdjMatrix g(n, n);
  g.setFromTriplets(entries.begin(), entries.end());
  g.makeCompressed();

  auto predictions = solve(g, r, seed_opt, verbose);

  const double error = errorRate(x0, predictions);

  cout << "error " << error << endl;

  return 0;
}
