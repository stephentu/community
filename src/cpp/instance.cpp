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
    double eta = 0.0,
    unsigned epochs = 10,
    double conv_tol = 1e-4)
{
  PRNG prng(seed);
  bernoulli_distribution faircoin(0.5);

  const unsigned n = numNodes(g);

  if (eta <= 0.0)
    eta = 2.0 * numEdges(g) / static_cast<double>(n * n);

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
    cout << "epoch " << epoch << ", maxdiff " << maxdiff << endl;
    if (maxdiff <= conv_tol)
      break;
  }

  // compute covariance matrix
  const Eigen::MatrixXd cov = 1.0 / static_cast<double>(n) * (spins * spins.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(cov);

  cout << "evals " << eigh.eigenvalues() << endl;

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

  const unsigned n = 128000;
  const double a = 20.0;
  const double b = 3.0;
  const unsigned r = 100;

  //const unsigned n = 5;
  //const double a = 3;
  //const double b = 2;
  //const unsigned r = 2;

  vector<int> x0;
  vector<Eigen::Triplet<unsigned>> entries;
  {
    scoped_timer t("sample graph");
    sampleGraphWithGroundTruth(x0, entries, n, a, b, 8323753);
  }

  AdjMatrix g(n, n);
  g.setFromTriplets(entries.begin(), entries.end());
  g.makeCompressed();

  //cout << "nonZeroDegree(0) " << nonZeroDegree(g, 0) << endl;
  //cout << "nonZeroDegree(1) " << nonZeroDegree(g, 1) << endl;
  //cout << "nonZeroDegree(2) " << nonZeroDegree(g, 2) << endl;
  //cout << "numEdges() " << numEdges(g) << endl;
  //cout << Eigen::MatrixXd(g) << endl;
  //cout << "degree(0) " << degree(g, 0) << endl;
  //cout << "degree(1) " << degree(g, 1) << endl;
  //cout << "degree(2) " << degree(g, 2) << endl;
  //cout << "degree(3) " << degree(g, 3) << endl;
  //cout << "degree(4) " << degree(g, 4) << endl;

  auto predictions = solve(g, r, 39243, 0.0, 1000);

  cout << "error " << errorRate(x0, predictions) << endl;

  return 0;
}
