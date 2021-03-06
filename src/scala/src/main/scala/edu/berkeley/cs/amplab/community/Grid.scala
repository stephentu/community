package edu.berkeley.cs.amplab.community

import collection.mutable.HashSet

case class GridPoint(n: Int, trial: Int, a: Double, b: Double, seedSample: Long, seedOpt: Long) {
  def toCLI(r: Int): String = 
    Seq("--n", n.toString,
        "--r", r.toString,
        "--a", a.toString,
        "--b", b.toString,
        "--seed-gen", seedSample.toString,
        "--seed-opt", seedOpt.toString).mkString(" ")
}

object Grid {

  val DefaultLambdas = Seq(1.0, 1.005, 1.01, 1.015, 1.02, 1.025, 1.03)

  private def asUnsignedLong(x: Int): Long = 
    x & 0x00000000ffffffffL

  private[this] class LongSamplerWithoutReplacement {
    private val seen = new HashSet[Long]
    def next(): Long = {
      var candidate = asUnsignedLong(util.Random.nextInt) 
      while (seen.contains(candidate)) {
        candidate = asUnsignedLong(util.Random.nextInt)
      }
      seen += candidate
      candidate
    }
  }

  def parameters(n: Int, d: Double, lam: Double, tol: Double = 1e-8): (Double, Double) = {
    val b = d - lam * Math.sqrt(d)
    val a = 2.0 * lam * Math.sqrt(d) + b

    if (a < 0.0 || a > n.toDouble)
      throw new RuntimeException("a is not in range")
    if (b < 0.0 || b > n.toDouble)
      throw new RuntimeException("b is not in range")

    if (Math.abs((a + b)/2.0 - d) > tol)
      throw new RuntimeException("does not reconstruct d")
    if (Math.abs((a - b)/Math.sqrt(2.0 * (a + b)) - lam) > tol)
      throw new RuntimeException("does not reconstruct lambda")

    (a, b)
  }

  def ntrials(n: Int): Int = Math.ceil(6.4e8 / n.toDouble).toInt

  def grid(n: Int, d: Double, lambdas: Seq[Double], numTrials: Int): Seq[GridPoint] = {
    val seedGenSamp = new LongSamplerWithoutReplacement
    val seedOptSamp = new LongSamplerWithoutReplacement

    lambdas.flatMap { lam =>
      val (a, b) = parameters(n, d, lam)
      (0 until numTrials).map { trial =>
        GridPoint(n, trial, a, b, seedGenSamp.next(), seedOptSamp.next())
      }
    }
  }

  def grid(n: Int, d: Double, lambdas: Seq[Double]): Seq[GridPoint] =
    grid(n, d, lambdas, ntrials(n))

}
