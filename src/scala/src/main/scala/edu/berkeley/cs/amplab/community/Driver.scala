package edu.berkeley.cs.amplab.community

import org.apache.spark.rdd.RDD
import org.apache.spark.broadcast.Broadcast
import org.apache.spark.SparkConf
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkFiles

object Driver {

  case class Config(
    instanceBinPath: String = "")

  def main(args: Array[String]) {
    val parser: scopt.OptionParser[Config] = new scopt.OptionParser[Config]("Driver") {
          head("Driver")
          opt[String]("instance-binary") required() action { (x, c) => c.copy(instanceBinPath=x) }
    }
    println(args)
    val config = parser.parse(args, Config()).getOrElse(sys.exit(1))

    val conf = new SparkConf().setAppName("Driver")

    val sc = new SparkContext(conf)
    sc.setCheckpointDir("/tmp/spark-checkpoint")

    val r = 100
    val d = 10.0
    val hdfsFileName = s"/community-results-d${d}"
    println(s"Driver r=${r}, d=${d}")
    println(s"Saving to hdfsFileName=${hdfsFileName}")

    val grid =
      Grid.grid(16000, d, Grid.DefaultLambdas) ++
      Grid.grid(32000, d, Grid.DefaultLambdas)
    println("Grid length: " + grid.size)

    val rdd = sc.parallelize(grid.map(_.toCLI(r)))
    val result = rdd.pipe(config.instanceBinPath)
    result.cache()
    println("result.count: " + result.count())

    result.saveAsTextFile(hdfsFileName)
  }

}
