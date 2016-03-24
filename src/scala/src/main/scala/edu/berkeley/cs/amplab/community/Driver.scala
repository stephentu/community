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

    val grid = Grid.grid(16000, 5.0, Grid.DefaultLambdas, 1).map(_.toCLI(100))
    println("Grid length: " + grid.size)
    val rdd = sc.parallelize(grid)
    val result = rdd.pipe(config.instanceBinPath)
    //val result = rdd.pipe("cat")
    result.cache()
    println("result.count: " + result.count())
    result.saveAsTextFile("/community-results")    

  }

}
