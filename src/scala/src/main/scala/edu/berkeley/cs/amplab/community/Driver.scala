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
    sc.addFile(config.instanceBinPath)

    val rdd = sc.parallelize(Grid.grid(16000, 5.0, Grid.DefaultLambdas, 1))
    val result = rdd.pipe(SparkFiles.get(config.instanceBinPath)).collect()

    println(result)

  }

}
