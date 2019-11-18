package org.chusj.test

import org.apache.spark.sql.SparkSession


object Test extends {

  def main(args: Array[String]): Unit = {
    val spark = SparkSession.builder()
      .appName("Joins")
      .config("spark.master", "local")
      .getOrCreate()


    //  val testing = new test()
    //  println(testing.parse("4"))


    val mutationsDF: Unit = spark.read
      .option("inferSchema", "true")
      .json("mutations.json")
      .show()
  }


//  class test {
//
//    case class Failure(reason: String)
//
//    def parse(numberAsString: String) : Either[Failure,Int] =
//      try {
//        val result = Integer.parseInt(numberAsString)
//        Right(result)
//      } catch {
//        case _ : Throwable => Left(Failure("Error when parsing number")) }
//
//  }





}
