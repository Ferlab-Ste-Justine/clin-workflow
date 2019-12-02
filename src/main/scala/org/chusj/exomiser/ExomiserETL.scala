package org.chusj.exomiser

import java.util
import java.util.{List, Map, Properties}

import org.chusj.VepHelper._
import org.apache.spark.sql.{Column, DataFrame, SparkSession}
import org.chusj.PatientHelper.loadPedigree
import org.chusj.{Patient, PatientHelper, Pedigree}
import org.apache.spark.sql.functions._


object ExomiserETL {

  def main(args: Array[String]): Unit = {

    val newArgs = if (args.length < 3)
       Array[String]("exomiser.json", "pedigree.properties", "pedigree.ped")
    else
      args
    val extractFile = newArgs(0)
    val pedigrePropsFile = newArgs(1)
    val pedFile = newArgs(2)

    val pedigreeProps = getPropertiesFromFile(pedigrePropsFile)
    val pedigrees = loadPedigree(pedFile)
    val patientMap = PatientHelper.preparePedigreeFromProps(pedigreeProps)

    pedigrees.forEach((ped : Pedigree) => println (ped))

    patientMap.forEach((k,v) => println(s"k=$k, v=$v"))

    case class Exomiser (combinedScore: Double,
                         compatibleInheritanceModes: Array[String],
                         geneIdentifier: String, geneScore:String, geneSymbol:String, priorityRes:String, priorityScore: String, variant:String, variantScore : String)

    val spark = SparkSession.builder()
      .appName("Joins")
      .config("spark.master", "local")
      .getOrCreate()

    val exomiserGenesDF: Unit = spark.read
      .option("inferSchema", "false")
      .json("exomiser/FAM_C3_92.json")
      .show()

    val exomiserGenes2DF: DataFrame = spark.read
      .option("inferSchema", "false")
      .option("multiLine",true)
      .json("exomiser.json")

    exomiserGenes2DF.show()
    exomiserGenes2DF.printSchema()


    import org.apache.spark.sql.DataFrame
    import org.apache.spark.sql.types.{DataType, StructType}

    implicit class DataFrameFlattener(df: DataFrame) {
      def flattenSchema: DataFrame = {
        df.select(flatten(Nil, df.schema): _*)
      }

      protected def flatten(path: Seq[String], schema: DataType): Seq[Column] = schema match {
        case s: StructType => s.fields.flatMap(f => flatten(path :+ f.name, f.dataType))
        //case ar: Array[] => ar.fields.DataFrameFlattener.this
        case other => col(path.map(n => s"`$n`").mkString(".")).as(path.mkString(".")) :: Nil
      }
    }

    exomiserGenes2DF.flattenSchema.flattenSchema.flattenSchema.show()
    exomiserGenes2DF.flattenSchema.show()

    val flattenDF = exomiserGenes2DF.flattenSchema
    flattenDF.show()
    flattenDF.printSchema()
    import spark.implicits._

    //val exomiserDS = exomiserGenes2DF.as[Exomiser]

    //exomiserDS.show()


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
