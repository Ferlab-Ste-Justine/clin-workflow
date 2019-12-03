package org.chusj

import java.time.LocalDate
import java.util
import java.util.{ArrayList, List}

import org.chusj.VepHelper._
import org.apache.spark.sql.{Column, DataFrame, Encoders, SparkSession}
import org.chusj.PatientHelper.loadPedigree
import org.apache.spark.sql.functions._
import org.chusj.VEPSparkDriverProgram.getMD5Hash
import org.json.JSONObject
import org.json.JSONArray

import scala.collection.mutable

//import scala.util.parsing.json.JSONArray


object ExomiserETL {

  case class Exomiser (
                        geneSymbol : String, geneId: String,
                        combinedScore: Double, variantQty: Int,
                        chrom: Seq[Long], position: Seq[Long],
                        ref: Seq[String], alt: Seq[String],
                        Inheritance: Seq[String],
                        qtyOfInheritance : Int
                      )


  def main(args: Array[String]): Unit = {


    val newArgs = if (args.length < 3)
       Array[String]("exomiser.json", "pedigree.properties", "pedigree.ped", "12", "50")
    else
      args
    val extractFile = newArgs(0)
    val pedigrePropsFile = newArgs(1)
    val pedFile = newArgs(2)
    val nbOfPartition = newArgs(3).toInt
    val bulkOpsQty = newArgs(4).toInt
    import collection.JavaConverters._
    val pedigreeProps = getPropertiesFromFile(pedigrePropsFile)
    val pedigrees: mutable.Seq[Pedigree] = loadPedigree(pedFile).asScala
    val patientMap: mutable.Map[String, Patient] = PatientHelper.preparePedigreeFromProps(pedigreeProps).asScala

    var totalCount = 0;

    pedigrees.foreach((ped:Pedigree) => println(ped))

    val proban = pedigrees(0)

    val build:String = pedigreeProps.get("assemblyVersion").toString

    for (elem <- patientMap) {
      println(s"ley=${elem._1},v=${elem._2}")
    }

    val spark = SparkSession.builder()
      .appName("Joins")
      .config("spark.master", "local")
      .getOrCreate()

    val exomiserGenesDF: DataFrame = spark.read
      .option("inferSchema", "false")
      .json("exomiser/FAM_C3_92.json")


//    val exomiserGenes2DF: DataFrame = spark.read
//      .option("inferSchema", "false")
//      .option("multiLine",true)
//      .json("exomiser.json")

//    exomiserGenes2DF.show()
//    exomiserGenes2DF.printSchema()

//    val selected2 = exomiserGenes2DF.select(
//      "geneSymbol",
//      "geneIdentifier.geneId",
//      "combinedScore",
//      "compatibleInheritanceModes",
//      "variantEvaluations.chromosome",
//      "variantEvaluations.position",
//      "variantEvaluations.ref",
//      "variantEvaluations.alt"
//    )

    val selected1 = exomiserGenesDF.select(
      col("geneSymbol"),
      col("geneIdentifier.geneId").as("geneId"),
      col("combinedScore"),
      size(col("variantEvaluations")).as("variantQty"),
      col("variantEvaluations.chromosome").as("chrom"),
      col("variantEvaluations.position").as("position"),
      col("variantEvaluations.ref").as("ref"),
      col("variantEvaluations.alt").as("alt"),
      col("variantEvaluations.compatibleInheritanceModes").as("Inheritance"),
      size(col("variantEvaluations.compatibleInheritanceModes")).as("qtyOfInheritance")
    )

    import spark.implicits._

    val exoDS = selected1.as[ExomiserETL.Exomiser]

    println(exoDS.count())

    val exoRDD = exoDS.rdd.repartition(nbOfPartition)

    exoRDD.foreachPartition(partitionOfRecords => {

    var jsonObjectList = new util.ArrayList[JSONObject]

    while (partitionOfRecords.hasNext) {
      val exomiser = partitionOfRecords.next()

      toJson(exomiser, build, jsonObjectList, proban.id)

      if (jsonObjectList.size() > bulkOpsQty) {
        VEPSparkDriverProgram.bulkStoreJsonObj(jsonObjectList, false, pedigreeProps, false, true)
        VEPSparkDriverProgram.TOTAL_COUNT += jsonObjectList.size()
        jsonObjectList = new util.ArrayList[JSONObject]

      }
    }
    // empty bucket
    VEPSparkDriverProgram.bulkStoreJsonObj(jsonObjectList, false, pedigreeProps, false, true)
    VEPSparkDriverProgram.TOTAL_COUNT += jsonObjectList.size()
    })

    println(s"proban is $proban")
    println(s"Total count=${VEPSparkDriverProgram.TOTAL_COUNT}")

  }


  def toJson(oneExo: Exomiser, build: String, objectList : util.ArrayList[JSONObject], specimenId: String): org.json.JSONArray = {
    val arrayOfVariant = new org.json.JSONArray
    //        JSONObject oneVariant = new JSONObject();
    var i = 0
    while ( {
      i < oneExo.variantQty
    }) {
      val oneVariant = new JSONObject

      val chrPos = oneExo.chrom(i)
      val position = oneExo.position(i)
      val reference = oneExo.ref(i)
      val alt = oneExo.alt(i)
      val localDate = LocalDate.now
      oneVariant.put("lastUpdate", localDate)
      oneVariant.put("chrom", chrPos)
      oneVariant.put("position", position)
      oneVariant.put("alt", alt )
      oneVariant.put("ref", reference)
      oneVariant.put("transmission", oneExo.Inheritance(i).replace("[", "").replace("]", ""))
      oneVariant.put("combinedScore", oneExo.combinedScore)
      oneVariant.put("specimenId", specimenId)
      arrayOfVariant.put(oneVariant)
      val dnaChange = reference + ">" + alt.split(",")(0)
      val mutationId = "chr" + chrPos + ":g." + position + dnaChange
      val uid = getMD5Hash(mutationId + "@" + build)
      oneVariant.put("id", uid)
      objectList.add(oneVariant)

      {
        i += 1; i - 1
      }
    }
    arrayOfVariant
  }


}
