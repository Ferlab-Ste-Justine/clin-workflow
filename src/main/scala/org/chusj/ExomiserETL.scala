package org.chusj

import java.time.LocalDate
import java.util
import java.util.{ArrayList, List}

import org.apache.http.HttpHost
import org.chusj.VepHelper._
import org.apache.spark.sql.{Column, DataFrame, Encoders, SparkSession}
import org.chusj.PatientHelper.loadPedigree
import org.apache.spark.sql.functions._
import org.chusj.VEPSparkDriverProgram.getMD5Hash
import org.elasticsearch.client.{RestClient, RestHighLevelClient}
import org.json.JSONObject

import scala.collection.mutable


object ExomiserETL {

  case class Exomiser (
                        geneSymbol : String, geneId: String,
                        combinedScore: Float, variantQty: Int,
                        chrom: Long, position: Long,
                        ref: String, alt: String,
                        Inheritance: Array[String],
                        qtyOfInheritance : Int
                      )


  def main(args: Array[String]): Unit = {


    val newArgs = if (args.length < 3)
       Array[String]("exomiser/FAM_C3_92.json", "pedigree.properties", "pedigree.ped", "6", "20")
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
      .json(extractFile)

    import spark.implicits._

    val exomiserGenesDFTransformed = exomiserGenesDF.select(
//      col("geneSymbol"),
//      col("geneIdentifier.geneId").as("geneId"),
//      col("combinedScore"),
//      size(col("variantEvaluations")).as("variantQty"),
//      col("variantEvaluations.chromosome").as("chrom"),
//      col("variantEvaluations.position").as("position"),
//      col("variantEvaluations.ref").as("ref"),
//      col("variantEvaluations.alt").as("alt"),
//      col("variantEvaluations.compatibleInheritanceModes").as("Inheritance"),
//      size(col("variantEvaluations.compatibleInheritanceModes")).as("qtyOfInheritance")
//    ).select(
      col("geneSymbol"),
      col("geneIdentifier.geneId").as("geneId"),
      col("combinedScore"),
      size(col("variantEvaluations")).as("variantQty"),
      explode(col("variantEvaluations")).as("variant"),
      size(col("variantEvaluations.compatibleInheritanceModes")).as("qtyOfInheritance")
    ).select(
      col("geneSymbol"),
      col("geneId"),
      expr("round (cast(combinedScore as float),5) combinedScore"),
      col("variantQty"),
      col("variant.chromosome").as("chrom"),
      col("variant.ref").as("ref"),
      col("variant.alt").as("alt"),
      col("variant.position").as("position"),
      col("variant.compatibleInheritanceModes").as("Inheritance"),
      col("qtyOfInheritance")
    )


//    selected3.show()
    println(exomiserGenesDFTransformed.count())

    val exoDS = exomiserGenesDFTransformed.as[ExomiserETL.Exomiser]
//    val filtered = exoDS.filter((variant) => {
//      variant.qtyOfInheritance != variant.variantQty
//    })
//    filtered.show()
    exoDS.show()
    println(exoDS.count())

    val clientTry = new RestHighLevelClient(
      RestClient.builder(
        new HttpHost("localhost", 9200, "http")))

    VEPSparkDriverProgram.client = clientTry

    val exoRDD = exoDS.rdd.repartition(nbOfPartition)
    exoRDD.foreachPartition(partitionOfRecords => {

      // bucket
      var jsonObjectList = new util.ArrayList[JSONObject]

      while (partitionOfRecords.hasNext) {
        val exomiser = partitionOfRecords.next()

        jsonObjectList.add(toJsonObj(exomiser,build,proban.id))

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
    VEPSparkDriverProgram.client.close()
  }


  def toJsonObj(oneExo: Exomiser, build: String, specimenId: String) = {
    val oneVariant = new JSONObject

    // Exomiser chromose X is 23 while VEP and previous ETL is using X
    val chrom = if (oneExo.chrom == 23) "X" else oneExo.chrom
    oneVariant.put("lastUpdate", LocalDate.now)
    oneVariant.put("chrom", chrom)
    oneVariant.put("position", oneExo.position)
    oneVariant.put("alt", oneExo.alt )
    oneVariant.put("ref", oneExo.ref)
    oneVariant.put("transmission", oneExo.Inheritance)
    oneVariant.put("combinedScore", oneExo.combinedScore)
    oneVariant.put("specimenId", specimenId)
    //arrayOfVariant.put(oneVariant)
    val mutationId = "chr" + chrom + ":g." + oneExo.position + oneExo.ref + ">" + oneExo.alt.split(",")(0)
    oneVariant.put("mutationId", mutationId)
    val uid = getMD5Hash(mutationId + "@" + build)
    oneVariant.put("id", uid)

    oneVariant
  }


}
