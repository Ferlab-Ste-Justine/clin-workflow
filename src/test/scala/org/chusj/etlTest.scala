package org.chusj

import java.util
import java.util.{ArrayList, Arrays, List}

import org.chusj.PatientHelper.loadPedigree
import org.scalatest.FunSuite

import collection.JavaConverters._
import scala.collection.mutable

class etlTest extends FunSuite {

  val genotypesFamily1: util.List[String] = new util.ArrayList[String](util.Arrays.asList("1/1", "0/0", "./."))
  val genotypesFamily2: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/1", "0/0", "0/0"))
  val genotypesFamily2b: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/1", "./.", "0/0"))
  val genotypesFamily3: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/0", "1/0", "0/0"))
  val genotypesFamily4: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/0", "0/0", "1/0"))
  val genotypesFamily5: util.List[String] = new util.ArrayList[String](util.Arrays.asList("1/0", "./.", "./."))
  val genotypesFamily6: util.List[String] = new util.ArrayList[String](util.Arrays.asList("1/1", "0/0", "0/0"))
  val genotypesFamily7: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/1", "1/0", "1/0"))
  val genotypesFamily8: util.List[String] = new util.ArrayList[String](util.Arrays.asList("1/1", "1/0", "1/0"))
  val genotypesFamily9: util.List[String] = new util.ArrayList[String](util.Arrays.asList("1/1", "0/1", "./."))
  val genotypesFamily10: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/1", "0/1", "0/0"))
  val genotypesFamily11: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/1", "0/1", "./."))
  val genotypesFamily12: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/1", "1/1", "0/1"))
  val genotypesFamily13: util.List[String] = new util.ArrayList[String](util.Arrays.asList("0/1", "0/0", "1/1"))
  val genotypesFamily14: util.List[String] = new util.ArrayList[String](util.Arrays.asList("1/1", "1/1", "0/1"))
  val genotypesFamily15: util.List[String] = new util.ArrayList[String](util.Arrays.asList("1/1", "./.", "./."))



  test("Testing zygosity") {

    assert(testZygo())

  }

  test("Testing count alternative allele") {

    assert(testCountAltAllele())

  }

  test("Testing count number of allele") {

    assert(testCountNumberOfAllele())

  }

  test("Testing impact score") {

    assert(
    VepHelper.getImpactScore("HIGH") == 4 &&
      VepHelper.getImpactScore("MODIFIER") == 1 &&
      VepHelper.getImpactScore("MODERATE") == 3 &&
      VepHelper.getImpactScore("LOW") == 2 &&
      VepHelper.getImpactScore("high") == 4 &&
      VepHelper.getImpactScore("BOB") == 0 &&
      VepHelper.getImpactScore("") == 0 &&
      VepHelper.getImpactScore(null) == 0)

  }

  test("Testing autosomal dominant transmission with a Solo") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest0.ped")
    // FA0002 SP00011 0 0 1 2

    assert(!testAutosomalDominant(genotypesFamily1, pedigrees))
    assert(testAutosomalDominant(genotypesFamily2, pedigrees))
    assert(testAutosomalDominant(genotypesFamily2b, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily3, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily4, pedigrees))
    assert(testAutosomalDominant(genotypesFamily5, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily6, pedigrees))
    assert(testAutosomalDominant(genotypesFamily7, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily8, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily9, pedigrees))
    assert(testAutosomalDominant(genotypesFamily10, pedigrees))
    assert(testAutosomalDominant(genotypesFamily11, pedigrees))
    assert(testAutosomalDominant(genotypesFamily12, pedigrees))
    assert(testAutosomalDominant(genotypesFamily13, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily14, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily15, pedigrees))

  }

  test("Testing autosomal dominant transmission with 0 affected parents") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest1.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 1
    FA0002 SP00036 0 0 2 1
     */

    assert(!testAutosomalDominant(genotypesFamily1, pedigrees))
    assert(testAutosomalDominant(genotypesFamily2, pedigrees))
    assert(testAutosomalDominant(genotypesFamily2b, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily3, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily4, pedigrees))
    assert(testAutosomalDominant(genotypesFamily5, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily6, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily7, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily8, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily9, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily10, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily11, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily12, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily13, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily14, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily15, pedigrees))

  }

  test("Testing autosomal dominant transmission with affected father @ pos 1") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest2.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 2
    FA0002 SP00036 0 0 2 1
     */

    assert(!testAutosomalDominant(genotypesFamily1, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily2, pedigrees))
    assert(testAutosomalDominant(genotypesFamily2b, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily3, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily4, pedigrees))
    assert(testAutosomalDominant(genotypesFamily5, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily6, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily7, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily8, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily9, pedigrees))
    assert(testAutosomalDominant(genotypesFamily10, pedigrees))
    assert(testAutosomalDominant(genotypesFamily11, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily12, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily13, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily14, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily15, pedigrees))


  }

  test("Testing autosomal dominant transmission with 2 affected parents") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest3.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 2
    FA0002 SP00036 0 0 2 HP:0000268
     */

    assert(!testAutosomalDominant(genotypesFamily1, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily2, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily2b, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily3, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily4, pedigrees))
    assert(testAutosomalDominant(genotypesFamily5, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily6, pedigrees))
    assert(testAutosomalDominant(genotypesFamily7, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily8, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily9, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily10, pedigrees))
    assert(testAutosomalDominant(genotypesFamily11, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily12, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily13, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily14, pedigrees))
    assert(!testAutosomalDominant(genotypesFamily15, pedigrees))

  }

  test("Testing autosomal recessif transmission with a solo") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest0.ped")
    // FA0002 SP00011 0 0 1 2

    assert(testAutosomalRecessif(genotypesFamily1, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2b, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily3, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily4, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily5, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily6, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily7, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily8, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily9, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily10, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily10, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily11, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily12, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily13, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily14, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily15, pedigrees))

  }

  test("Testing autosomal recessif transmission with 0 affected parents") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest1.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 1
    FA0002 SP00036 0 0 2 1
     */

    assert(!testAutosomalRecessif(genotypesFamily1, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2b, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily3, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily4, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily5, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily6, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily7, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily8, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily9, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily10, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily11, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily12, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily13, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily14, pedigrees))
    assert(testAutosomalRecessif(genotypesFamily15, pedigrees))


  }

  test("Testing autosomal recessif transmission with affected father @ pos 1") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest2.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 2
    FA0002 SP00036 0 0 2 1
     */

    assert(!testAutosomalRecessif(genotypesFamily1, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2b, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily3, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily4, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily5, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily6, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily7, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily8, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily9, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily10, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily11, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily12, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily13, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily14, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily15, pedigrees))


  }

  test("Testing autosomal recessif transmission with 2 affected parents") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest3.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 2
    FA0002 SP00036 0 0 2 HP:0000268
     */

    assert(!testAutosomalRecessif(genotypesFamily1, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily2b, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily3, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily4, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily5, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily6, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily7, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily8, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily9, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily10, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily11, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily12, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily13, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily14, pedigrees))
    assert(!testAutosomalRecessif(genotypesFamily15, pedigrees))
  }

  test("Testing denovo transmission with a Solo") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest0.ped")
    // FA0002 SP00011 0 0 1 2

    assert(testDenovo(genotypesFamily1, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily2, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily2b, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily3, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily4, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily5, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily6, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily7, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily8, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily9, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily10, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily11, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily12, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily13, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily14, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily15, pedigrees) == "NO")

  }

  test("Testing denovo transmission with 0 affected parents") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest1.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 1
    FA0002 SP00036 0 0 2 1
     */

    assert(testDenovo(genotypesFamily1, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily2, pedigrees) == "DENOVO")
    assert(testDenovo(genotypesFamily2b, pedigrees)== "PROBABLY DENOVO")
    assert(testDenovo(genotypesFamily3, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily4, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily5, pedigrees)== "PROBABLY DENOVO")
    assert(testDenovo(genotypesFamily6, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily7, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily8, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily9, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily10, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily11, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily12, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily13, pedigrees) == "DENOVO")
    assert(testDenovo(genotypesFamily14, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily15, pedigrees) == "NO")

  }

  test("Testing denovo transmission with affected father @ pos 1") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest2.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 2
    FA0002 SP00036 0 0 2 1
     */

    assert(testDenovo(genotypesFamily1, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily2, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily2b, pedigrees) == "PROBABLY DENOVO")
    assert(testDenovo(genotypesFamily3, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily4, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily5, pedigrees) == "PROBABLY DENOVO")
    assert(testDenovo(genotypesFamily6, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily7, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily8, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily9, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily10, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily11, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily12, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily13, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily14, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily15, pedigrees) == "NO")

  }

  test("Testing denovo transmission with 2 affected parents") {
    val pedigrees = loadPedigree("src/test/resources/pedigreeTest3.ped")
    /*
    FA0002 SP00011 SP00061 SP00036 1 2
    FA0002 SP00061 0 0 1 2
    FA0002 SP00036 0 0 2 HP:0000268
     */

    assert(testDenovo(genotypesFamily1, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily2, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily2b, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily3, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily4, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily5, pedigrees) == "PROBABLY DENOVO")
    assert(testDenovo(genotypesFamily6, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily7, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily8, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily9, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily10, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily11, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily12, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily13, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily14, pedigrees) == "NO")
    assert(testDenovo(genotypesFamily15, pedigrees) == "NO")

  }

  def testZygo(): Boolean = {

      VepHelper.zygosity("0/0") == "HOM REF" &&
      VepHelper.zygosity("0|0") == "HOM REF" &&
      VepHelper.zygosity("0/1") == "HET" &&
      VepHelper.zygosity("1/2") == "HET" &&
      VepHelper.zygosity("1/1") == "HOM" &&
      VepHelper.zygosity("2/2") == "HOM"

  }

  def testCountAltAllele() : Boolean = {

      VepHelper.countAlternativeAllele("0/0") == 0 &&
      VepHelper.countAlternativeAllele("0|0") == 0 &&
      VepHelper.countAlternativeAllele("0/1") == 1 &&
      VepHelper.countAlternativeAllele("1/0") == 1 &&
      VepHelper.countAlternativeAllele("1/1") == 2 &&
      VepHelper.countAlternativeAllele("1|1") == 2 &&
      VepHelper.countAlternativeAllele("./.") == 0

  }

  def testCountNumberOfAllele() : Boolean = {

    VepHelper.countNumberOfAllele("0/0") == 2 &&
      VepHelper.countNumberOfAllele("0|0") == 2 &&
      VepHelper.countNumberOfAllele("0/1") == 2 &&
      VepHelper.countNumberOfAllele("1/0") == 2 &&
      VepHelper.countNumberOfAllele("1/1") == 2 &&
      VepHelper.countNumberOfAllele("1|1") == 2 &&
      VepHelper.countNumberOfAllele("./.") == 0

  }


  def testAutosomalDominant(genotypesFamily: util.List[String], familyPed: util.List[Pedigree]) : Boolean = {
    VepHelper.isAutosomalDominant(genotypesFamily,familyPed)
  }

  def testAutosomalRecessif(genotypesFamily: util.List[String], familyPed: util.List[Pedigree]) : Boolean = {
    VepHelper.isAutosomalRecessif(genotypesFamily,familyPed)
  }

  def testDenovo(genotypesFamily: util.List[String], familyPed: util.List[Pedigree]) : String = {
    VepHelper.isDeNovo(genotypesFamily,familyPed)
  }


}
