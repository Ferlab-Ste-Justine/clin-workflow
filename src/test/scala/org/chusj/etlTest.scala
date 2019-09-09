package org.chusj

import org.scalatest.FunSuite

class etlTest extends FunSuite  {

  test("Testing zygosity") {

    assert(testZygo() == true)

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


  def testZygo(): Boolean = {

    VepHelper.zygosity("0/0") == "HOM REF" &&
      VepHelper.zygosity("0/1") == "HET REF" &&
      VepHelper.zygosity("1/2") == "HET ALT" &&
      VepHelper.zygosity("1/1") == "HOM ALT" &&
      VepHelper.zygosity("2/2") == "HOM ALT"

  }

}
