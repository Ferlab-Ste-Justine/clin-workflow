package org.chusj

import org.scalatest.FunSuite

class etlTest extends FunSuite  {

  test("testing zygosity") {

    assert(testZygo() == true)

  }


  def testZygo(): Boolean = {

    VepHelper.zygosity("0/0") == "HOM REF" &&
      VepHelper.zygosity("0/1") == "HET REF" &&
      VepHelper.zygosity("1/2") == "HET ALT" &&
      VepHelper.zygosity("1/1") == "HOM ALT" &&
      VepHelper.zygosity("2/2") == "HOM ALT"

  }

}
