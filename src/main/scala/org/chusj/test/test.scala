package org.chusj.test



object test extends App {


  class test {

    case class Failure(reason: String)

    def parse(numberAsString: String) : Either[Failure,Int] =
      try {
        val result = Integer.parseInt(numberAsString)
        Right(result)
      } catch {
        case _ : Throwable => Left(Failure("Error when parsing number")) }

  }


  val testing = new test()
  println(testing.parse("4"))

}
