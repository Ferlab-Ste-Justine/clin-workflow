#!/usr/bin/env bash
set -ex
time ~/bin/spark-2.4.3/bin/spark-submit --class org.chusj.VEPSparkDriverProgram --deploy-mode client --master 'local[*]' target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar out10000-3.txt 11180 S12 p123 chujs WXS P Y local 'local[8]' 10g 8


#--spark-runner SPARK --spark-master local[2] --num-executors 2 --executor-cores 2 --executor-memory 2g