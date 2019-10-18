#!/usr/bin/env bash
set -ex
time ~/bin/spark-2.4.3/bin/spark-submit --class org.chusj.VEPSparkDriverProgram --deploy-mode client --master 'local[*]' \
target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar out10000.txt pedigree.properties true local 'local[12]' 12g 12 60


#--spark-runner SPARK --spark-master local[2] --num-executors 2 --executor-cores 2 --executor-memory 2g -- out10000-3.txt out5000.txt