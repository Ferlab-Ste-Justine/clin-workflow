#!/usr/bin/env bash
set -ex
time ~/bin/spark-2.4.3/bin/spark-submit --class org.chusj.VEPSparkDriverProgram --deploy-mode client --master 'local[*]' \
target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar FAM_C3_92_new.txt etl.properties true local 'local[12]' 12g 12 51 pedigree.ped 9200


#--spark-runner SPARK --spark-master local[2] --num-executors 2 --executor-cores 2 --executor-memory 2g -- out10000-3.txt out5000.txt
#FAM_C3_92_new.txt etl.properties true local 'local[12]' 12g 12 51 pedigree.ped 9201
#batch.txt etl.properties true local 'local[12]' 12g 12 51 pedigreeTest1.ped 9201
