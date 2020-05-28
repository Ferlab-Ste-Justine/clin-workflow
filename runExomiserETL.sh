#!/usr/bin/env bash
set -ex
time ~/bin/spark-2.4.3/bin/spark-submit --class org.chusj.ExomiserETL --deploy-mode client --master 'local[*]' \
target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar exomiser/FAM_C3_92.json etl.properties SP00011 6 45 9201


#--spark-runner SPARK --spark-master local[2] --num-executo rs 2 --executor-cores 2 --executor-memory 2g
