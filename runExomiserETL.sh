#!/usr/bin/env bash
set -ex
time ~/bin/spark-2.4.3/bin/spark-submit --class org.chusj.ExomiserETL --deploy-mode client --master 'local[*]' target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar exomiser.json


#--spark-runner SPARK --spark-master local[2] --num-executors 2 --executor-cores 2 --executor-memory 2g