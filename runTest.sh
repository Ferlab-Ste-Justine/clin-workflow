#!/usr/bin/env bash
set -ex
time ~/bin/spark-2.4.3/bin/spark-submit --class org.chusj.test.Test --deploy-mode client --master 'local[*]' target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar out11180_S12_L008avtdf.txt 11180 S12 p123 chujs WXS P Y local 'local[24]' 12g 24


#--spark-runner SPARK --spark-master local[2] --num-executors 2 --executor-cores 2 --executor-memory 2g