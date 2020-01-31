#!/usr/bin/env bash
set -ex
time ~/bin/spark-2.4.3/bin/spark-submit --class org.chusj.VEPSparkDriverProgram --deploy-mode client --master 'local[*]' \
target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar batch.txt pedigree.properties true local 'local[12]' 12g 12 51 pedigreeTest1.ped


#--spark-runner SPARK --spark-master local[2] --num-executors 2 --executor-cores 2 --executor-memory 2g -- out10000-3.txt out5000.txt
#FAM_C3_92_new.txt pedigree.properties false local 'local[12]' 12g 12 51 pedigree.ped
#batch.txt pedigree.properties true local 'local[12]' 12g 12 51 pedigreeTest1.ped
