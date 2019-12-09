#!/usr/bin/env bash
set -ex
time java -jar target/ExtractTLoad-1.0-SNAPSHOT-jar-with-dependencies.jar ${1} ${2} ${3} validate

# without redis
# real    0m4.124s
# With Redis calls
# real    0m4.951s
