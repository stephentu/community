#!/bin/bash

FWDIR="$(cd `dirname $0`; pwd)"
pushd $FWDIR

spark-submit \
  --master yarn \
  --deploy-mode client \
  --num-executors 12 \
  --executor-cores 16 \
  --class edu.berkeley.cs.amplab.community.Driver \
  --driver-class-path $FWDIR/target/scala-2.10/community-assembly-0.1.0-SNAPSHOT.jar \
  --conf spark.executor.extraClassPath=$FWDIR/target/scala-2.10/community-assembly-0.1.0-SNAPSHOT.jar \
  --conf spark.serializer=org.apache.spark.serializer.JavaSerializer \
  --driver-memory 100g \
  --jars $FWDIR/target/scala-2.10/community-assembly-0.1.0-SNAPSHOT.jar \
  target/scala-2.10/community-assembly-0.1.0-SNAPSHOT.jar \
  --instance-binary /home/eecs/sltu/bin/instance-pipe.sh \
  --partition-factor 10
