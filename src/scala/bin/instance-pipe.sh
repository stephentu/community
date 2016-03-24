#!/bin/bash
while read args ; do
    if [ -n "$args" ]; then
        /home/eecs/sltu/bin/instance $args
    fi
done
