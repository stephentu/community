#!/bin/bash
read args
if [ -n "$args" ]; then
    /home/eecs/sltu/bin/instance $args
fi
