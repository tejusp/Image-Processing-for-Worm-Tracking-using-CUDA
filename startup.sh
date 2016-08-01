#!/bin/sh
make clean
make -j 48
cat main.cpp | grep "frames("
date
./main
date
