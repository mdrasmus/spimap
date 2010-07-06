#!/bin/sh


rm -rf test/output/install
mkdir -p test/output/install
rm -rf build

make install prefix=test/output/install





