#! /usr/bin/env bash

mkdir CorrSpectra

#default
root -l -b -q YieldEffCorr.C\(\)

#systematics
for i in {0..500}
do
    root -l -b -q YieldEffCorr.C\(\"$i\"\)
done