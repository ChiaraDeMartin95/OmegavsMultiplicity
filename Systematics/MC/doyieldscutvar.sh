#! /usr/bin/env bash

mkdir Eff

root -l -b -q effCalc.C\(\"results/MC-LHC22oapass4.root\",\"Eff/Eff-LHC22oapass4-Default.root\"\)

#systematics
for i in {200..500}
do
    root -l -b -q effCalc.C\(\"results/MC-LHC22oapass4-$i.root\",\"Eff/Eff-LHC22oapass4-$i.root\"\)
done