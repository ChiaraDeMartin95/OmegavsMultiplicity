#! /usr/bin/env bash

mkdir Yields

#default
root -l -b -q YieldsVsPt.C\(\)


#systematics
for i in {0..500}
do
    root -l -b -q YieldsVsPt.C\(\"$i\",\"results/Data-LHC22oapass4-$i.root\"\)
done

