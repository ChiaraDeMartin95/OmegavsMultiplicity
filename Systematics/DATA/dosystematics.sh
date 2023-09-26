#! /usr/bin/env bash

mkdir results

#default
o2-analysis-lf-cascpostprocessing -b --configuration json://config/configDEF.json
mv AnalysisResults.root results/Data-LHC22oapass4.root


#systematics
for i in {0..500}
do
    o2-analysis-lf-cascpostprocessing -b  --configuration json://config/systjson$i.json
    mv AnalysisResults.root results/Data-LHC22oapass4-$i.root
done
