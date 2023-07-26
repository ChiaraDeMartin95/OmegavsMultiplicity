#!/bin/bash

Sel="_Train100720"
year="LHC22o_pass4_Train99659"
SPathIn="RunByRunComparison/AnalysisResults"
SPathInEvt="RunByRunComparison/AnalysisResultsEvts"

for i in `seq 0 7` #8 runs
do
    if [ $i ==  0 ]
    then Sel="_Train100720_Run528534" StringPathInEvt="_Run528534"
    fi
    if [ $i ==  1 ]
    then Sel="_Train100720_Run528531" StringPathInEvt="_Run528531"
    fi
    if [ $i ==  2 ]
    then Sel="_Train100720_Run528463" StringPathInEvt="_Run528463"
    fi
    if [ $i ==  3 ]
    then Sel="_Train100720_Run528461" StringPathInEvt="_Run528461"
    fi
    if [ $i ==  4 ]
    then Sel="_Train100720_Run528448" StringPathInEvt="_Run528448"
    fi
    if [ $i ==  5 ]
    then Sel="_Train100720_Run528381" StringPathInEvt="_Run528381"
    fi
    if [ $i ==  6 ]
    then Sel="_Train100720_Run528379" StringPathInEvt="_Run528379"
    fi
    if [ $i ==  7 ]
    then Sel="_Train100720_Run528232" StringPathInEvt="_Run528232"
    fi
        root -l -b -q "YieldsVsPt.C(1, \"${StringPathInEvt}\", \"${Sel}\", 8, 0, 2, 1, 1, 0, \"${year}\", \"${SPathIn}\", \"${SPathInEvt}\")" #omegasum
done
