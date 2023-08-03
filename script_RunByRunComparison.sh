#!/bin/bash

SPathIn="RunByRunComparison/AnalysisResults"
SPathInEvt="RunByRunComparison/AnalysisResultsEvts"

#old ones (LHC22o_pass4)
#year="LHC22o_pass4_Train99659"
#RunN=( 528534 528531 528463 528461 528448 528381 528379 528232 )
#CascPPTrain="_Train100720"

#new ones (LHC22o_pass4_MinBias) 
year="LHC22o_pass4_MinBias_Train108123"
RunN=( 528531 528461 528292 527899 527895 527871 527850 527240 527109 527057 527041 526964 526641 )
#CascPPTrain="_Train109354" #omega
CascPPTrain="_Train109827" #xi

#for i in `seq 0 7` #8 runs (LHC22o_pass4)
for i in `seq 0 12` #12 runs (LHC22o_pass4_MinBias)
do
Sel=${CascPPTrain}"_Run"${RunN[i]}
StringPathInEvt="_Run"${RunN[i]}
echo ${Sel}
echo ${StringPathInEvt}
    #root -l -b -q "YieldsVsPt.C(1, \"${StringPathInEvt}\", \"${Sel}\", 8, 0, 2, 1, 1, 0, \"${year}\", \"${SPathIn}\", \"${SPathInEvt}\")" #omegasum
    root -l -b -q "YieldsVsPt.C(1, \"${Sel}\", \"${StringPathInEvt}\", 5, 0, 2, 1, 1, 0, \"${year}\", \"${SPathIn}\", \"${SPathInEvt}\")" #omegasum
done
