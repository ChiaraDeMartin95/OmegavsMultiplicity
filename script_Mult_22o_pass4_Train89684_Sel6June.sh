#!/bin/bash

#MultType = 1 for FT0M, isMB = 0, mul = i
for i in `seq 0 11`
#for i in `seq 0 1`
    do
    echo "ciao"
        root -l -b -q "YieldsVsPt.C(1, \"\", 6, 0, 2, 1, 0, $i)" #omeganeg
        root -l -b -q "YieldsVsPt.C(1, \"\", 7, 0, 2, 1, 0, $i)" #omegaplus
        root -l -b -q "YieldsVsPt.C(1, \"\", 8, 0, 2, 1, 0, $i)" #omegasum
        root -l -b -q "YieldEffCorr.C(6, \"eff6June.root\", \"\", \"Yields\", \"LHC22o_pass4_Train89684\", 1, 1, 1, 0, $i)"
        root -l -b -q "YieldEffCorr.C(7, \"eff6June.root\", \"\", \"Yields\", \"LHC22o_pass4_Train89684\", 1, 1, 1, 0, $i)"
        root -l -b -q "YieldEffCorr.C(8, \"eff6June.root\", \"\", \"Yields\", \"LHC22o_pass4_Train89684\", 1, 1, 1, 0, $i)"
    done
#MB ones
root -l -b -q "YieldsVsPt.C(1, \"\", 6, 0, 2, 1, 1)"
root -l -b -q "YieldsVsPt.C(1, \"\", 7, 0, 2, 1, 1)" 
root -l -b -q "YieldsVsPt.C(1, \"\", 8, 0, 2, 1, 1)" 
root -l -b -q "YieldEffCorr.C(6, \"eff6June.root\", \"\", \"Yields\", \"LHC22o_pass4_Train89684\", 1, 1, 1, 1)"
root -l -b -q "YieldEffCorr.C(7, \"eff6June.root\", \"\", \"Yields\", \"LHC22o_pass4_Train89684\", 1, 1, 1, 1)"
root -l -b -q "YieldEffCorr.C(8, \"eff6June.root\", \"\", \"Yields\", \"LHC22o_pass4_Train89684\", 1, 1, 1, 1)"