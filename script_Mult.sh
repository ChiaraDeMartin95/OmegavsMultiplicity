#!/bin/bash

#MultType = 1 for FT0M, isMB = 0, mul = i
for i in `seq 0 9`
    do
        root -l -b -q "YieldsVsPt.C(1, \"_Sel6June\", 5, 2, 1, 0, $i)"
        root -l -b -q "YieldsVsPt.C(1, \"_Sel6June\", 6, 2, 1, 0, $i)" 
        root -l -b -q "YieldEffCorr.C(5, \"eff6June.root\", \"_Sel6June\", \"Yields\", \"LHC22m_pass4_Train79153\", 1, 1, 0, $i)"
        root -l -b -q "YieldEffCorr.C(6, \"eff6June.root\", \"_Sel6June\", \"Yields\", \"LHC22m_pass4_Train79153\", 1, 1, 0, $i)"
    done
#MB ones
root -l -b -q "YieldsVsPt.C(1, \"_Sel6June\", 5, 2, 1, 1)"
root -l -b -q "YieldsVsPt.C(1, \"_Sel6June\", 6, 2, 1, 1)" 
root -l -b -q "YieldEffCorr.C(5, \"eff6June.root\", \"_Sel6June\", \"Yields\", \"LHC22m_pass4_Train79153\", 1, 1, 1)"
root -l -b -q "YieldEffCorr.C(6, \"eff6June.root\", \"_Sel6June\", \"Yields\", \"LHC22m_pass4_Train79153\", 1, 1, 1)"