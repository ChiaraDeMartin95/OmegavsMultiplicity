#!/bin/bash

#MultType = 1 for FT0M, isMB = 0, mul = i
for i in `seq 0 11`
    do
        #root -l -b -q "YieldsVsPt.C(3, 0, $i)" #particle type
        #root -l -b -q "YieldsVsPt.C(4, 0, $i)"
        #root -l -b -q "YieldsVsPt.C(5, 0, $i)"
        #root -l -b -q "YieldEffCorr.C(3, 0, $i)"
        #root -l -b -q "YieldEffCorr.C(3, 0, $i)"
        #root -l -b -q "YieldEffCorr.C(3, 0, $i)"
        root -l -b -q "YieldsVsPt.C(6, 0, $i)" #particle type
        root -l -b -q "YieldsVsPt.C(7, 0, $i)"
        root -l -b -q "YieldsVsPt.C(8, 0, $i)"
        root -l -b -q "YieldEffCorr.C(6, 0, $i)"
        root -l -b -q "YieldEffCorr.C(7, 0, $i)"
        root -l -b -q "YieldEffCorr.C(8, 0, $i)"
    done
#MB ones
#root -l -b -q "YieldsVsPt.C(3, 1)"
#root -l -b -q "YieldsVsPt.C(4, 1)"
#root -l -b -q "YieldsVsPt.C(5, 1)"
#root -l -b -q "YieldEffCorr.C(3, 1)"
#root -l -b -q "YieldEffCorr.C(4, 1)"
#root -l -b -q "YieldEffCorr.C(5, 1)"
root -l -b -q "YieldsVsPt.C(6, 1)"
root -l -b -q "YieldsVsPt.C(7, 1)"
root -l -b -q "YieldsVsPt.C(8, 1)"
root -l -b -q "YieldEffCorr.C(6, 1)"
root -l -b -q "YieldEffCorr.C(7, 1)"
root -l -b -q "YieldEffCorr.C(8, 1)"