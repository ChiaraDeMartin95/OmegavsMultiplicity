#!/bin/bash

#casccospa
#for i in 0 1 2 3 4 5 6 7
for i in 0 1 2 3 4 5
#for i in 0
do
#root -l -b -q "YieldsVsPt.C(1, \"_casccospa$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_dcacascdau$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_dcabachtopv$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_dcapostopv$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_dcanegtopv$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_lambdamasswin$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_rejcomp$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_nsigmatpcKa$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_cascradius$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_v0radius$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_dcav0dau$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_v0cospa$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_casclifetime$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_cosbachbar$i\")"
#root -l -b -q "YieldsVsPt.C(1, \"_dcabachbar$i\")"
root -l -b -q "YieldsVsPt.C(1, \"_dcav0topv$i\")"
done
for i in 0 1 2 3 4
do
#root -l -b -q "CompareYields.C(2, $i)" #bachtopv
#root -l -b -q "CompareYields.C(3, $i)" #postopv
#root -l -b -q "CompareYields.C(4, $i)" #negtopv
#root -l -b -q "CompareYields.C(5, $i)" #lambdamass
#root -l -b -q "CompareYields.C(6, $i)" #rejcomp
#root -l -b -q "CompareYields.C(7, $i)" #nsigmaTPCKa
#root -l -b -q "CompareYields.C(8, $i)" #cascradius
#root -l -b -q "CompareYields.C(9, $i)" #v0radius
#root -l -b -q "CompareYields.C(10, $i)" #dcav0dau
#root -l -b -q "CompareYields.C(11, $i)" #v0cospa
#root -l -b -q "CompareYields.C(12, $i)" #casclifetime
#root -l -b -q "CompareYields.C(13, $i)" #cosbachbar
#root -l -b -q "CompareYields.C(14, $i)" #dcabachbar
root -l -b -q "CompareYields.C(15, $i)" #dcav0topv
done
