#mult=( "0.0-1.0" "1.0-5.0" "5.0-10.0" "10.0-15.0" "15.0-20.0" "20.0-25.0" "25.0-30.0" "30.0-35.0" "35.0-40.0" "40.0-50.0" "50.0-70.0" "70.0-100.0" "0-100" )
mult=( "0.0-1.0" "1.0-5.0" "5.0-10.0" "10.0-20.0" "20.0-30.0" "30.0-40.0" "40.0-50.0" "50.0-70.0" "70.0-100.0" "0-100" )
for i in `seq 0 11`
#for i in 0
do
fileYield="Yields/YieldEffCorrLHC22o_pass4_MinBias_Train108123_Omega_BkgParab_Mult"${mult[i]}"_Train110357_INELgt0.root"
echo ${fileYield}
scp ${fileYield} cdemart@gr3srv.ts.infn.it://data/dataalice/cdemart/Run3Analyses/OmegavsMult/Yields/.
done