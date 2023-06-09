#!/bin/bash 

o2-analysis-timestamp -b --aod-memory-rate-limit 300000000 --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-event-selection -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-multiplicity-table -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-centrality-table -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-track-propagation -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-pid-tof -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-pid-tof-base -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-pid-tpc-base -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-pid-tpc -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-lf-lambdakzerobuilder -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-lf-cascadebuilder -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
| o2-analysis-lf-cascqaanalysis -b  --configuration json://${PWD}/cascqaanalysis_Run3Data.json --time-limit 100 --resources-monitoring 2 --resources-monitoring-dump-interval 600 --aod-writer-json table_cascqaanalysis.json --fairmq-ipc-prefix .
#| o2-analysis-pid-tpc-full -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
#| o2-analysis-collision-converter -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \
#| o2-analysis-qa-event-track -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json #
#| o2-analysis-trackselection -b --configuration json://${PWD}/cascqaanalysis_Run3Data.json \








