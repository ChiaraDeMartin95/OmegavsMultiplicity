# Systematics

## Multi trial

### Data:
The macro ```generatejson_systematics.cpp``` generates 500 different jsons with 500 different cut settings when given as input the limits for each cut

## DATA:
0. ``` >> cd DATA/ ```
1. ``` >> . dosystematics.sh```
2. ``` >> . doyieldscutvar.sh```

## MC:
3. ``` >> cd MC/ ```
4. ``` >> . dosystematics.sh```
5. ``` >> . doyieldscutvar.sh```


6. ``` >> . docorryieldscutvar.sh ```

7. ``` >> root MultiTrial.cxx ```

## Total spectra syst
``` root CalculateTotalSystematics.cxx ```

## pT uncorrelared

1. ```root CalculateYields.cxx```

