# Systematics

## Multi trial

The macro ```generatejson_systematics.cpp``` generates 500 different jsons with 500 different cut settings when given as input the limits for each cut

1. ``` >> cd DATA/ ```
2. ``` >> . dosystematics.sh```
3. ``` >> . doyieldscutvar.sh```
4. ``` >> cd MC/ ```
5. ``` >> . dosystematics.sh```
6. ``` >> . doyieldscutvar.sh```
7. in this folder:
   ``` >> . docorryieldscutvar.sh ```
8. ``` >> root MultiTrial.cxx ```

### Total spectra syst
``` root CalculateTotalSystematics.cxx ```

### pT uncorrelared
```root CalculateYields.cxx```

