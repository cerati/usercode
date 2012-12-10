In order to use this tool you should log into a lxplus machine (or anything with access to LSF) 
and be in an AFS area (or anything readable from the batch machines). 

You should have one directory with the default cards (e.g. cards_def) and one with toys (e.g. cards_inj_stat)
which should contain the root files with toys (e.g. hwwsf_0j_shape_7TeV_PseudoData_sb.root).
Toys can be produced with LandS command: 
lands.exe -d card.txt -M Hybrid  -m 125 --minuitSTRATEGY 0 --bWriteToys 1 -n Name --nToysForCLsb 1000 --nToysForCLb 1 --singlePoint 1 --seed 12344  -rMin 0 -rMax 5

Edit directory names in getPseudoData.C. This is where the cards for pseudo eperiments are built.

Two scripts to submit jobs are available:
- dleGetSignificanceCombine.sh: test of limits, significance and signal strength
- gcRunNormTest.sh: test of post-fit normalizations

Results are then analyzed with other scripts:
- mergeLogs.sh, averageResultsNew.C: analyze the output of dleGetSignificanceCombine.sh
- plotNorm.C: analyze the output of gcRunNormTest.sh
