# pyFRESCO
This is a python implementation of the DistributeFoldx/DistributeRosetta scripts written by Hein Wijma.
This implementation includes some quality of life improvments that reduce the amounts of inputs needed, 
standardizes the input format between the two scrips, and includes a few additional warnings and error messages. 
This program is intended to distribute the calculation of a large number of mutations by FoldX and Rosetta ddg_monomer.
For the first phase of preparing the files for calculations, use a command like:

```
python DistributeFoldX.py Phase1 MyPrecious.pdb MyPreciousSelection.tab 100 /home/wijma/Fold_X/FoldX.linux64
```
or 
```
python DistributeRosettaddg.py Phase1 MyPrecious.pdb MyPreciousSelection.tab 100 FLAGrow3 /home/wijma/mini20101104/mini/bin/fix_bb_monomer_ddg.linuxgccrelease
```

This should PREPARE mutations of MyPrecious.pdb for residues specified in the File MyPreciousSelection.tab and 
distribute them over directories with each a 100 different mutations. For the second phase of collecting the results into a list, use a command like:

```
python DistributeFoldX.py Phase2 MyPrecious.pdb -5
```
or 
```
python DistributeRosetta.py Phase2 MyPrecious.pdb -5
```

This should COLLECT the mutations of MyPrecious.pdb assuming a cutoff of -5 KJ mol-1,resulting in FOUR lists:
- a complete list of all mutations
- a list of the mutations that are less than -5 kJ mol-1
- a list with the best mutation per position
- a list with the best mutation per position that are less than -5 kJ mol-1

The script makes use of `numpy` python package. the code was written and tested with version 1.23.3 but earlier or later versions will work. 
If needed, the script accepts the original input format as specified in the step-by-step protocol published at doi.org/10.1007/978-1-4939-7366-8_5 by setting the variable `LEGACY` to `True` inside the script.
