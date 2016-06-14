### INTRODUCTION

AlFinder reads the .mzML MS data, and extracts the MS/MS information of SLGGDSIMGIQLVSR peptide and its modified analogue according to the tag ions.

The tag ions are the y ions of precursors peptide's HCD MS/MS, 2 strategies for tag ions selection are used respectively:
1. The hit ions contain y9, y8, y7, y5, and y4, the precursor will be considered as a target
2. Count the hits of y1-y9, if the number of hits >=6, the precursor will be considered as a target


### AUTHOR
B. Pang (pangbo@sioc.ac.cn)
Shanghai Institute of Organic Chemistry, Chinese Academy of Sciences

### PYTHON AND LIBRARIES DEPENDENCE

Python 3.5
pymzML
