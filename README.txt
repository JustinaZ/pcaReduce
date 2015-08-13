
Depends: R (>= 3.0.1), pcaMethods (>= 1.50.0), mnormt (>= 1.5-1), mclust (>= 4.3)

==================INSTALL (on Mac’s)================================

to install, open terminal and type:
R

then type:
install.packages(<path_to_package>, repos = NULL, type="source")

p.s. this <pathtopackage> must be in “”

load the package by typing

library("pcaReduce")

p.p.s. type ?PCAreduce to see description

=================EXAMPLE===========================================

DATA:
is in Pollen2014.txt, cluster labels and additional information is in SupplementaryLabels.txt

CODE:

pcaReduce_example.R — an example code how to use pcaReduce package

===================================================================
J. 
