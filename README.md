# Compressed-sensing-project

Files related to the project I did for the compressed sensing class during my masters.

These are separated into two categories: the report and the codes. 

## Report
Under 'Project - final version' is the report. Located here are the pdf files and the tex document. As is, the tex document cannot be compiled; missing are the figures and bib file. One may use the codes to generate figures that are essentially the same (as data may be randomly generated) as the ones in the report.

The title of the project is 'Survey of Results in Compressed Sensing using Total Variation Minimization'. In the report, a number of recent theoretical results are discussed, and supporting numerical experiments are conducted. Three key ideas covered are the use of the TV-norm (justification, results), effectiveness of common sampling strategies (uniform, etc.), and implementation using the split Bregman algorithm. 

## Codes
Under 'Project codes' are a number of Matlab scripts/functions. These are to test the stability and robustness of the signal recovery technique with examples with 1D signals and 2D images. For starters, the codes to run would be stable_n_robust_1D.m and stable_n_robust_2D.m. The initial comment blurb covers one way to call these two pieces of code.
