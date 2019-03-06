# logrank-macro
This is a SAS macro that performs a clustered logrank test.  The manuscript that describes this software can be found here: https://www.ncbi.nlm.nih.gov/pubmed/21496938 

simsurvival.sas -- simulates example data to demonstrate the macro.

logrank_macro -- is the SAS macro to perform the clustered logrank test

By clustered data, I mean the randomization unit is the cluster. Example clusters include: patients clustered within physicians, students clustered within classrooms. This macro applies to designs where the physician/ classroom is randomized to the treatment of interest and outcomes are measured at the patient/ student level.  Please note that this is a different design than randomization within cluster.
