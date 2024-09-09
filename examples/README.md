# Folder contents (toy example)

The folder contains a CSV file that could be used to generate an analysis sequence for exploring Scorecard library functions.

The dataset was created artificially, and it is not constituted by computations on biological samples: each column contains fake Fold Change and p-values representing ten bacteria cultured in the laboratory and exposed to a certain experimental condition (versus a control or baseline). The ten bacteria cultures could be compared to each other. Gene IDs are alphanumeric codes built by randomly mixing numbers and letters.

The column names could be interpreted as post-processed RNA-seq data from bacteria cultured on Petri dishes:

| Abbreviations | Bacteria |
|-----:|---------------|
|MRSA| Methicillin-Resistant Staphylococcus aureus|
|VRE| Vancomycin-Resistant Enterococcus|
|EHEC| Enterohemorrhagic Escherichia coli|
|VRSA| Vancomycin-Resistant Staphylococcus aureus|
|CRE| Carbapenem-Resistant Enterobacteriaceae|
|MDR_TB | Multidrug-Resistant Mycobacterium tuberculosis|
|XDR_TB | Extensively Drug-Resistant Mycobacterium tuberculosis|
|GAS| Group A Streptococcus|
|GBS| Group B Streptococcus|
|NTHi| Nontypeable Haemophilus influenzae|

As a side note, the Scorecard library contains a utility function to generate a list of combinations (without repetitions) for iterating over all samples.
