# A scorecard to compare wet-lab experimental conditions 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13808354.svg)](https://doi.org/10.5281/zenodo.13808354)

This collection of functions creates a scorecard to compare wet-lab experimental conditions, useful to identify extreme or moderately changing expression values based on fold change and statistical significance. The goal is to narrow down the information to the relevant or extremely changing OMICS expression values by comparing several pairwise experimental conditions. On this reduced set of entries, the library could produce additional data-driven insights for investigating the experimental findings.

![Scorecard example](examples/Scorecard_examples.png?raw=True)

## Quick introduction
At the beginning of the file include the import statement, for example:
```
import scorecard_functions
```
Initially, create the parameter set to tune the scorecard features. 
The following line of code could be used to generate the parameter dictionary with default values: 
```
param_dict=scorecard_functions.generate_parameters()
```
Afterward, one can modify the default parameters with users preferences. For example, one should define the information regarding name and location of the CSV file and dataset which will be imported as a Pandas DataFrame:
*    'base_dir' # Path of the directory where the CSV resides
*   'filename'# filename of the CSV with the log2 Fold change and Adjusted p-values. Should include gene names or gene ID as separate column
*    'zscore' # Perform standardization on log2 F.C. values
*    'Treatment1 name' # custom name of Treatment 1
*    'Treatment2 name' # custom name of Treatment 2
*    'Control name' # custom name of Control/Baseline
*    'fig_size' # figure size and limits
*    'incl aver' # Include intermediate regions (between F.C. thresholds)  
*    'save_dir' # path of the directory where saving scorecards and JSON outcomes
*    'colors' # Color codes of the markers in each area of the scorecard (inside the regions of interest)
*    'other_colors' # Color of the markers outside the regions of interest
*    'th_fold_change' # log2 fold change standard threshold (another threshold will be automatically added based on 'multiplication factor')
*    'th_significance'  # Significance threshold i.e. 0.05 or lower
*    'font_size_quadrants'  # font size of the genes over the scorecard
*    'font_size_examples'  # Scorecard legend font size when typing the areas of interest
*    'marker_trasp'  # transparency of markers
*    'rect_trasp'  # transparency of shaded areas of the scorecard
*    'rect_colors'  # colors of shaded areas of the scorecard
*    'markers'  #  markers of the scorecard
*    'markers_sizes'  # Size of the markers plotted on the scorecard
*    'is_example'  # in scorecard legend use color names or charaters as codes
*    'use_notation'  # Use custom names inserted in 'Treatment1 name','Treatment2 name','Control name' when plotting axes labels
*    'gene_name'  # Name of the column of the DataFrame containing genes or proteins names or IDs (one per row)
*    'multiplication factor' # Factor to multiply the log2 Fold change threshold to detect extreme values (setting it to 1 will produce a four-way plot, not the scorecard)
*    'CSV delimiter'  # delimiter in the CSV file (usually ',')
*    'Log Epsilon' # tiny value to adjust the log calculation (during Volcano plots)
*    'Scorecard title' # add some short description regarding the experimental conditions on the scorecard title

Other important parameters are the column names in the CSV file containing the log2 fold change and the adjusted p-values:
*    'FC cond x' # replace with exact column name containing log2 fold change Treatment1 vs Control
*    'FC cond y' # replace with exact column name containing log2 fold change Treatment2 vs Control
*    'padj cond x' # replace with exact column name containing adj p-values Treatment1 vs Control
*    'padj cond y' # replace with exact column name containing adj p-values Treatment2 vs Control

Indeed, your dataset should contain pre-computed fold change values (and related adjusted p-values) as columns in the CSV file.

For example the notation could be "Treatment1vsCtrl_Log2FoldChange", "Treatment2vsCtrl_Log2FoldChange", "Treatment4vsCtrl_Log2FoldChange", etc... (with related p-values in a similar record). 
Remember that the CSV file should have the first row as header, reporting the column names matching the information provided in the dictionary of parameters.
| Gene ID | "Treatment1vsCtrl_Log2FoldChange" | "Treatment2vsCtrl_Log2FoldChange" | "Treatment3vsCtrl_Log2FoldChange" |"Treatment1vsCtrl_padj" |"Treatment2vsCtrl_padj" |"Treatment3vsCtrl_padj" |
|--------:|-----------|-----------|-----------|-----------|-----------|-----------|
|     GN1| 1.3 | 2.5 | 1.6 | 0.02 | 0.008 | 0.1 |
|     GN2| 5.9    | 9.0 | 2.5 | 0.56 | 0.036 | 0.07 |
|     GN3| -0.2       | -6.1 | -5.7 | 0.001 | 0.03 | 0.09 |
| ... | ... | ... |... |... |... |... |

Now load the CSV data into the computer memory as Pandas Dataframe by typing:
```
df=scorecard_functions.data_loading(param_dict)
```

As a side note, the Scorecard library contains a utility function to generate a list of combinations (without repetitions) for iterating over all samples.
It is suggested to create a legend of the color codes and the regions of interest by typing:
```
scorecard_functions.scorecard_legend(param_dict)
```

Initially, run scorecard single quandrant calculations using the following code:
```
scorecard_functions.scorecard(df,param_dict)
```
This function will store statistically significant entries for each quadrant together with plots of individual quandrants.
Individual plots could be useful in case of crowded scorecards.

After creating all four quadrants with previous function, one should execute scorecard reconstruction from individual quandrants using this function call:
```
scorecard_functions.reconstruct_scorecard(df,param_dict)
```
The parameters "use_figsize" and "figsize_factor" allow for personalization of the image dimension.
Default dimension is set by "fig_size" in the initial parameters, otherwise set "use_figsize" as False per customization.

Repeat the procedure for scorecard reconstruction on all comparisons requested. Theoretically one should have a main folder and several subfolders each one comparing two experimental conditions.
The folders structure should be as the following:

            main_folder
             |
             |
             |-------Exp. Comparison 1
             |
             |-------Exp. Comparison 2
             |
             |-------Exp. Comparison 3

One could create a general overview of all wet-lab experimental conditions by calling the function (only for the scorecard, not for the four-way plot):
```
scorecard_functions.multiple_view(main_folder)
```
Remember that "main_folder" is a string containing the full path of the folder hosting all sub-directories with experimental comparisions.

To build Volcano plots and highlight genes or entries belonging to the regions of interest of the scorecard type:
```
scorecard_functions.make_volcano(main_folder)
```

Another function provides a single gene or entry analysis as barplot
```
scorecard_functions.multiple_bars(main_folder)
```
To get a summary of the number of entries/genes indentified by the scorecard type:
```
scorecard_functions.count_frequencies(main_folder)
```
Also, the summary of the number of entries/genes indentified by the scorecard can be visualized as heatmap:
```
scorecard_functions.quadrants_heatmap(main_folder)
```

To create a final report on single genes/entries run the following function. Remember to call it after *count_frequencies* because
it re-uses information previously created by this function. Output file could be in CSV or Excel format showing total occurences.
```
scorecard_functions.common_entries(main_folder)
```
This function could also include a separate report about uncommon entries of the scorecard.

In case of experimental comparisons are over time or one wants to sum-up more experiments, this function can track modifications of the expression values over the scorecard.
The input folder should contain sub-folders with pre-computed data regarding the comparisons (common_entries should be called first).
```
scorecard_functions.track_over_exper(main_folder)
```
