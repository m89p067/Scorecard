import pandas as pd
import scorecard_functions

scorecard_functions.help()

param_dict=scorecard_functions.generate_parameters()
print(param_dict) # standard parameters
test_col_folders=False
unq = scorecard_functions.identified_comparisons(['VAN','DEL','CPT','LZD'])

for indice, valore in enumerate(unq):
    e_x,e_y=valore[0],valore[1]
    print('----------------> WORKING ON ',indice+1,' OVER ',len(unq))
    param_dict['FC cond x']="48hr "+e_x+" vs 48 hr NO ABX log2FC" # x-axis
    param_dict['FC cond y']="48hr "+e_y+" vs 48 hr NO ABX log2FC" #y-axis
    param_dict['padj cond x']="48hr "+e_x+" vs 48 hr NO ABX Adj. p-value"
    param_dict['padj cond y']="48hr "+e_y+" vs 48 hr NO ABX Adj. p-value"
    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset1"
    param_dict['filename']="Cartel1.CSV"
    param_dict['Treatment1 name']=e_x # treatment1: Vancomycin 400 mcg/ml
    param_dict['Treatment2 name']=e_y # treatment2: Linezolid 20 mcg/ml
    param_dict['Control name']="NO ABX" # No antibiotic (baseline condition)
    param_dict['gene_name']="Orf"
    param_dict['CSV delimiter']=";"
    param_dict['Scorecard title']="After 48hr"
    param_dict['multiplication factor']=1.5 # STANDARD SCORECARD
    param_dict['th_significance']=0.05
    if test_col_folders:
        param_dict['is_example']=True
        mydir="D:/GED_VIS_PROJECT/Dataset1/Scorecard_STANDARD2"
        param_dict['save_dir']=mydir
        param_dict['colors']=['green','lime','blue','red','magenta']
    else:
        mydir="D:/GED_VIS_PROJECT/Dataset1/Scorecard_STANDARD"
        param_dict['save_dir']=mydir
    print(param_dict) # Updated parameters


    df=scorecard_functions.data_loading(param_dict)

    scorecard_functions.scorecard_legend(param_dict)
    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
mydir="D:/GED_VIS_PROJECT/Dataset1/Scorecard_STANDARD"
scorecard_functions.reconstruct_scorecard(mydir,use_figsize=False,figsize_factor=2)
scorecard_functions.multiple_view(mydir,fs_size=14,single_quadr=True)
scorecard_functions.ranking_bars(mydir,title_size=26,remove_string='SAUSA300_')
scorecard_functions.make_volcano(mydir)
scorecard_functions.multiple_bars(mydir, try_adj_test=False,text_adj_x=0.0,text_adj_y=0.0,remove_string='SAUSA300_')
scorecard_functions.largest_diff(mydir,top_entries=0,do_excel=True)
scorecard_functions.count_frequencies(mydir)
scorecard_functions.common_entries(mydir,do_excel=True)
scorecard_functions.quadrants_heatmap(mydir)
scorecard_functions.track_over_exper(mydir,is_time=False)



scorecard_functions.help()
param_dict=scorecard_functions.generate_parameters()
print(param_dict) # standard parameters


for indice, valore in enumerate(unq):
    e_x,e_y=valore[0],valore[1]
    print('----------------> WORKING ON ',indice+1,' OVER ',len(unq))
    param_dict['FC cond x']="48hr "+e_x+" vs 48 hr NO ABX log2FC" # x-axis
    param_dict['FC cond y']="48hr "+e_y+" vs 48 hr NO ABX log2FC" #y-axis
    param_dict['padj cond x']="48hr "+e_x+" vs 48 hr NO ABX Adj. p-value"
    param_dict['padj cond y']="48hr "+e_y+" vs 48 hr NO ABX Adj. p-value"
    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset1"
    param_dict['filename']="Cartel1.CSV"
    param_dict['Treatment1 name']=e_x # treatment1: Vancomycin 400 mcg/ml
    param_dict['Treatment2 name']=e_y # treatment2: Linezolid 20 mcg/ml
    param_dict['Control name']="NO ABX" # No antibiotic (baseline condition)
    param_dict['gene_name']="Orf"
    param_dict['CSV delimiter']=";"
    param_dict['Scorecard title']="After 48hr"
    param_dict['multiplication factor']=1.5 
    param_dict['th_significance']=0.05
    param_dict['incl aver']=True  # <-------------------------- FULL SCORECARD
    if test_col_folders:
        param_dict['is_example']=True
        mydir="D:/GED_VIS_PROJECT/Dataset1/Scorecard_FULL2"
        param_dict['save_dir']=mydir
        param_dict['colors']=['green','lime','blue','red','magenta']
    else:
        mydir="D:/GED_VIS_PROJECT/Dataset1/Scorecard_FULL"
        param_dict['save_dir']=mydir
    print(param_dict) # Updated parameters


    df=scorecard_functions.data_loading(param_dict)

    scorecard_functions.scorecard_legend(param_dict)
    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
mydir="D:/GED_VIS_PROJECT/Dataset1/Scorecard_FULL"
scorecard_functions.reconstruct_scorecard(mydir,use_figsize=False,figsize_factor=2)
scorecard_functions.count_frequencies(mydir)
scorecard_functions.common_entries(mydir,do_excel=True)
scorecard_functions.multiple_view(mydir,single_quadr=True,confirm_ave=True)
scorecard_functions.ranking_bars(mydir,title_size=26,remove_string='SAUSA300_')
scorecard_functions.make_volcano(mydir)
scorecard_functions.multiple_bars(mydir,remove_string='SAUSA300_')
scorecard_functions.largest_diff(mydir,top_entries=0,do_excel=True)
scorecard_functions.quadrants_heatmap(mydir)
scorecard_functions.track_over_exper(mydir,is_time=False)

param_dict=scorecard_functions.generate_parameters()
print(param_dict) # standard parameters




