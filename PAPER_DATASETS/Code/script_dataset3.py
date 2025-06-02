import pandas as pd
import scorecard_functions

param_dict=scorecard_functions.generate_parameters()
print(param_dict)
unq = scorecard_functions.identified_comparisons(['Lac','Glu','Suc','Xyl'])


mydir2='D:/GED_VIS_PROJECT/Dataset5/Scorecard_STANDARD'
param_dict=scorecard_functions.generate_parameters()
for indice, valore in enumerate(unq):
    e_x,e_y=valore[0],valore[1]
    print('----------------> WORKING ON ',indice+1,' OVER ',len(unq))
    param_dict['FC cond x']="TSB_versus_"+e_x+"_log2FoldChange" # x-axis
    param_dict['FC cond y']="TSB_versus_"+e_y+"_log2FoldChange" #y-axis    
    param_dict['padj cond x']="TSB_versus_"+e_x+"_padj"
    param_dict['padj cond y']="TSB_versus_"+e_y+"_padj"
    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset5"
    param_dict['filename']="preproc_data_red.CSV"
    param_dict['Treatment1 name']=e_x # treatment1
    param_dict['Treatment2 name']=e_y # treatment2
    param_dict['Control name']="TSB" # Baseline condition
    param_dict['gene_name']="Name"
    param_dict['CSV delimiter']=";"
    param_dict['Scorecard title']="Gene expr. of S. mutans"
    param_dict['multiplication factor']=2
    param_dict['th_fold_change']=2
    param_dict['th_significance']=0.001
    param_dict['use_notation']=True
    param_dict['is_example']=False
    param_dict['save_dir']=mydir2
    param_dict['colors']=['tomato','darkorange','olivedrab','powderblue','hotpink']
    param_dict['other_colors']=['peru','khaki','lavender']
    print(param_dict) # Updated parameters


    df=scorecard_functions.data_loading(param_dict)

    scorecard_functions.scorecard_legend(param_dict)
    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
scorecard_functions.reconstruct_scorecard(mydir2)
scorecard_functions.multiple_view(mydir2,fs_size=14,single_quadr=True)
scorecard_functions.make_volcano(mydir2)
scorecard_functions.multiple_bars(mydir2)
scorecard_functions.count_frequencies(mydir2)
scorecard_functions.common_entries(mydir2,do_excel=True,fs_size=4)
scorecard_functions.quadrants_heatmap(mydir2)
scorecard_functions.track_over_exper(mydir2,is_time=False)
scorecard_functions.ranking_bars(mydir2,title_size=26)
scorecard_functions.largest_diff(mydir2,top_entries=0,do_excel=True)
mydir4='D:/GED_VIS_PROJECT/Dataset5/Scorecard_FULL'
param_dict=scorecard_functions.generate_parameters()
for indice, valore in enumerate(unq):
    e_x,e_y=valore[0],valore[1]
    print('----------------> WORKING ON ',indice+1,' OVER ',len(unq))
    param_dict['FC cond x']="TSB_versus_"+e_x+"_log2FoldChange" # x-axis
    param_dict['FC cond y']="TSB_versus_"+e_y+"_log2FoldChange" #y-axis    
    param_dict['padj cond x']="TSB_versus_"+e_x+"_padj"
    param_dict['padj cond y']="TSB_versus_"+e_y+"_padj"
    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset5"
    param_dict['filename']="preproc_data_red.CSV"
    param_dict['Treatment1 name']=e_x # treatment1
    param_dict['Treatment2 name']=e_y # treatment2
    param_dict['Control name']="TSB" # Baseline condition
    param_dict['gene_name']="Name"
    param_dict['CSV delimiter']=";"
    param_dict['Scorecard title']="Gene expr. of S. mutans"
    param_dict['multiplication factor']=2
    param_dict['th_fold_change']=2
    param_dict['th_significance']=0.001
    param_dict['use_notation']=True
    param_dict['is_example']=False
    param_dict['incl aver']=True #<------------------------------
    param_dict['save_dir']=mydir4
    param_dict['colors']=['tomato','darkorange','olivedrab','powderblue','hotpink']
    param_dict['other_colors']=['peru','khaki','lavender']
    print(param_dict) # Updated parameters


    df=scorecard_functions.data_loading(param_dict)

    scorecard_functions.scorecard_legend(param_dict)
    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
scorecard_functions.reconstruct_scorecard(mydir4)
scorecard_functions.multiple_view(mydir4,single_quadr=True,confirm_ave=True)
scorecard_functions.make_volcano(mydir4)
scorecard_functions.multiple_bars(mydir4)
scorecard_functions.count_frequencies(mydir4)
scorecard_functions.common_entries(mydir4,do_excel=True,fs_size=4)
scorecard_functions.quadrants_heatmap(mydir4)
scorecard_functions.track_over_exper(mydir4,is_time=False)
scorecard_functions.ranking_bars(mydir4,title_size=26)
scorecard_functions.largest_diff(mydir4,top_entries=0,do_excel=True)
