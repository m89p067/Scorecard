import pandas as pd
import scorecard_functions

exp_list=[
"OxvsUn",
"LorvsUn",
"Cmpd8vsUn",
"Cmpd8vsUn",
"Cmpd8_OxvsUn",
"Lor_OxvsUn",
"OxvsCmpd8",
"Cmpd8_OxvsOx",
"Cmpd8_OxvsOx",
"Lor_OxvsOx",
"OxvsLor",
"Cmpd8vsLor",
"Cmpd8_OxvsLor",
"Lor_OxvsLor",
"Cmpd8_OxvsCmpd8",
"Lor_OxvsCmpd8",
"Lor_OxvsCmpd8_Ox"
]

unq = scorecard_functions.identified_comparisons(exp_list)
test_col_folders=False

param_dict=scorecard_functions.generate_parameters()
#print(param_dict) # standard parameters

for indice, valore in enumerate(unq):
    print('----------------> WORKING ON ',indice+1,' OVER ',len(unq))
    e_x,e_y=valore[0],valore[1]
    param_dict['FC cond x']=e_x+"_log2FoldChange" # x-axis
    param_dict['FC cond y']=e_y+"_log2FoldChange" #y-axis
    param_dict['padj cond x']=e_x+"_padj"
    param_dict['padj cond y']=e_y+"_padj"
    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset3"
    param_dict['filename']="USA100_comp.CSV"
    param_dict['Treatment1 name']=e_x # treatment1
    param_dict['Treatment2 name']=e_y # treatment2
    param_dict['Control name']="Un" # Baseline condition
    param_dict['gene_name']="gene_id"
    param_dict['CSV delimiter']=";"
    param_dict['Scorecard title']="MRSA USA100 cultures"
    param_dict['multiplication factor']=2 # Increasing the THRESHOLD of F.C.
    param_dict['th_significance']=0.05
    param_dict['use_notation']=False
    param_dict['incl aver']=True # <----------------------------
    if test_col_folders:
        mydir="D:/GED_VIS_PROJECT/Dataset3/Scorecard_FULL2"
        param_dict['is_example']=True
        param_dict['save_dir']=mydir
        param_dict['colors']=['green','lime','blue','red','magenta']
    else:
        mydir="D:/GED_VIS_PROJECT/Dataset3/Scorecard_FULL"
        param_dict['save_dir']=mydir
        param_dict['colors']=['goldenrod','limegreen','dodgerblue','sandybrown','crimson']
        param_dict['other_colors']=['darkcyan','indigo','orangered']
    print(param_dict) # Updated parameters


    df=scorecard_functions.data_loading(param_dict)

    scorecard_functions.scorecard_legend(param_dict)
    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
scorecard_functions.reconstruct_scorecard(mydir)
scorecard_functions.multiple_view(mydir,fs_size=6)
scorecard_functions.make_volcano(mydir)
scorecard_functions.multiple_bars(mydir)
scorecard_functions.count_frequencies(mydir)
scorecard_functions.quadrants_heatmap(mydir,above_lab_rot=-90,horiz_alig="center",font_counts=2.75)
scorecard_functions.common_entries(mydir,do_excel=True,linewidth=0.6,fs_size=3)
scorecard_functions.track_over_exper(mydir,is_time=False)

param_dict=scorecard_functions.generate_parameters()
#print(param_dict) # standard parameters

for indice, valore in enumerate(unq):
    print('----------------> WORKING ON ',indice+1,' OVER ',len(unq))
    e_x,e_y=valore[0],valore[1]
    param_dict['FC cond x']=e_x+"_log2FoldChange" # x-axis
    param_dict['FC cond y']=e_y+"_log2FoldChange" #y-axis
    param_dict['padj cond x']=e_x+"_padj"
    param_dict['padj cond y']=e_y+"_padj"
    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset3"
    param_dict['filename']="USA100_comp.CSV"
    param_dict['Treatment1 name']=e_x # treatment1
    param_dict['Treatment2 name']=e_y # treatment2
    param_dict['Control name']="Un" # Baseline condition
    param_dict['gene_name']="gene_id"
    param_dict['CSV delimiter']=";"
    param_dict['Scorecard title']="MRSA USA100 cultures"
    param_dict['multiplication factor']=1.5 # Fixing the THRESHOLD of F.C. as the first example
    param_dict['th_significance']=0.05
    param_dict['use_notation']=False
    param_dict['incl aver']=True # <----------------------------
    if test_col_folders:
        mydir="D:/GED_VIS_PROJECT/Dataset3/Scorecard_FULL3"
        param_dict['is_example']=True
        param_dict['save_dir']=mydir
        param_dict['colors']=['green','lime','blue','red','magenta']
    else:
        mydir="D:/GED_VIS_PROJECT/Dataset3/Scorecard_FULL4"
        param_dict['save_dir']=mydir
        param_dict['colors']=['goldenrod','limegreen','dodgerblue','sandybrown','crimson']
        param_dict['other_colors']=['darkcyan','indigo','orangered']
    print(param_dict) # Updated parameters


    df=scorecard_functions.data_loading(param_dict)

    scorecard_functions.scorecard_legend(param_dict)
    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
scorecard_functions.reconstruct_scorecard(mydir)
scorecard_functions.multiple_view(mydir,fs_size=6)
scorecard_functions.make_volcano(mydir)
scorecard_functions.multiple_bars(mydir)
scorecard_functions.count_frequencies(mydir)
scorecard_functions.quadrants_heatmap(mydir,above_lab_rot=-90,horiz_alig="center",font_counts=2.75)
scorecard_functions.common_entries(mydir,do_excel=True,linewidth=0.6,fs_size=3)
scorecard_functions.track_over_exper(mydir,th_sel=6,is_time=False)
