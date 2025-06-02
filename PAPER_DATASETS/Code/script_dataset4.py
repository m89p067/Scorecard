import pandas as pd
import scorecard_functions

##mydir="D:/GED_VIS_PROJECT/Dataset7/Scorecard_STANDARD"
##
##param_dict=scorecard_functions.generate_parameters()
##for e_x,e_y in zip(['Ag_t03','Ag_t03','Ag_t90','Ag_t90','V2A_t03','V2A_t90'],
##                   ['AgXX_t03','V2A_t03','AgXX_t90','V2A_t90','AgXX_t03','AgXX_t90']):
##    param_dict['FC cond x']=e_x+"_logFC" # x-axis
##    param_dict['FC cond y']=e_y+"_logFC" #y-axis    
##    param_dict['padj cond x']=e_x+"_padj"
##    param_dict['padj cond y']=e_y+"_padj"
##    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset7"
##    param_dict['filename']="out.csv"
##    param_dict['Treatment1 name']=e_x # treatment1
##    param_dict['Treatment2 name']=e_y # treatment2
##    param_dict['Control name']="Ctrl" # Baseline condition
##    param_dict['gene_name']="GN"
##    param_dict['CSV delimiter']=","
##    param_dict['Scorecard title']="Nosocomial Ent. faec. under metal stress"
##    param_dict['multiplication factor']=2.5
##    param_dict['th_fold_change']=2.5
##    param_dict['th_significance']=0.001
##    param_dict['use_notation']=True
##    param_dict['is_example']=False
##    param_dict['save_dir']=mydir
##    param_dict['colors']=['turquoise','crimson','orangered','navy','greenyellow'] # Color codes of the markers in each area of the scorecard
##    param_dict['other_colors']=['lavenderblush','mistyrose','plum'] # colors of points outside the regions of interest
##    print(param_dict) # Updated parameters
##
##
##    df=scorecard_functions.data_loading(param_dict)
##
##    scorecard_functions.scorecard_legend(param_dict)
##    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
##    scorecard_functions.reconstruct_scorecard(mydir)
##
##scorecard_functions.multiple_view(mydir,single_quadr=True)
##scorecard_functions.make_volcano(mydir)
##
##
##scorecard_functions.multiple_bars(mydir,remove_string='EF.peg.')
##scorecard_functions.ranking_bars(mydir,title_size=26,remove_string='EF.peg.')
##scorecard_functions.count_frequencies(mydir)
##scorecard_functions.quadrants_heatmap(mydir)
##scorecard_functions.common_entries(mydir,do_excel=True,fs_size=4)
##scorecard_functions.track_over_exper(mydir,is_time=False)
##scorecard_functions.largest_diff(mydir,top_entries=0,do_excel=True)
##
##mydir="D:/GED_VIS_PROJECT/Dataset7/Scorecard_FULL"
##
##param_dict=scorecard_functions.generate_parameters()
##for e_x,e_y in zip(['Ag_t03','Ag_t03','Ag_t90','Ag_t90','V2A_t03','V2A_t90'],
##                   ['AgXX_t03','V2A_t03','AgXX_t90','V2A_t90','AgXX_t03','AgXX_t90']):
##    param_dict['FC cond x']=e_x+"_logFC" # x-axis
##    param_dict['FC cond y']=e_y+"_logFC" #y-axis    
##    param_dict['padj cond x']=e_x+"_padj"
##    param_dict['padj cond y']=e_y+"_padj"
##    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset7"
##    param_dict['filename']="out.csv"
##    param_dict['Treatment1 name']=e_x # treatment1
##    param_dict['Treatment2 name']=e_y # treatment2
##    param_dict['Control name']="Ctrl" # Baseline condition
##    param_dict['gene_name']="GN"
##    param_dict['CSV delimiter']=","
##    param_dict['Scorecard title']="Nosocomial Ent. faec. under metal stress"
##    param_dict['multiplication factor']=2
##    param_dict['th_fold_change']=2
##    param_dict['th_significance']=0.001
##    param_dict['use_notation']=True
##    param_dict['is_example']=False
##    param_dict['save_dir']=mydir
##    param_dict['incl aver']=True
##    param_dict['colors']=['turquoise','crimson','orangered','navy','greenyellow'] # Color codes of the markers in each area of the scorecard
##    param_dict['other_colors']=['lavenderblush','mistyrose','plum'] # colors of points outside the regions of interest
##    print(param_dict) # Updated parameters
##
##
##    df=scorecard_functions.data_loading(param_dict)
##
##    scorecard_functions.scorecard_legend(param_dict)
##    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
##    scorecard_functions.reconstruct_scorecard(mydir)
##
##scorecard_functions.multiple_view(mydir,single_quadr=True)
##scorecard_functions.make_volcano(mydir)
##scorecard_functions.multiple_bars(mydir,remove_string='EF.peg.')
##scorecard_functions.ranking_bars(mydir,title_size=26,remove_string='EF.peg.')
##scorecard_functions.count_frequencies(mydir)
##scorecard_functions.quadrants_heatmap(mydir)
##scorecard_functions.common_entries(mydir,do_excel=True,fs_size=4)
##scorecard_functions.track_over_exper(mydir,is_time=False)
##scorecard_functions.largest_diff(mydir,top_entries=0,do_excel=True)



mydir="D:/GED_VIS_PROJECT/Dataset7/Scorecard_FULL_time_comp"

param_dict=scorecard_functions.generate_parameters()
for e_x,e_y in zip(['V2A_t90','V2A_t90','V2A_t90','V2A_t90','V2A_t90'],
                   ['V2A_t06','V2A_t12','V2A_t24','V2A_t60','V2A_t03']):
    param_dict['FC cond x']=e_x+"_logFC" # x-axis
    param_dict['FC cond y']=e_y+"_logFC" #y-axis    
    param_dict['padj cond x']=e_x+"_padj"
    param_dict['padj cond y']=e_y+"_padj"
    param_dict['base_dir']="D:/GED_VIS_PROJECT/Dataset7"
    param_dict['filename']="out.csv"
    param_dict['Treatment1 name']=e_x # treatment1
    param_dict['Treatment2 name']=e_y # treatment2
    param_dict['Control name']="Ctrl" # Baseline condition
    param_dict['gene_name']="GN"
    param_dict['CSV delimiter']=","
    param_dict['Scorecard title']="Nosocomial Ent. faec. under metal stress"
    param_dict['multiplication factor']=2
    param_dict['th_fold_change']=2
    param_dict['th_significance']=0.001
    param_dict['use_notation']=True
    param_dict['is_example']=False
    param_dict['incl aver']=True
    param_dict['save_dir']=mydir
    param_dict['colors']=['indigo','darkorange','chocolate','lightcoral','forestgreen'] # Color codes of the markers in each area of the scorecard
    param_dict['other_colors']=['mintcream','pink','bisque'] # colors of points outside the regions of interest
    print(param_dict) # Updated parameters


    df=scorecard_functions.data_loading(param_dict)

    scorecard_functions.scorecard_legend(param_dict)
    scorecard_functions.scorecard(df,param_dict) # Intermediate data, use individual quandrants in case of crowded plots
scorecard_functions.reconstruct_scorecard(mydir)

scorecard_functions.multiple_view(mydir,single_quadr=True)
scorecard_functions.make_volcano(mydir)


scorecard_functions.multiple_bars(mydir,remove_string='EF.peg.')
scorecard_functions.ranking_bars(mydir,title_size=26,remove_string='EF.peg.')
scorecard_functions.count_frequencies(mydir)
scorecard_functions.quadrants_heatmap(mydir)
scorecard_functions.common_entries(mydir,do_excel=True,fs_size=4)
scorecard_functions.track_over_exper(mydir)
scorecard_functions.largest_diff(mydir,top_entries=0,do_excel=True)
