import numpy as np
import pandas as pd
from os.path import join,isdir
from scipy.stats import zscore
from os import listdir,getcwd,makedirs,scandir
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.markers import MarkerStyle
import json
import ctypes
from statsmodels.stats.multitest import multipletests,fdrcorrection
from scipy.stats import ttest_ind,levene,linregress
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
import pdb
def help():
    with open("help.txt", "r", encoding="utf-8") as file:
        for line in file:
            stringa=line.strip()
            if stringa.startswith(">>"):
                print('\n '+stringa+' \n')
            else:
                print(stringa)
    print("\n\n")
def generate_parameters():
    the_dict={}
    # Other information, assuming CSV file contains log2 fold change already computed
    the_dict['FC cond x']='log2 fold change T1 vs T0' # replace with exact column name containing log2 fold change Treatment1 vs Control
    the_dict['FC cond y']='log2 fold change T2 vs T0' # replace with exact column name containing log2 fold change Treatment1 vs Control
    the_dict['padj cond x']='padj T1 vs T0' # replace with exact column name containing adj p-values Treatment1 vs Control
    the_dict['padj cond y']='padj T2 vs T0' # replace with exact column name containing adj p-values Treatment1 vs Control
 
    the_dict['base_dir']=getcwd() # Path of the directory where the CSV resides
    the_dict['filename']='*.csv' # filename of the CSV with the log2 Fold change and Adjusted p-values. Should include gene names or gene ID as separate column
    the_dict['zscore']=False # Perform standardization on log2 F.C. values
    the_dict['Treatment1 name']='T1' # custom name of Treatment 1
    the_dict['Treatment2 name']='T2' # custom name of Treatment 2
    the_dict['Control name']='T0' # custom name of Control/Baseline
    the_dict['fig_size']=10 # figure size and limits
    the_dict['save_dir']=getcwd() # path of the directory where saving scorecards and JSON outcomes
    #the_dict['colors']=['seagreen','darkturquoise','dodgerblue','red','fuchsia'] # Color codes of the markers in each area of the scorecard
    the_dict['colors']=['green','teal','blue','red','magenta'] # Color codes of the markers in each area of the scorecard
    the_dict['other_colors']=['darkorange','gold','lightpink'] # colors of points outside the regions of interest
    the_dict['th_fold_change']=2 # log2 fold change standard threshold (it will be automatically added another threshold based on 'multiplication factor')
    the_dict['th_significance']=0.01 # Significance threshold i.e. 0.05 or lower
    the_dict['font_size_quadrants']=8 # font size of the genes over the scorecard
    the_dict['font_size_examples']=16 # Scorecard legend font size when typing the areas of interest
    the_dict['marker_trasp']=[0.95,0.8,0.65,0.45] # transparency of markers
    the_dict['rect_trasp']=[0.25,0.15,0.4] # transparency of shaded areas of the scorecard
    the_dict['rect_colors']=['silver','lightgrey','darkgrey'] # colors of shaded areas of the scorecard
    the_dict['markers']=['o','^','s','*','X','p','.'] #  markers of the scorecard
    the_dict['markers_sizes']= [14,12,10,8] # Size of the markers
    the_dict['is_example']=False # in scorecard legend use color names or charaters as codes for regions of interest (False means use letters A,B,C,D,E for the regions of inetrest)
    the_dict['use_notation']=True # Use cutom names inserted in 'Treatment1 name','Treatment2 name','Control name' when plotting axes labels
    the_dict['gene_name']='GN' # Name of the column name containing genes or proteins (one per row)
    the_dict['multiplication factor']=2 # Factor to multiply the log2 Fold change threshold for detecting extreme values
    the_dict['CSV delimiter']=';' # delimiter in the CSV file (usually ',')
    the_dict['Log Epsilon']=1e-10 # tiny value to adjust the log calculation
    the_dict['Scorecard title']='' # add some short description regarding the experimental conditions on the scorecard title    
    return the_dict
def check_point(a,b,xp,yp):
    # line of equation y=f(x)=a*x+b
    if yp - ( a * xp + b ) < 0: #BELOW LINE
        return 1 
    else:
        return 0
def log2_with_epsilon(arr,EPSILON=1e-10):
    result = np.log2(arr + EPSILON)
    return result
def log10_with_epsilon(arr,EPSILON=1e-10):
    result = np.log10(arr + EPSILON)
    return result
def flatten(xss):
    return [x for xs in xss for x in xs]
def data_loading(info_dict):
    '''
    After creating a dict of parameters, load the CSV file containing Fold Change and adjusted p-values
    The routine will check for NaN in the data, and output a DataFrame for further usage (building the scorecard).
    '''
    load_folder=info_dict['base_dir']
    if load_folder[-1]!="/":
        load_folder=load_folder+"/"
    full_path =join(load_folder,info_dict['filename'])
    my_df=pd.read_csv(full_path, header=0,index_col=False,delimiter=info_dict['CSV delimiter'])
    my_df = my_df.rename(columns={info_dict['FC cond x']: "FC cond x", info_dict['FC cond y']: "FC cond y",
                               info_dict['padj cond x']: "padj cond x", info_dict['padj cond y']: "padj cond y"})
    if info_dict['zscore']==True :
        my_df['FC cond x']=zscore(my_df['FC cond x'])
        my_df['FC cond y']=zscore(my_df['FC cond y'])

    if (my_df['FC cond x'].dtype != np.float64 or my_df['FC cond x'].dtype != np.int64):
        print('Data in column ',info_dict['FC cond x'],' contains not numeric data [ERROR!]')
        print('Attempt to remove not numeric rows')
        my_df=my_df.dropna(subset=['FC cond x','FC cond y','padj cond x','padj cond y'],how='any')
        my_df['FC cond x'] = pd.to_numeric(my_df['FC cond x'], errors='coerce')
    if (my_df['FC cond y'].dtype != np.float64 or my_df['FC cond y'].dtype != np.int64):
        print('Data in column ',info_dict['FC cond y'],' contains not numeric data [ERROR!]')
        print('Attempt to remove not numeric rows')
        my_df=my_df.dropna(subset=['FC cond x','FC cond y','padj cond x','padj cond y'],how='any')
        my_df['FC cond y'] = pd.to_numeric(my_df['FC cond y'], errors='coerce')
    if (my_df['padj cond x'].dtype != np.float64 or my_df['padj cond x'].dtype != np.int64):
        print('Data in column ',info_dict['FC cond x'],' contains not numeric data [ERROR!]')
        print('Attempt to remove not numeric rows')
        my_df=my_df.dropna(subset=['FC cond x','FC cond y','padj cond x','padj cond y'],how='any')
        my_df['padj cond x'] = pd.to_numeric(my_df['padj cond x'], errors='coerce')
    if (my_df['padj cond y'].dtype != np.float64 or my_df['padj cond y'].dtype != np.int64):
        print('Data in column ',info_dict['FC cond y'],' contains not numeric data [ERROR!]')
        print('Attempt to remove not numeric rows')
        my_df=my_df.dropna(subset=['FC cond x','FC cond y','padj cond x','padj cond y'],how='any')
        my_df['padj cond y'] = pd.to_numeric(my_df['padj cond y'], errors='coerce')
    return my_df
def scorecard_legend(info_dict3):
    '''
    This function creates a reference scorecard with colors and shaded areas.
    It could be used to guide data analysis. As input pass the parameters dictionary.
    '''
    mf=info_dict3['multiplication factor']
    save_folder1=info_dict3['save_dir']   
    if save_folder1[-1]!="/":
        save_folder1=save_folder1+"/"+info_dict3['Treatment1 name']+" "+info_dict3['Treatment2 name']+"/"
    else:
        save_folder1=save_folder1+info_dict3['Treatment1 name']+" "+info_dict3['Treatment2 name']+"/"
    if not isdir(save_folder1):
        makedirs(save_folder1)
    colori=info_dict3['colors']
    other_colori=info_dict3['other_colors']
    th_fold_change=info_dict3['th_fold_change']
    th_significance=info_dict3['th_significance']
    font_size1=info_dict3['font_size_examples']
    trasp=info_dict3['marker_trasp']
    trasp_rect=info_dict3['rect_trasp']
    col_rect=info_dict3['rect_colors']
    markers=info_dict3['markers']
    sizes= info_dict3['markers_sizes']
    use_notation=info_dict3['use_notation'],
    IS_EXAMPLE=info_dict3['is_example']
    fig_size=info_dict3['fig_size']
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    minimo,massimo=-fig_size,fig_size    
    ax.set_xlim(left=minimo,right=massimo)
    ax.set_ylim(bottom=minimo,top=massimo)
    ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)    
    ax.axhline(y=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    ax.axvline(x=0,color='k',linestyle='solid',lw=2.0)
    ax.axhline(y=0,color='k',linestyle='solid',lw=2.0)
     
    if use_notation:   
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict3['Treatment1 name']+" vs "+info_dict3['Control name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict3['Treatment2 name']+" vs "+info_dict3['Control name']+")")
    else:
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict3['Treatment1 name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict3['Treatment2 name']+")")        
    labels=[xc.upper() for xc in colori]
    
    ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change*mf), (massimo-th_fold_change*mf), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
    ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change*mf), (-massimo+th_fold_change*mf), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
    ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change*mf), (-massimo+th_fold_change*mf), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
    ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change*mf), (massimo-th_fold_change*mf), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))

    ax.add_patch(Rectangle((th_fold_change, th_fold_change*mf), (th_fold_change*mf-th_fold_change), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change), (massimo-th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change), (-massimo+th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change, -th_fold_change*mf), (th_fold_change*mf-th_fold_change), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change*mf,- th_fold_change), (massimo-th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change), (massimo-th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change, -th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change), (-massimo+th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))

    ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change), (massimo-th_fold_change*mf), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (th_fold_change*2), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    ax.add_patch(Rectangle((minimo,-th_fold_change), (-minimo-th_fold_change*mf), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    ax.add_patch(Rectangle((-th_fold_change,minimo), (th_fold_change*2), (-minimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))

    mid_green=(massimo+(th_fold_change*mf))/2
    if IS_EXAMPLE:
        ax.text(mid_green,mid_green, labels[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,-mid_green, labels[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(mid_green,-mid_green, labels[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,mid_green, labels[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )

        ax.text(mid_green,(th_fold_change*mf+th_fold_change)/2, labels[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[1] )
        ax.text(-mid_green,(th_fold_change*mf+th_fold_change)/2, labels[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[1] )
        ax.text(-mid_green,-(th_fold_change*mf+th_fold_change)/2, labels[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[1] )
        ax.text(mid_green,-(th_fold_change*mf+th_fold_change)/2, labels[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[1] )

        ax.text((th_fold_change*mf+th_fold_change)/2,-mid_green, labels[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[2],rotation=90 )
        ax.text((th_fold_change*mf+th_fold_change)/2,mid_green, labels[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[2],rotation=90 )
        ax.text(-(th_fold_change*mf+th_fold_change)/2,-mid_green, labels[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[2],rotation=90 )
        ax.text(-(th_fold_change*mf+th_fold_change)/2,mid_green, labels[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[2],rotation=90 )

        ax.text(mid_green,th_fold_change/2, labels[3]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(mid_green,-th_fold_change/2, labels[3]+'D (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,th_fold_change/2, labels[3]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,-th_fold_change/2, labels[3]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )

        ax.text(th_fold_change/2,-mid_green, labels[4]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.text(th_fold_change/2,mid_green, labels[4]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.text(-th_fold_change/2,-mid_green, labels[4]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.text(-th_fold_change/2,mid_green, labels[4]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.set_title('Regions of interest and gene color scheme for p<'+str(th_significance)+' diff. expr. entries')
    else:
        
        ax.text(mid_green,mid_green, 'A',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,-mid_green, 'A',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(mid_green,-mid_green, 'A',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,mid_green, 'A',size=font_size1, ha='center', va='center',color=colori[0]  )

        ax.text(mid_green,(th_fold_change*mf+th_fold_change)/2, 'B',size=font_size1, ha='center', va='center',color=colori[1] )
        ax.text(-mid_green,(th_fold_change*mf+th_fold_change)/2, 'B',size=font_size1, ha='center', va='center',color=colori[1] )
        ax.text(-mid_green,-(th_fold_change*mf+th_fold_change)/2, 'B',size=font_size1, ha='center', va='center',color=colori[1] )
        ax.text(mid_green,-(th_fold_change*mf+th_fold_change)/2, 'B',size=font_size1, ha='center', va='center',color=colori[1] )

        ax.text((th_fold_change*mf+th_fold_change)/2,-mid_green, 'C',size=font_size1, ha='center', va='center',color=colori[2],rotation=0 )
        ax.text((th_fold_change*mf+th_fold_change)/2,mid_green, 'C',size=font_size1, ha='center', va='center',color=colori[2],rotation=0 )
        ax.text(-(th_fold_change*mf+th_fold_change)/2,-mid_green, 'C',size=font_size1, ha='center', va='center',color=colori[2],rotation=0 )
        ax.text(-(th_fold_change*mf+th_fold_change)/2,mid_green, 'C',size=font_size1, ha='center', va='center',color=colori[2],rotation=0 )

        ax.text(mid_green,th_fold_change/2, 'D',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(mid_green,-th_fold_change/2, 'D',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,th_fold_change/2, 'D',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,-th_fold_change/2, 'D',size=font_size1, ha='center', va='center',color=colori[3]  )

        ax.text(th_fold_change/2,-mid_green, 'E',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.text(th_fold_change/2,mid_green, 'E',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.text(-th_fold_change/2,-mid_green, 'E',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.text(-th_fold_change/2,mid_green, 'E',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.set_title('Regions of interest and color scheme for p<'+str(th_significance)+' diff. expr. entries')
    if IS_EXAMPLE:
        plt.savefig(save_folder1+'EXAMPLE_colors.png',dpi=300,bbox_inches='tight')
    else:
        plt.savefig(save_folder1+'EXAMPLE_letters.png',dpi=300,bbox_inches='tight')
    plt.close()      
def scorecard(the_df,info_dict2):
    '''
    The DataFrame containing the CSV data will be used to build the four quadrants of the scorecard as single images (useful in case of
    a crowded graph) into a specific subfolder. Also the entries or genes on each region of interest will be saved on the hard disk as json
    file for further use.

    As input pass the Dataframe containing CSV data (previously loaded) and the dictionary of parameters.
    Repeat this function call for all experimental conditions (check examples on how to speed up and automate computations in case of
    multiple experimental conditions).
    '''
    IS_EXAMPLE=info_dict2['is_example']
    titolo=info_dict2['Scorecard title']
    fig_size=info_dict2['fig_size']
    save_folder=info_dict2['save_dir']    
    if save_folder[-1]!="/":
        save_folder=save_folder+"/"+info_dict2['Treatment1 name']+" "+info_dict2['Treatment2 name']+"/"
    else:
        save_folder=save_folder+info_dict2['Treatment1 name']+" "+info_dict2['Treatment2 name']+"/"
    if not isdir(save_folder):
        makedirs(save_folder)
    colori=info_dict2['colors']
    other_colori=info_dict2['other_colors']
    th_fold_change=info_dict2['th_fold_change']
    th_significance=info_dict2['th_significance']
    font_size1=info_dict2['font_size_quadrants']
    trasp=info_dict2['marker_trasp']
    trasp_rect=info_dict2['rect_trasp']
    col_rect=info_dict2['rect_colors']
    markers=info_dict2['markers']
    sizes= info_dict2['markers_sizes']
    gene_name= info_dict2['gene_name']    
    mf=info_dict2['multiplication factor']
    use_notation=info_dict2['use_notation']
    print('The dataset includes ',the_df.shape[0],' entries in total')
    labels=[xc.upper() for xc in colori]
    ###############################################################################################################################################
    ###############################################################################################################################################
    ###############################################################################################################################################
    fig, ax = plt.subplots(figsize=(info_dict2['fig_size'], info_dict2['fig_size']))
    texts1,texts2,texts3,texts4,texts5=[],[],[],[],[]
    p_val_x1,p_val_x2,p_val_x3,p_val_x4,p_val_x5=[],[],[],[],[]
    p_val_y1,p_val_y2,p_val_y3,p_val_y4,p_val_y5=[],[],[],[],[]
    
    common_up= the_df[(the_df['FC cond x'] > 0) & (the_df['FC cond y'] > 0) ] # FIRST QUADRANT 
    print('First quadrant will host ',common_up.shape[0],' entries')
    lista=common_up[gene_name].tolist()    
    for my_gene in lista:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y'] 
        if fch_x>th_fold_change*mf and fch_y>th_fold_change*mf:
            ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])
            if pval_x<th_significance and pval_y<th_significance:
                texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))
                p_val_x1.append(pval_x)
                p_val_y1.append(pval_y)
        elif (fch_x>th_fold_change and fch_x<=th_fold_change*mf) and (fch_y>th_fold_change and fch_y<=th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
        elif (fch_x>th_fold_change*mf) and (fch_y>th_fold_change and fch_y<=th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts2.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[1] ))
                p_val_x2.append(pval_x)
                p_val_y2.append(pval_y)
        elif (fch_x>th_fold_change and fch_x<=th_fold_change*mf) and (fch_y>th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts3.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[2]  ))
                p_val_x3.append(pval_x)
                p_val_y3.append(pval_y)
        elif (fch_x>th_fold_change*mf) and (fch_y<=th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                p_val_x4.append(pval_x)
                p_val_y4.append(pval_y)
        elif (fch_x<=th_fold_change ) and (fch_y>th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                p_val_x5.append(pval_x)
                p_val_y5.append(pval_y)
        elif (fch_x>th_fold_change and fch_x<th_fold_change*mf) and (fch_y<=th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])        
        elif (fch_x<=th_fold_change ) and (fch_y>th_fold_change and fch_y<th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
        else:
            ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        

    ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    #adjust_text(texts1, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts2, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts3, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts4, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts5, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    if IS_EXAMPLE:
        quadr[labels[0]]=[el.get_text() for el in texts1]
        quadr[labels[1]]=[el.get_text() for el in texts2]
        quadr[labels[2]]=[el.get_text() for el in texts3]
        quadr[labels[3]]=[el.get_text() for el in texts4]
        quadr[labels[4]]=[el.get_text() for el in texts5]
        quadr[labels[0]+'_x']=[el.get_position()[0] for el in texts1]
        quadr[labels[1]+'_x']=[el.get_position()[0] for el in texts2]
        quadr[labels[2]+'_x']=[el.get_position()[0] for el in texts3]
        quadr[labels[3]+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4]+'_x']=[el.get_position()[0] for el in texts5]
        quadr[labels[0]+'_y']=[el.get_position()[1] for el in texts1]
        quadr[labels[1]+'_y']=[el.get_position()[1] for el in texts2]
        quadr[labels[2]+'_y']=[el.get_position()[1] for el in texts3]
        quadr[labels[3]+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4]+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista
        quadr[labels[0]+'_pval_x']=p_val_x1
        quadr[labels[0]+'_pval_y']=p_val_y1
        quadr[labels[1]+'_pval_x']=p_val_x2
        quadr[labels[1]+'_pval_y']=p_val_y2
        quadr[labels[2]+'_pval_x']=p_val_x3
        quadr[labels[2]+'_pval_y']=p_val_y3
        quadr[labels[3]+'_pval_x']=p_val_x4
        quadr[labels[3]+'_pval_y']=p_val_y4
        quadr[labels[4]+'_pval_x']=p_val_x5
        quadr[labels[4]+'_pval_y']=p_val_y5
        
    else:
        quadr['A']=[el.get_text() for el in texts1]
        quadr['B']=[el.get_text() for el in texts2]
        quadr['C']=[el.get_text() for el in texts3]
        quadr['D']=[el.get_text() for el in texts4]
        quadr['E']=[el.get_text() for el in texts5]
        quadr['A_x']=[el.get_position()[0] for el in texts1]
        quadr['B_x']=[el.get_position()[0] for el in texts2]
        quadr['C_x']=[el.get_position()[0] for el in texts3]
        quadr['D_x']=[el.get_position()[0] for el in texts4]
        quadr['E_x']=[el.get_position()[0] for el in texts5]
        quadr['A_y']=[el.get_position()[1] for el in texts1]
        quadr['B_y']=[el.get_position()[1] for el in texts2]
        quadr['C_y']=[el.get_position()[1] for el in texts3]
        quadr['D_y']=[el.get_position()[1] for el in texts4]
        quadr['E_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista
        quadr['A_pval_x']=p_val_x1
        quadr['A_pval_y']=p_val_y1
        quadr['B_pval_x']=p_val_x2
        quadr['B_pval_y']=p_val_y2
        quadr['C_pval_x']=p_val_x3
        quadr['C_pval_y']=p_val_y3
        quadr['D_pval_x']=p_val_x4
        quadr['D_pval_y']=p_val_y4
        quadr['E_pval_x']=p_val_x5
        quadr['E_pval_y']=p_val_y5
    quadr['COLORS']=colori
    quadr['params']=info_dict2
    if use_notation:   
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+" vs "+info_dict2['Control name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment2 name']+" vs "+info_dict2['Control name']+")")
        ax.set_title('Both Exp. Conditions up-regulated ('+titolo+')')
    else:
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment2 name']+")") 
        ax.set_title('Both Exp. Conditions up-regulated')
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    
    ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change*mf), (xmax-th_fold_change*mf), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
    ax.add_patch(Rectangle((th_fold_change, th_fold_change*mf), (th_fold_change*mf-th_fold_change), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change), (xmax-th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change*mf, ymin), (xmax-th_fold_change*mf), (th_fold_change-ymin),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    ax.add_patch(Rectangle((xmin, th_fold_change*mf), (th_fold_change-xmin), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))

    plt.savefig(save_folder+'Quadrant1.png',dpi=300,bbox_inches='tight')
    plt.close()  
    with open(save_folder+'Quadrant1.json', 'w') as json_file:
        json.dump(quadr,json_file,  indent = 4)

    ###############################################################################################################################################
    ###############################################################################################################################################
    ###############################################################################################################################################
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    texts1,texts2,texts3,texts4,texts5=[],[],[],[],[]
    p_val_x1,p_val_x2,p_val_x3,p_val_x4,p_val_x5=[],[],[],[],[]
    p_val_y1,p_val_y2,p_val_y3,p_val_y4,p_val_y5=[],[],[],[],[]
    common_down=the_df[(the_df['FC cond x'] < 0) & (the_df['FC cond y'] < 0) ]
    print('Third quadrant will host ',common_down.shape[0],' entries')
    lista2=common_down[gene_name].tolist() 
    for my_gene in lista2:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y']

        if fch_x<-th_fold_change*mf and fch_y<-th_fold_change*mf:
            ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
            if pval_x<th_significance and pval_y<th_significance:
                texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                p_val_x1.append(pval_x)
                p_val_y1.append(pval_y)
        elif (fch_x<-th_fold_change and fch_x>=-th_fold_change*mf) and (fch_y<-th_fold_change and fch_y>=-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
        elif (fch_x<-th_fold_change*mf) and (fch_y<-th_fold_change and fch_y>=-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts2.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[1] ))
                p_val_x2.append(pval_x)
                p_val_y2.append(pval_y)
        elif (fch_x<-th_fold_change and fch_x>=-th_fold_change*mf) and (fch_y<-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])          
            if pval_x<th_significance and pval_y<th_significance:
                texts3.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[2]  ))
                p_val_x3.append(pval_x)
                p_val_y3.append(pval_y)
        elif (fch_x<-th_fold_change*mf) and (fch_y>=-th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                p_val_x4.append(pval_x)
                p_val_y4.append(pval_y)
        elif (fch_x>=-th_fold_change ) and (fch_y<-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                p_val_x5.append(pval_x)
                p_val_y5.append(pval_y)
        elif (fch_x<-th_fold_change and fch_x>-th_fold_change*mf) and (fch_y>=-th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])        
        elif (fch_x>=-th_fold_change ) and (fch_y<-th_fold_change and fch_y>-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
        else:
            ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        
        
    ax.axhline(y=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    #adjust_text(texts1, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts2, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts3, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts4, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts5, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    if IS_EXAMPLE:
        quadr[labels[0]]=[el.get_text() for el in texts1]
        quadr[labels[1]]=[el.get_text() for el in texts2]
        quadr[labels[2]]=[el.get_text() for el in texts3]
        quadr[labels[3]]=[el.get_text() for el in texts4]
        quadr[labels[4]]=[el.get_text() for el in texts5]
        quadr[labels[0]+'_x']=[el.get_position()[0] for el in texts1]
        quadr[labels[1]+'_x']=[el.get_position()[0] for el in texts2]
        quadr[labels[2]+'_x']=[el.get_position()[0] for el in texts3]
        quadr[labels[3]+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4]+'_x']=[el.get_position()[0] for el in texts5]
        quadr[labels[0]+'_y']=[el.get_position()[1] for el in texts1]
        quadr[labels[1]+'_y']=[el.get_position()[1] for el in texts2]
        quadr[labels[2]+'_y']=[el.get_position()[1] for el in texts3]
        quadr[labels[3]+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4]+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista2
        quadr[labels[0]+'_pval_x']=p_val_x1
        quadr[labels[0]+'_pval_y']=p_val_y1
        quadr[labels[1]+'_pval_x']=p_val_x2
        quadr[labels[1]+'_pval_y']=p_val_y2
        quadr[labels[2]+'_pval_x']=p_val_x3
        quadr[labels[2]+'_pval_y']=p_val_y3
        quadr[labels[3]+'_pval_x']=p_val_x4
        quadr[labels[3]+'_pval_y']=p_val_y4
        quadr[labels[4]+'_pval_x']=p_val_x5
        quadr[labels[4]+'_pval_y']=p_val_y5
        
    else:
        quadr['A']=[el.get_text() for el in texts1]
        quadr['B']=[el.get_text() for el in texts2]
        quadr['C']=[el.get_text() for el in texts3]
        quadr['D']=[el.get_text() for el in texts4]
        quadr['E']=[el.get_text() for el in texts5]
        quadr['A_x']=[el.get_position()[0] for el in texts1]
        quadr['B_x']=[el.get_position()[0] for el in texts2]
        quadr['C_x']=[el.get_position()[0] for el in texts3]
        quadr['D_x']=[el.get_position()[0] for el in texts4]
        quadr['E_x']=[el.get_position()[0] for el in texts5]
        quadr['A_y']=[el.get_position()[1] for el in texts1]
        quadr['B_y']=[el.get_position()[1] for el in texts2]
        quadr['C_y']=[el.get_position()[1] for el in texts3]
        quadr['D_y']=[el.get_position()[1] for el in texts4]
        quadr['E_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista2
        quadr['A_pval_x']=p_val_x1
        quadr['A_pval_y']=p_val_y1
        quadr['B_pval_x']=p_val_x2
        quadr['B_pval_y']=p_val_y2
        quadr['C_pval_x']=p_val_x3
        quadr['C_pval_y']=p_val_y3
        quadr['D_pval_x']=p_val_x4
        quadr['D_pval_y']=p_val_y4
        quadr['E_pval_x']=p_val_x5
        quadr['E_pval_y']=p_val_y5
    quadr['COLORS']=colori
    quadr['params']=info_dict2
    if use_notation:   
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+" vs "+info_dict2['Control name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment2 name']+" vs "+info_dict2['Control name']+")")
        ax.set_title('Both Exp. Conditions down-regulated ('+titolo+')')
    else:
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+")")   
        ax.set_title('Both Exp. Conditions down-regulated')
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change*mf), (xmin+th_fold_change*mf), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
    ax.add_patch(Rectangle((-th_fold_change, -th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change), (xmin+th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change*mf, ymax), (xmin+th_fold_change*mf), -(ymax+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    ax.add_patch(Rectangle((xmax, -th_fold_change*mf), -(xmax+th_fold_change), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))


    plt.savefig(save_folder+'Quadrant3.png',dpi=300,bbox_inches='tight')
    plt.close()  
    with open(save_folder+'Quadrant3.json', 'w') as json_file:
        json.dump(quadr,json_file,  indent = 4)
    ###############################################################################################################################################
    ###############################################################################################################################################
    ###############################################################################################################################################
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    texts1,texts2,texts3,texts4,texts5=[],[],[],[],[]
    p_val_x1,p_val_x2,p_val_x3,p_val_x4,p_val_x5=[],[],[],[],[]
    p_val_y1,p_val_y2,p_val_y3,p_val_y4,p_val_y5=[],[],[],[],[]
    common_1=the_df[(the_df['FC cond x'] < 0) & (the_df['FC cond y'] > 0) ] # Quadrant II
    print('Second quadrant will host ',common_1.shape[0],' entries')
    lista3=common_1[gene_name].tolist()
    for my_gene in lista3:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y']
        if fch_x<-th_fold_change*mf and fch_y>th_fold_change*mf:
            ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
            if pval_x<th_significance and pval_y<th_significance:
                texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                p_val_x1.append(pval_x)
                p_val_y1.append(pval_y)
        elif (fch_x<-th_fold_change and fch_x>=-th_fold_change*mf) and (fch_y>th_fold_change and fch_y<=th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
        elif (fch_x<-th_fold_change*mf) and (fch_y>th_fold_change and fch_y<=th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts2.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[1] ))
                p_val_x2.append(pval_x)
                p_val_y2.append(pval_y)
        elif (fch_x<-th_fold_change and fch_x>=-th_fold_change*mf) and (fch_y>th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])          
            if pval_x<th_significance and pval_y<th_significance:
                texts3.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[2]  ))
                p_val_x3.append(pval_x)
                p_val_y3.append(pval_y)
        elif (fch_x<-th_fold_change*mf) and (fch_y<=th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                p_val_x4.append(pval_x)
                p_val_y4.append(pval_y)
        elif (fch_x>=-th_fold_change ) and (fch_y>th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                p_val_x5.append(pval_x)
                p_val_y5.append(pval_y)
        elif (fch_x<-th_fold_change and fch_x>-th_fold_change*mf) and (fch_y<=th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])        
        elif (fch_x>=-th_fold_change ) and (fch_y>th_fold_change and fch_y<th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
        else:
            ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        
        
    ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    #adjust_text(texts1, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts2, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts3, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts4, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts5, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    
    if IS_EXAMPLE:
        quadr[labels[0]]=[el.get_text() for el in texts1]
        quadr[labels[1]]=[el.get_text() for el in texts2]
        quadr[labels[2]]=[el.get_text() for el in texts3]
        quadr[labels[3]]=[el.get_text() for el in texts4]
        quadr[labels[4]]=[el.get_text() for el in texts5]
        quadr[labels[0]+'_x']=[el.get_position()[0] for el in texts1]
        quadr[labels[1]+'_x']=[el.get_position()[0] for el in texts2]
        quadr[labels[2]+'_x']=[el.get_position()[0] for el in texts3]
        quadr[labels[3]+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4]+'_x']=[el.get_position()[0] for el in texts5]
        quadr[labels[0]+'_y']=[el.get_position()[1] for el in texts1]
        quadr[labels[1]+'_y']=[el.get_position()[1] for el in texts2]
        quadr[labels[2]+'_y']=[el.get_position()[1] for el in texts3]
        quadr[labels[3]+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4]+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista3
        quadr[labels[0]+'_pval_x']=p_val_x1
        quadr[labels[0]+'_pval_y']=p_val_y1
        quadr[labels[1]+'_pval_x']=p_val_x2
        quadr[labels[1]+'_pval_y']=p_val_y2
        quadr[labels[2]+'_pval_x']=p_val_x3
        quadr[labels[2]+'_pval_y']=p_val_y3
        quadr[labels[3]+'_pval_x']=p_val_x4
        quadr[labels[3]+'_pval_y']=p_val_y4
        quadr[labels[4]+'_pval_x']=p_val_x5
        quadr[labels[4]+'_pval_y']=p_val_y5
        
    else:
        quadr['A']=[el.get_text() for el in texts1]
        quadr['B']=[el.get_text() for el in texts2]
        quadr['C']=[el.get_text() for el in texts3]
        quadr['D']=[el.get_text() for el in texts4]
        quadr['E']=[el.get_text() for el in texts5]
        quadr['A_x']=[el.get_position()[0] for el in texts1]
        quadr['B_x']=[el.get_position()[0] for el in texts2]
        quadr['C_x']=[el.get_position()[0] for el in texts3]
        quadr['D_x']=[el.get_position()[0] for el in texts4]
        quadr['E_x']=[el.get_position()[0] for el in texts5]
        quadr['A_y']=[el.get_position()[1] for el in texts1]
        quadr['B_y']=[el.get_position()[1] for el in texts2]
        quadr['C_y']=[el.get_position()[1] for el in texts3]
        quadr['D_y']=[el.get_position()[1] for el in texts4]
        quadr['E_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista3
        quadr['A_pval_x']=p_val_x1
        quadr['A_pval_y']=p_val_y1
        quadr['B_pval_x']=p_val_x2
        quadr['B_pval_y']=p_val_y2
        quadr['C_pval_x']=p_val_x3
        quadr['C_pval_y']=p_val_y3
        quadr['D_pval_x']=p_val_x4
        quadr['D_pval_y']=p_val_y4
        quadr['E_pval_x']=p_val_x5
        quadr['E_pval_y']=p_val_y5
    quadr['COLORS']=colori
    quadr['params']=info_dict2
    if use_notation:   
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+" vs "+info_dict2['Control name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment2 name']+" vs "+info_dict2['Control name']+")")
        ax.set_title('Exp. Condition on X down-reg. and Exp. Condition on Y up-reg. ('+titolo+')')
    else:
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment2 name']+")")   
        ax.set_title(info_dict2['Treatment1 name']+' down-reg. and '+info_dict2['Treatment2 name']+' up-reg.')

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change*mf), (xmin+th_fold_change*mf), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
    ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change), (xmin+th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((-th_fold_change*mf, ymin), (xmin+th_fold_change*mf), (th_fold_change-ymin),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    ax.add_patch(Rectangle((xmax, th_fold_change*mf), -(xmax+th_fold_change), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    
    plt.savefig(save_folder+'Quadrant2.png',dpi=300,bbox_inches='tight')
    plt.close()  
    with open(save_folder+'Quadrant2.json', 'w') as json_file:
        json.dump(quadr,json_file,  indent = 4)
    ###############################################################################################################################################
    ###############################################################################################################################################
    ###############################################################################################################################################
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    texts1,texts2,texts3,texts4,texts5=[],[],[],[],[]
    p_val_x1,p_val_x2,p_val_x3,p_val_x4,p_val_x5=[],[],[],[],[]
    p_val_y1,p_val_y2,p_val_y3,p_val_y4,p_val_y5=[],[],[],[],[]
    common_2=the_df[(the_df['FC cond x'] > 0) & (the_df['FC cond y'] < 0) ] # Quadrant IV
    print('Forth quadrant will host ',common_2.shape[0],' entries')
    lista4=common_2[gene_name].tolist()
    for my_gene in lista4:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y']
        if fch_x>th_fold_change*mf and fch_y<-th_fold_change*mf:
            ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
            if pval_x<th_significance and pval_y<th_significance:
                texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                p_val_x1.append(pval_x)
                p_val_y1.append(pval_y)
        elif (fch_x>th_fold_change and fch_x<=th_fold_change*mf) and (fch_y<-th_fold_change and fch_y>=-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
        elif (fch_x>th_fold_change*mf) and (fch_y<-th_fold_change and fch_y>=-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts2.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[1] ))
                p_val_x2.append(pval_x)
                p_val_y2.append(pval_y)
        elif (fch_x>th_fold_change and fch_x<=th_fold_change*mf) and (fch_y<-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])          
            if pval_x<th_significance and pval_y<th_significance:
                texts3.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[2]  ))
                p_val_x3.append(pval_x)
                p_val_y3.append(pval_y)
        elif (fch_x>th_fold_change*mf) and (fch_y>=-th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                p_val_x4.append(pval_x)
                p_val_y4.append(pval_y)
        elif (fch_x<=th_fold_change ) and (fch_y<-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
            if pval_x<th_significance and pval_y<th_significance:
                texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                p_val_x5.append(pval_x)
                p_val_y5.append(pval_y)
        elif (fch_x>th_fold_change and fch_x<th_fold_change*mf) and (fch_y>=-th_fold_change) :
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])        
        elif (fch_x<=th_fold_change ) and (fch_y<-th_fold_change and fch_y>-th_fold_change*mf):
            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
        else:
            ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        

    ax.axhline(y=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    #adjust_text(texts1, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts2, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts3, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts4, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts5, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    if IS_EXAMPLE:
        quadr[labels[0]]=[el.get_text() for el in texts1]
        quadr[labels[1]]=[el.get_text() for el in texts2]
        quadr[labels[2]]=[el.get_text() for el in texts3]
        quadr[labels[3]]=[el.get_text() for el in texts4]
        quadr[labels[4]]=[el.get_text() for el in texts5]
        quadr[labels[0]+'_x']=[el.get_position()[0] for el in texts1]
        quadr[labels[1]+'_x']=[el.get_position()[0] for el in texts2]
        quadr[labels[2]+'_x']=[el.get_position()[0] for el in texts3]
        quadr[labels[3]+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4]+'_x']=[el.get_position()[0] for el in texts5]
        quadr[labels[0]+'_y']=[el.get_position()[1] for el in texts1]
        quadr[labels[1]+'_y']=[el.get_position()[1] for el in texts2]
        quadr[labels[2]+'_y']=[el.get_position()[1] for el in texts3]
        quadr[labels[3]+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4]+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista4
        quadr[labels[0]+'_pval_x']=p_val_x1
        quadr[labels[0]+'_pval_y']=p_val_y1
        quadr[labels[1]+'_pval_x']=p_val_x2
        quadr[labels[1]+'_pval_y']=p_val_y2
        quadr[labels[2]+'_pval_x']=p_val_x3
        quadr[labels[2]+'_pval_y']=p_val_y3
        quadr[labels[3]+'_pval_x']=p_val_x4
        quadr[labels[3]+'_pval_y']=p_val_y4
        quadr[labels[4]+'_pval_x']=p_val_x5
        quadr[labels[4]+'_pval_y']=p_val_y5
        
    else:
        quadr['A']=[el.get_text() for el in texts1]
        quadr['B']=[el.get_text() for el in texts2]
        quadr['C']=[el.get_text() for el in texts3]
        quadr['D']=[el.get_text() for el in texts4]
        quadr['E']=[el.get_text() for el in texts5]
        quadr['A_x']=[el.get_position()[0] for el in texts1]
        quadr['B_x']=[el.get_position()[0] for el in texts2]
        quadr['C_x']=[el.get_position()[0] for el in texts3]
        quadr['D_x']=[el.get_position()[0] for el in texts4]
        quadr['E_x']=[el.get_position()[0] for el in texts5]
        quadr['A_y']=[el.get_position()[1] for el in texts1]
        quadr['B_y']=[el.get_position()[1] for el in texts2]
        quadr['C_y']=[el.get_position()[1] for el in texts3]
        quadr['D_y']=[el.get_position()[1] for el in texts4]
        quadr['E_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista4
        quadr['A_pval_x']=p_val_x1
        quadr['A_pval_y']=p_val_y1
        quadr['B_pval_x']=p_val_x2
        quadr['B_pval_y']=p_val_y2
        quadr['C_pval_x']=p_val_x3
        quadr['C_pval_y']=p_val_y3
        quadr['D_pval_x']=p_val_x4
        quadr['D_pval_y']=p_val_y4
        quadr['E_pval_x']=p_val_x5
        quadr['E_pval_y']=p_val_y5
    quadr['COLORS']=colori
    quadr['params']=info_dict2
    if use_notation:   
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+" vs "+info_dict2['Control name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment2 name']+" vs "+info_dict2['Control name']+")")
        ax.set_title('Exp. Condition on X up-reg. and Exp. Condition on Y down-reg. ('+titolo+')')
    else:
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict2['Treatment1 name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict2['Treatment2 name']+")")   
        ax.set_title(info_dict2['Treatment1 name']+' up-reg. and '+info_dict2['Treatment2 name']+' down-reg.')

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change*mf), (xmax-th_fold_change*mf), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
    ax.add_patch(Rectangle((th_fold_change, -th_fold_change*mf), (th_fold_change*mf-th_fold_change), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change*mf,-th_fold_change), (xmax-th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
    ax.add_patch(Rectangle((th_fold_change*mf, ymax), (xmax-th_fold_change*mf), -(ymax+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    ax.add_patch(Rectangle((xmin, -th_fold_change*mf), (th_fold_change-xmin), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))

    plt.savefig(save_folder+'Quadrant4.png',dpi=300,bbox_inches='tight')
    plt.close()  
    with open(save_folder+'Quadrant4.json', 'w') as json_file:
        json.dump(quadr,json_file,  indent = 4)    

def reconstruct_scorecard(my_directory):
    '''
    The functions loads each Scorecard quandrant in memory and builds the Cartesian plane with the Scorecard for a global view of the dataset.
    It will be stored on the hard disk into the specified folder.

    As input pass a string with the folder address containing the subfolders with each experimental conditions.
    For example, only pass the address of the main_folder as string:

        main_folder
             |
             |
             |-------Exp. Comparison 1
             |
             |-------Exp. Comparison 2
             |
             |-------Exp. Comparison 3
    '''
    all_dir=[ f.path for f in scandir(my_directory) if f.is_dir() ]

    for the_folder in all_dir:
        nome=the_folder.split('/')[-1]
        print('Scorecard reconstruction on '+nome+' folder')
        results=[]
        quadr_list=[]
        my_data={}
        results += [each for each in listdir(the_folder) if each.endswith('.json')]
        for file in results:
            quadrante=file.split('_')[0]
            quadr_list.append(quadrante)
            with open(the_folder+'/'+file) as f:
                my_data[quadrante]=json.load(f)
        use_notation=my_data[quadrante]['params']['use_notation']
        if use_notation:
            trt1=my_data[quadrante]['params']['Treatment1 name']
            trt2=my_data[quadrante]['params']['Treatment2 name']
            ctrl=my_data[quadrante]['params']['Control name']
        else:
            trt1=my_data[quadrante]['params']['Treatment1 name']
            trt2=my_data[quadrante]['params']['Treatment2 name']            
        IS_EXAMPLE=my_data[quadrante]['params']['is_example']
        colori=my_data[quadrante]['params']['colors']
        if IS_EXAMPLE:
            etichette=[xc.upper() for xc in colori]
        else:
            etichette=['A','B','C','D','E']
        titolo=my_data[quadrante]['params']['Scorecard title']
        fig_size=my_data[quadrante]['params']['fig_size']
        save_folder=my_data[quadrante]['params']['save_dir']    
        if save_folder[-1]!="/":
            save_folder=save_folder+"/"+trt1+" "+trt2+"/"
        else:
            save_folder=save_folder+trt1+" "+trt2+"/"
        if not isdir(save_folder):
            print('Run scorecard quadrants creation before attempting reconstruction')
        th_fold_change=my_data[quadrante]['params']['th_fold_change']
        th_significance=my_data[quadrante]['params']['th_significance']
        font_size1=my_data[quadrante]['params']['font_size_quadrants']
        trasp=my_data[quadrante]['params']['marker_trasp']
        trasp_rect=my_data[quadrante]['params']['rect_trasp']
        col_rect=my_data[quadrante]['params']['rect_colors']
        markers=my_data[quadrante]['params']['markers']
        sizes= my_data[quadrante]['params']['markers_sizes']
        gene_name= my_data[quadrante]['params']['gene_name']        
        mf=my_data[quadrante]['params']['multiplication factor']
        
          
            
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        minimo,massimo=-fig_size,fig_size    
        ax.set_xlim(left=minimo,right=massimo)
        ax.set_ylim(bottom=minimo,top=massimo)
        ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
        ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)    
        ax.axhline(y=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axhline(y=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
        ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

        ax.axvline(x=0,color='k',linestyle='solid',lw=2.0)
        ax.axhline(y=0,color='k',linestyle='solid',lw=2.0)
         
        if use_notation:   
            ax.set_xlabel("$log_2$ Fold Change ("+trt1+" vs "+ctrl+")")
            ax.set_ylabel("$log_2$ Fold Change ("+trt2+" vs "+ctrl+")")
        else:
            ax.set_xlabel("$log_2$ Fold Change (Treatment 1 vs Control)")
            ax.set_ylabel("$log_2$ Fold Change (Treatment 2 vs Control)")        
        
        texts1,texts2,texts3,texts4,texts5=[],[],[],[],[]
        ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change*mf), (massimo-th_fold_change*mf), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change*mf), (-massimo+th_fold_change*mf), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change*mf), (-massimo+th_fold_change*mf), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change*mf), (massimo-th_fold_change*mf), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))

        ax.add_patch(Rectangle((th_fold_change, th_fold_change*mf), (th_fold_change*mf-th_fold_change), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change), (massimo-th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change), (-massimo+th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change, -th_fold_change*mf), (th_fold_change*mf-th_fold_change), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change*mf,- th_fold_change), (massimo-th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change), (massimo-th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change, -th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (-massimo+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change), (-massimo+th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))

        ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change), (massimo-th_fold_change*mf), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (th_fold_change*2), (massimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((minimo,-th_fold_change), (-minimo-th_fold_change*mf), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((-th_fold_change,minimo), (th_fold_change*2), (-minimo-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        for quadrante in quadr_list:
            for eti in etichette:
                lista_tmp=my_data[quadrante][eti]
                for idx_gene,my_gene in enumerate(lista_tmp):
                    fch_x=  my_data[quadrante][eti+'_x'][idx_gene]
                    fch_y=  my_data[quadrante][eti+'_y'][idx_gene]
                    if eti==etichette[0]:
                        ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])
                        texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))
                    elif eti==etichette[1]:
                        ax.scatter( fch_x,fch_y, facecolors = colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
                        texts2.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[1] ))
                    elif eti==etichette[2]:
                        ax.scatter( fch_x,fch_y, facecolors = colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
                        texts3.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[2]  ))
                    elif eti==etichette[3]:
                        ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                        texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                    elif eti==etichette[4]:
                        ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                        texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
        
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        if use_notation:   
            ax.set_xlabel("$log_2$ Fold Change ("+trt1+" vs "+ctrl+")")
            ax.set_ylabel("$log_2$ Fold Change ("+trt2+" vs "+ctrl+")")
            ax.set_title('Scorecard ('+titolo+')')
        else:
            ax.set_xlabel("$log_2$ Fold Change ("+trt1+")")
            ax.set_ylabel("$log_2$ Fold Change ("+trt2+")")   
            ax.set_title('Scorecard and regions of interest')
        plt.savefig(save_folder+'Scorecard.png',dpi=300,bbox_inches='tight')
        plt.close()

def calc_scarto(radians_r,valore_r):
    d_angle=radians_r* 180.0 / np.pi    
    return (d_angle+valore_r)* np.pi / 180.0    
def multiple_view(my_directory,scarto=3,marker_size=100): # "scarto" adjusts the jitter of the points +/- the radial axes
    '''
    The function will generate a view of all experimental comparisons, highlighting genes or entries beloning to the regions
    of interest of the scorecard. Pass as input the string of the main_folder containing all the subfolder with the experimental
    tests previously created by the scorecard function.

            main_folder
             |
             |
             |-------Exp. Comparison 1
             |
             |-------Exp. Comparison 2
             |
             |-------Exp. Comparison 3
    '''
    if my_directory[-1]!="/":
        my_directory=my_directory+"/"
    GREY_LIGHT = "#f2efe8"    
    all_dir=[ f.path for f in scandir(my_directory) if f.is_dir() ]
    all_data={}    
    VARIABLES=[]
    for the_folder in all_dir:
        nome=the_folder.split('/')[-1]
        results=[]
        quadr_list=[]
        my_data={}
        results += [each for each in listdir(the_folder) if each.endswith('.json')]
        for file in results:
            quadrante=file.split('.')[0]
            quadr_list.append(quadrante)
            with open(the_folder+'/'+file) as f:
                my_data[quadrante]=json.load(f)
        all_data[nome]=my_data
        VARIABLES.append(nome)
    VARIABLES_N = len(VARIABLES)
    ANGLES = [n / VARIABLES_N * 2 * np.pi for n in range(VARIABLES_N)]
    ANGLES += ANGLES[:1]
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, polar=True)
    #fig.patch.set_facecolor(BG_WHITE)
    #ax.set_facecolor(BG_WHITE)
    all_x=[]
    all_y=[]
    VAR_NAMES=[]    
    for ind_i,i_folder in enumerate(VARIABLES):
        nomi=i_folder.split(' ')
        VAR_NAMES.append(nomi[1]+'   '+nomi[0])
        tmp=all_data[i_folder]
        for i_key, key_name in enumerate(['Quadrant1', 'Quadrant2', 'Quadrant3', 'Quadrant4']):# Quadrants
            tmp2=tmp[key_name]
            all_gruppo=[*tmp2]
            etichette=[xc.upper() for xc in tmp2['COLORS']]
            colore=tmp2['COLORS']
            for i_gruppo, gruppo in enumerate(all_gruppo):                
                if gruppo in ['A', 'B', 'C', 'D', 'E']:                    
                    valore1=tmp2[gruppo+'_x']
                    valore2=tmp2[gruppo+'_y']                    
                    for x,y in zip(valore1,valore2):
                        ax.scatter(calc_scarto(ANGLES[ind_i],-scarto), x, s=marker_size, c=colore[i_gruppo], zorder=10)
                        ax.scatter(calc_scarto(ANGLES[ind_i],scarto), y, s=marker_size, c=colore[i_gruppo], zorder=10)
                        ax.plot([calc_scarto(ANGLES[ind_i],-scarto),calc_scarto(ANGLES[ind_i],scarto)], [x,y], c=colore[i_gruppo], linewidth=0.5, label=gruppo)
                        all_x.append(x)
                        all_y.append(y)
                elif gruppo in etichette:                    
                    valore1=tmp2[gruppo+'_x']
                    valore2=tmp2[gruppo+'_y']
                    for x,y in zip(valore1,valore2):
                        ax.scatter(calc_scarto(ANGLES[ind_i],-scarto), x, s=marker_size, c=colore[i_gruppo], zorder=10)
                        ax.scatter(calc_scarto(ANGLES[ind_i],scarto), y, s=marker_size, c=colore[i_gruppo], zorder=10)
                        ax.plot([calc_scarto(ANGLES[ind_i],-scarto),calc_scarto(ANGLES[ind_i],scarto)], [x,y], c=colore[i_gruppo], linewidth=0.5, label=gruppo)
                        all_x.append(x)
                        all_y.append(y)
    labels = []
    r_min,r_max=np.amin(all_x+all_y)-0.5,np.amax(all_x+all_y)+0.5    
    ax.set_ylim(r_min,r_max)
    plt.title('Exp. cond. comparison', size=20, y=1.05)
    ax.set_xticks(ANGLES[:-1])
    ax.set_xticklabels(VAR_NAMES, size=14)
    angles = np.rad2deg(ANGLES[:-1])-90
    labels = []
    for label, angle in zip(ax.get_xticklabels(), angles):
        x,y = label.get_position()
        l_txt=label.get_text()
        n=len(l_txt)
        if n>6:
            if n%2 == 0:
              s1 = slice(0,n//2)
              s2 = slice(n//2,n)
            else:
              s1 = slice(0,n//2)
              s2 = slice(n//2,n)
            lab = ax.text(x,y, chr(8592)+l_txt[s1]+"\n"+l_txt[s2]+chr(8594), transform=label.get_transform(),ha=label.get_ha(), va=label.get_va(),fontsize=10)
        else    :
            lab = ax.text(x,y, l_txt, transform=label.get_transform(),ha=label.get_ha(), va=label.get_va(),fontsize=16)
        lab.set_rotation(angle)
        labels.append(lab)
    ax.set_xticklabels([])
    HANGLES = np.linspace(0, 2 * np.pi)
    H0 = np.zeros(len(HANGLES))
    H1 = np.ones(len(HANGLES)) * -1
    H2 = np.ones(len(HANGLES))
    ax.fill(HANGLES, H0, GREY_LIGHT)

    plt.savefig(my_directory+'Exp_Comp.png',dpi=300,bbox_inches='tight')
    plt.close()

def make_volcano(my_directory):
    '''
    This routine saves standard Volcano plots of all the experimental comparisons carried out by the scorecard function.
    The volcano plots will only contain the genes or entries belonging to the regions of interest of the scorecard.
    The regions of interest of the scorecard will be shaded also on the Volcano plots.

    As input pass the main_folder containing all the subfolders with the experimental comparisons created by the scorecard function call.

            main_folder
             |
             |
             |-------Exp. Comparison 1
             |
             |-------Exp. Comparison 2
             |
             |-------Exp. Comparison 3
    
    '''
    all_dir=[ f.path for f in scandir(my_directory) if f.is_dir() ]

    for the_folder in all_dir:
        nome=the_folder.split('/')[-1]
        print('Scorecard reconstruction on '+nome+' folder')
        results=[]
        quadr_list=[]
        my_data={}
        results += [each for each in listdir(the_folder) if each.endswith('.json')]
        for file in results:
            quadrante=file.split('_')[0]
            quadr_list.append(quadrante)
            with open(the_folder+'/'+file) as f:
                my_data[quadrante]=json.load(f)
        epi=my_data[quadrante]['params']['Log Epsilon']
        ctrl=my_data[quadrante]['params']['Control name']
        IS_EXAMPLE=my_data[quadrante]['params']['is_example']
        colori=my_data[quadrante]['params']['colors']
        other_colori=my_data[quadrante]['params']['other_colors']
        if IS_EXAMPLE:
            etichette=[xc.upper() for xc in colori]
        else:
            etichette=['A','B','C','D','E']
        titolo=my_data[quadrante]['params']['Scorecard title']
        fig_size=my_data[quadrante]['params']['fig_size']
        save_folder=my_data[quadrante]['params']['save_dir']
        trt1=my_data[quadrante]['params']['Treatment1 name']
        trt2=my_data[quadrante]['params']['Treatment2 name']
        if save_folder[-1]!="/":
            save_folder=save_folder+"/"+trt1+" "+trt2+"/"
        else:
            save_folder=save_folder+trt1+" "+trt2+"/"
        if not isdir(save_folder):
            print('Run scorecard quadrants creation before attempting reconstruction')
        th_fold_change=my_data[quadrante]['params']['th_fold_change']
        th_significance=my_data[quadrante]['params']['th_significance']
        font_size1=my_data[quadrante]['params']['font_size_quadrants']
        trasp=my_data[quadrante]['params']['marker_trasp']
        trasp_rect=my_data[quadrante]['params']['rect_trasp']
        col_rect=my_data[quadrante]['params']['rect_colors']
        markers=my_data[quadrante]['params']['markers']
        sizes= my_data[quadrante]['params']['markers_sizes']
        gene_name= my_data[quadrante]['params']['gene_name']
        
        mf=my_data[quadrante]['params']['multiplication factor']
        use_notation=my_data[quadrante]['params']['use_notation']

            
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(18, 8),sharey=True)
        texts1x,texts2x,texts3x,texts4x,texts5x=[],[],[],[],[]
        texts1y,texts2y,texts3y,texts4y,texts5y=[],[],[],[],[]
        for quadrante in quadr_list:
            for eti in etichette:
                lista_tmp=my_data[quadrante][eti]
                if len(lista_tmp)>0:
                    for idx_gene,my_gene in enumerate(lista_tmp):
                        fch_x=  my_data[quadrante][eti+'_x'][idx_gene]
                        fch_y=  my_data[quadrante][eti+'_y'][idx_gene]
                        pval_x=  -log10_with_epsilon(my_data[quadrante][eti+'_pval_x'][idx_gene],epi)
                        pval_y=  -log10_with_epsilon(my_data[quadrante][eti+'_pval_y'][idx_gene],epi)
                        if eti==etichette[0] :
                            ax1.scatter( fch_x,pval_x, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])
                            texts1x.append( ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))
                            ax2.scatter( fch_y,pval_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])
                            texts1y.append( ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))
                        elif eti==etichette[1]:
                            ax1.scatter( fch_x,pval_x, facecolors = colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
                            texts2x.append(ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=colori[1] ))
                            ax2.scatter( fch_y,pval_y, facecolors = colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
                            texts2y.append(ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=colori[1] ))
                        elif eti==etichette[2]:
                            ax1.scatter( fch_x,pval_x, facecolors = colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
                            texts3x.append(ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=colori[2]  ))
                            ax2.scatter( fch_y,pval_y, facecolors = colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[1])
                            texts3y.append(ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=colori[2]  ))
                        elif eti==etichette[3]:
                            ax1.scatter( fch_x,pval_x, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                            texts4x.append(ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                            ax2.scatter( fch_y,pval_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                            texts4y.append(ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                        elif eti==etichette[4]:
                            ax1.scatter( fch_x,pval_x, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                            texts5x.append(ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                            ax2.scatter( fch_y,pval_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                            texts5y.append(ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
        out_x=flatten([texts1x,texts2x,texts3x,texts4x,texts5x])
        out_y=flatten([texts1y,texts2y,texts3y,texts4y,texts5y])
        if len(out_x)>0:
            adjust_text(out_x, ax=ax1,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        if len(out_y)>0:
            adjust_text(out_y, ax=ax2,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

        for ax,tmp_str in zip([ax1,ax2],[trt1,trt2]):    
            
            ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
            ax.axvline(x=0,color='grey',linestyle='dotted',lw=0.5)
            ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

            ax.axhline(y=-np.log10(th_significance),color='grey',linestyle='dashdot',lw=1.5)
            #ax.axhline(y=-np.log10(th_significance/mf),color='grey',linestyle='dashdot',lw=1.0)

            ax.set_xlabel("log2 Fold Change ("+tmp_str+" vs "+ctrl+")")
            ax.set_ylabel("-log10 Adjusted p-value")
            ax.set_title('Data: '+tmp_str+' versus '+ctrl)

            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()

            ax.add_patch(Rectangle((xmin, -np.log10(th_significance)), -(xmin+th_fold_change*mf), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
            ax.add_patch(Rectangle((th_fold_change*mf, -np.log10(th_significance)), (xmax-th_fold_change*mf), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2])                         )
            ax.add_patch(Rectangle((-th_fold_change*mf, -np.log10(th_significance)), -(-th_fold_change*mf+th_fold_change), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
            ax.add_patch(Rectangle((th_fold_change, -np.log10(th_significance)), (th_fold_change*mf-th_fold_change), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
            ax.add_patch(Rectangle((-th_fold_change, -np.log10(th_significance)), (th_fold_change*2), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
            

        if the_folder[-1]!="/":
            the_folder=the_folder+"/"
        plt.savefig(the_folder+'Volcano.png',dpi=300,bbox_inches='tight')
        plt.close()  
def multiple_bars(my_directory,height=0.4, try_adj_test=False,text_adj_x=0.1,text_adj_y=0.6):
    '''
    The graph created by this function shows each gene or entry previously extracted with the Scorecard, displaying the Log2 Fold Change
    as bar plots. Paired bars depict the experimental conditions being compared. Colors reflect the color scheme employed in the Scorecard.

    As input pass the string containing the main_folder hosting all subfolders with each experimental comparison created by the Scorecard function.

            main_folder
             |
             |
             |-------Exp. Comparison 1
             |
             |-------Exp. Comparison 2
             |
             |-------Exp. Comparison 3
    '''
    all_dir=[ f.path for f in scandir(my_directory) if f.is_dir() ]
    for the_folder in all_dir:
        nome=the_folder.split('/')[-1]
        print('Barplot on '+nome+' folder')
        results=[]
        quadr_list=[]
        my_data={}
        results += [each for each in listdir(the_folder) if each.endswith('.json')]
        for file in results:
            quadrante=file.split('_')[0]
            quadr_list.append(quadrante)
            with open(the_folder+'/'+file) as f:
                my_data[quadrante]=json.load(f)        
        ctrl=my_data[quadrante]['params']['Control name']
        IS_EXAMPLE=my_data[quadrante]['params']['is_example']
        colori=my_data[quadrante]['params']['colors']
        other_colori=my_data[quadrante]['params']['other_colors']
        if IS_EXAMPLE:
            etichette=[xc.upper() for xc in colori]
        else:
            etichette=['A','B','C','D','E']
        titolo=my_data[quadrante]['params']['Scorecard title']
        fig_size=my_data[quadrante]['params']['fig_size']
        save_folder=my_data[quadrante]['params']['save_dir']
        trt1=my_data[quadrante]['params']['Treatment1 name']
        trt2=my_data[quadrante]['params']['Treatment2 name']
        if save_folder[-1]!="/":
            save_folder=save_folder+"/"+trt1+" "+trt2+"/"
        else:
            save_folder=save_folder+trt1+" "+trt2+"/"
        if not isdir(save_folder):
            print('Run scorecard quadrants creation before attempting reconstruction')
        th_fold_change=my_data[quadrante]['params']['th_fold_change']
        th_significance=my_data[quadrante]['params']['th_significance']
        font_size1=my_data[quadrante]['params']['font_size_quadrants']
        trasp=my_data[quadrante]['params']['marker_trasp']
        trasp_rect=my_data[quadrante]['params']['rect_trasp']
        col_rect=my_data[quadrante]['params']['rect_colors']
        markers=my_data[quadrante]['params']['markers']
        sizes= my_data[quadrante]['params']['markers_sizes']
        gene_name= my_data[quadrante]['params']['gene_name']
        offset=height/2
        mf=my_data[quadrante]['params']['multiplication factor']
        use_notation=my_data[quadrante]['params']['use_notation']
        Position=0
        fig = plt.figure()
        ax = fig.add_subplot(111)
        all_x,all_y,all_genes,all_colors,labels_x,labels_y=[],[],[],[],[],[]
        for quadrante in quadr_list:
            for eti in etichette:
                lista_tmp=my_data[quadrante][eti]
                if len(lista_tmp)>0:                    
                    for idx_gene,my_gene in enumerate(lista_tmp):
                        fch_x= round( my_data[quadrante][eti+'_x'][idx_gene],2)
                        fch_y= round( my_data[quadrante][eti+'_y'][idx_gene],2)
                        if eti==etichette[0] :
                            all_x.append( fch_x)
                            all_y.append( fch_y)
                            all_genes.append( my_gene)
                            all_colors.append(colori[0])
                            Position = Position+1
                            labels_x.append(trt1)
                            labels_y.append(trt2)
                        elif eti==etichette[1]:
                            all_x.append( fch_x)
                            all_y.append( fch_y)
                            all_genes.append( my_gene)
                            all_colors.append(colori[1])
                            Position = Position+1
                            labels_x.append(trt1)
                            labels_y.append(trt2) 
                        elif eti==etichette[2]:
                            all_x.append( fch_x)
                            all_y.append( fch_y)
                            all_genes.append( my_gene)
                            all_colors.append(colori[2])
                            Position = Position+1
                            labels_x.append(trt1)
                            labels_y.append(trt2)
                        elif eti==etichette[3]:
                            all_x.append( fch_x)
                            all_y.append( fch_y)
                            all_genes.append( my_gene)
                            all_colors.append(colori[3])
                            Position = Position+1
                            labels_x.append(trt1)
                            labels_y.append(trt2)
                        elif eti==etichette[4]:
                            all_x.append( fch_x)
                            all_y.append( fch_y)
                            all_genes.append( my_gene)
                            all_colors.append(colori[4])
                            Position = Position+1
                            labels_x.append(trt1)
                            labels_y.append(trt2)
        str_x=[str(x) for x in all_x]
        str_y=[str(x) for x in all_y]
        lx=list(map(' '.join, zip(labels_x, str_x)))
        ly=list(map(' '.join, zip(labels_y, str_y)))
        if Position>0:
            y_pos = np.arange(Position)
            if try_adj_test:
                bars1 =ax.barh( y_pos-offset,all_x, height=height, color= all_colors, align='center', label=lx)
                bars2 =ax.barh( y_pos+offset,all_y, height=height, color= all_colors, align='center', label=ly)
                text1=[]
                text2=[]
            else:
                ax.barh( y_pos-offset,all_x, height=height, color= all_colors, align='center', label=lx)
                ax.barh( y_pos+offset,all_y, height=height, color= all_colors, align='center', label=ly)
            for i, v in enumerate(all_x):
                if v>0:
                    if try_adj_test:
                        text1.append(ax.text(v +text_adj_x, y_pos[i]-offset, lx[i], color='k', fontsize=4, verticalalignment='center'))
                    else:
                        ax.text(v +text_adj_x, y_pos[i]-offset, lx[i], color='k', fontsize=4, verticalalignment='center')
                else:
                    if try_adj_test:
                        text1.append(ax.text(v -text_adj_y, y_pos[i]-offset, lx[i], color='k', fontsize=4, verticalalignment='center'))
                    else:
                        ax.text(v -text_adj_y, y_pos[i]-offset, lx[i], color='k', fontsize=4, verticalalignment='center')
            for i, v in enumerate(all_y):
                if v>0:
                    if try_adj_test:
                        text2.append(ax.text(v + text_adj_x, y_pos[i]+offset, ly[i], color='k', fontsize=4, verticalalignment='center'))
                    else:
                        ax.text(v + text_adj_x, y_pos[i]+offset, ly[i], color='k', fontsize=4, verticalalignment='center')
                else:
                    if try_adj_test:
                        text2.append( ax.text(v - text_adj_y, y_pos[i]+offset, ly[i], color='k', fontsize=4, verticalalignment='center'))
                    else:
                        ax.text(v - text_adj_y, y_pos[i]+offset, ly[i], color='k', fontsize=4, verticalalignment='center')
            ax.set_yticks(y_pos, labels=all_genes)
            ax.invert_yaxis()  # labels read top-to-bottom
            ax.set_xlabel('$log_2$ Fold Change',fontsize=11)
            ax.set_title(titolo)
            if try_adj_test:
                adjust_text(flatten([text1,text2]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
            ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
            ax.axvline(x=0,color='grey',linestyle='dotted',lw=0.5)
            ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
            ax.tick_params(axis='y', labelsize=5)
            if the_folder[-1]!="/":
                the_folder=the_folder+"/"
            plt.savefig(the_folder+'Bars.png',dpi=300,bbox_inches='tight')
            plt.close()