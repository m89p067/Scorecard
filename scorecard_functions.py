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
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
from collections import Counter
import matplotlib

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
    # General information, assuming CSV file contains log2 fold change already computed
    the_dict['FC cond x']='log2 fold change T1 vs T0' # replace with exact column name containing log2 fold change Treatment1 vs Control
    the_dict['FC cond y']='log2 fold change T2 vs T0' # replace with exact column name containing log2 fold change Treatment1 vs Control
    the_dict['padj cond x']='padj T1 vs T0' # replace with exact column name containing adj p-values Treatment1 vs Control
    the_dict['padj cond y']='padj T2 vs T0' # replace with exact column name containing adj p-values Treatment1 vs Control
    the_dict['incl aver']=False # Include intermediate regions (between F.C. thresholds)    
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
    the_dict['other_colors']=['orange','gold','plum','lightgrey'] # colors of points outside the regions of interest
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
    the_dict['use_notation']=True # Use custom names inserted in 'Treatment1 name','Treatment2 name','Control name' when plotting axes labels
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
    if info_dict['multiplication factor']<1:
        print('ERROR: Change multiplication factor to a number >=1')
    elif info_dict['multiplication factor']==1 and info_dict['incl aver']==False:
        print('This will prodice a four-way plot, not the scorecard')
        info_dict['note']='Four Way Plot'
    elif info_dict['multiplication factor']==1 and info_dict['incl aver']==True:
        print('ERROR: Mismatch in selection, place incl_aver as FALSE otherwise multiplication factor will be raised automatically')
        info_dict['multiplication factor']=1.5
    return my_df
def scorecard_legend(info_dict3):
    '''
    This function creates a reference scorecard with colors and shaded areas.
    It could be used to guide data analysis. As input pass the parameters dictionary.
    '''
    mf=info_dict3['multiplication factor']
    if mf<1:
        print('ERROR: Change multiplication factor to a number >=1')

    elif mf==1:
        print('This will prodice a four-way plot, not the scorecard')
        
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
    incl_ave=info_dict3['incl aver']
    
    ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)    
    ax.axhline(y=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    if mf > 1:
        ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)        
        ax.axhline(y=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axvline(x=0,color='k',linestyle='solid',lw=2.0)
    ax.axhline(y=0,color='k',linestyle='solid',lw=2.0)
     
    if use_notation:   
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict3['Treatment1 name']+" vs "+info_dict3['Control name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict3['Treatment2 name']+" vs "+info_dict3['Control name']+")")
    else:
        ax.set_xlabel("$log_2$ Fold Change ("+info_dict3['Treatment1 name']+")")
        ax.set_ylabel("$log_2$ Fold Change ("+info_dict3['Treatment2 name']+")")
    
    labels=[xc.upper() for xc in colori]
    if incl_ave:        
        labels_ave=[xc.upper() for xc in other_colori]
    if mf > 1:
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
    elif mf==1:
        ax.add_patch(Rectangle((th_fold_change, th_fold_change), (massimo-th_fold_change), (massimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change, -th_fold_change), (-massimo+th_fold_change), (-massimo+th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change, th_fold_change), (-massimo+th_fold_change), (massimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((th_fold_change, -th_fold_change), (massimo-th_fold_change), (-massimo+th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((th_fold_change, -th_fold_change), (massimo-th_fold_change), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((-th_fold_change, th_fold_change), (th_fold_change*2), (massimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((minimo,-th_fold_change), (-minimo-th_fold_change), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((-th_fold_change,minimo), (th_fold_change*2), (-minimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))        
    mid_green=(massimo+(th_fold_change*mf))/2
    if IS_EXAMPLE==True and mf > 1:
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
        if incl_ave:
            ax.text((th_fold_change*mf+th_fold_change)/2,(th_fold_change*mf+th_fold_change)/2, labels_ave[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[0] )
            ax.text(-(th_fold_change*mf+th_fold_change)/2,-(th_fold_change*mf+th_fold_change)/2, labels_ave[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[0] )
            ax.text(-(th_fold_change*mf+th_fold_change)/2,(th_fold_change*mf+th_fold_change)/2, labels_ave[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[0] )
            ax.text((th_fold_change*mf+th_fold_change)/2,-(th_fold_change*mf+th_fold_change)/2, labels_ave[0]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[0] )

            ax.text(th_fold_change/2,(th_fold_change*mf+th_fold_change)/2, labels_ave[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[1] )
            ax.text(-th_fold_change/2,(th_fold_change*mf+th_fold_change)/2, labels_ave[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[1] )
            ax.text(th_fold_change/2,-(th_fold_change*mf+th_fold_change)/2, labels_ave[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[1] )
            ax.text(-th_fold_change/2,-(th_fold_change*mf+th_fold_change)/2, labels_ave[1]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[1] )

            ax.text((th_fold_change*mf+th_fold_change)/2,th_fold_change/2, labels_ave[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[2]  )
            ax.text((th_fold_change*mf+th_fold_change)/2,-th_fold_change/2, labels_ave[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[2]  )
            ax.text(-(th_fold_change*mf+th_fold_change)/2, th_fold_change/2,labels_ave[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[2]  )
            ax.text(-(th_fold_change*mf+th_fold_change)/2,-th_fold_change/2, labels_ave[2]+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=other_colori[2] )
    elif IS_EXAMPLE==False and mf > 1:
        
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
        if incl_ave:
            ax.text((th_fold_change*mf+th_fold_change)/2,(th_fold_change*mf+th_fold_change)/2, 'M',size=font_size1, ha='center', va='center',color=other_colori[0] )
            ax.text(-(th_fold_change*mf+th_fold_change)/2,-(th_fold_change*mf+th_fold_change)/2,'M',size=font_size1, ha='center', va='center',color=other_colori[0] )
            ax.text(-(th_fold_change*mf+th_fold_change)/2,(th_fold_change*mf+th_fold_change)/2, 'M',size=font_size1, ha='center', va='center',color=other_colori[0] )
            ax.text((th_fold_change*mf+th_fold_change)/2,-(th_fold_change*mf+th_fold_change)/2,'M',size=font_size1, ha='center', va='center',color=other_colori[0] )

            ax.text(th_fold_change/2,(th_fold_change*mf+th_fold_change)/2, 'R',size=font_size1, ha='center', va='center',color=other_colori[1] )
            ax.text(-th_fold_change/2,(th_fold_change*mf+th_fold_change)/2, 'R',size=font_size1, ha='center', va='center',color=other_colori[1] )
            ax.text(th_fold_change/2,-(th_fold_change*mf+th_fold_change)/2, 'R',size=font_size1, ha='center', va='center',color=other_colori[1] )
            ax.text(-th_fold_change/2,-(th_fold_change*mf+th_fold_change)/2,'R',size=font_size1, ha='center', va='center',color=other_colori[1] )

            ax.text((th_fold_change*mf+th_fold_change)/2,th_fold_change/2, 'S',size=font_size1, ha='center', va='center',color=other_colori[2] )
            ax.text((th_fold_change*mf+th_fold_change)/2,-th_fold_change/2, 'S',size=font_size1, ha='center', va='center',color=other_colori[2] )
            ax.text(-(th_fold_change*mf+th_fold_change)/2, th_fold_change/2,'S',size=font_size1, ha='center', va='center',color=other_colori[2] )
            ax.text(-(th_fold_change*mf+th_fold_change)/2,-th_fold_change/2, 'S',size=font_size1, ha='center', va='center',color=other_colori[2] )
    elif IS_EXAMPLE==True and mf == 1:
        ax.text(mid_green,mid_green, labels[0].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,-mid_green, labels[0].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(mid_green,-mid_green, labels[0].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,mid_green, labels[0].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(mid_green,th_fold_change/2, labels[3].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(mid_green,-th_fold_change/2, labels[3].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,th_fold_change/2, labels[3].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,-th_fold_change/2, labels[3].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[3]  )

        ax.text(th_fold_change/2,-mid_green, labels[4].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.text(th_fold_change/2,mid_green, labels[4].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.text(-th_fold_change/2,-mid_green, labels[4].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.text(-th_fold_change/2,mid_green, labels[4].lower()+' (p<'+str(th_significance)+')',size=font_size1, ha='center', va='center',color=colori[4],rotation=90 )
        ax.set_title('Four-way plot (p<'+str(th_significance)+')')
    elif IS_EXAMPLE==False and mf == 1:
        ax.text(mid_green,mid_green, 'a',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,-mid_green, 'a',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(mid_green,-mid_green, 'a',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(-mid_green,mid_green, 'a',size=font_size1, ha='center', va='center',color=colori[0]  )
        ax.text(mid_green,th_fold_change/2, 'd',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(mid_green,-th_fold_change/2, 'd',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,th_fold_change/2, 'd',size=font_size1, ha='center', va='center',color=colori[3]  )
        ax.text(-mid_green,-th_fold_change/2, 'd',size=font_size1, ha='center', va='center',color=colori[3]  )

        ax.text(th_fold_change/2,-mid_green, 'e',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.text(th_fold_change/2,mid_green, 'e',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.text(-th_fold_change/2,-mid_green, 'e',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.text(-th_fold_change/2,mid_green, 'e',size=font_size1, ha='center', va='center',color=colori[4],rotation=0 )
        ax.set_title('Four-way plot (p<'+str(th_significance)+')')
    if IS_EXAMPLE==True and mf>1:
        plt.savefig(save_folder1+'SCORECARD_colors.png',dpi=300,bbox_inches='tight')
    elif IS_EXAMPLE==False and mf>1:
        plt.savefig(save_folder1+'SCORECARD_letters.png',dpi=300,bbox_inches='tight')
    elif IS_EXAMPLE==True and mf==1:
        plt.savefig(save_folder1+'FOURWAYPLOT_colors.png',dpi=300,bbox_inches='tight')
    elif IS_EXAMPLE==False and mf==1:
        plt.savefig(save_folder1+'FOURWAYPLOT_letters.png',dpi=300,bbox_inches='tight')
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
    if mf<1:
        print('ERROR: Change multiplication factor to a number >=1')
        return
    elif mf==1:
        print('This will prodice a four-way plot, not the scorecard')
        nota=info_dict2['note']

    use_notation=info_dict2['use_notation']
    print('The dataset includes ',the_df.shape[0],' entries in total')
    labels=[xc.upper() for xc in colori]
    incl_ave=info_dict2['incl aver']
      
    labels_ave=[xc.upper() for xc in other_colori]
    ###############################################################################################################################################
    ###############################################################################################################################################
    ###############################################################################################################################################
    fig, ax = plt.subplots(figsize=(info_dict2['fig_size'], info_dict2['fig_size']))
    texts1,texts2,texts3,texts4,texts5=[],[],[],[],[]
    p_val_x1,p_val_x2,p_val_x3,p_val_x4,p_val_x5=[],[],[],[],[]
    p_val_y1,p_val_y2,p_val_y3,p_val_y4,p_val_y5=[],[],[],[],[]

    texts100,texts200,texts300=[],[],[]
    p_val_x100,p_val_x200,p_val_x300=[],[],[]
    p_val_y100,p_val_y200,p_val_y300=[],[],[]
    
    info_dict2['Total entries']=the_df.shape[0]
    common_up= the_df[(the_df['FC cond x'] > 0) & (the_df['FC cond y'] > 0) ] # FIRST QUADRANT 
    print('First quadrant will host ',common_up.shape[0],' entries')
    info_dict2['Quadrant entries']=common_up.shape[0]
    lista=common_up[gene_name].tolist()    
    for my_gene in lista:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y']
        if mf > 1:
            if fch_x>th_fold_change*mf and fch_y>th_fold_change*mf:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)
            elif (fch_x>th_fold_change and fch_x<=th_fold_change*mf) and (fch_y>th_fold_change and fch_y<=th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts100.append( ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[0]  ))
                    p_val_x100.append(pval_x)
                    p_val_y100.append(pval_y)
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
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts200.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[1]  ))
                    p_val_x200.append(pval_x) #S
                    p_val_y200.append(pval_y)
            elif (fch_x<=th_fold_change ) and (fch_y>th_fold_change and fch_y<th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts300.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[2]  ))
                    p_val_x300.append(pval_x) #R
                    p_val_y300.append(pval_y)
            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        
        elif mf ==1:
            if fch_x>th_fold_change and fch_y>th_fold_change:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)
            elif (fch_x>th_fold_change) and (fch_y<=th_fold_change) :
                ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                    p_val_x4.append(pval_x)
                    p_val_y4.append(pval_y)
            elif (fch_x<=th_fold_change ) and (fch_y>th_fold_change):
                ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                    p_val_x5.append(pval_x)
                    p_val_y5.append(pval_y)
            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])            
    if mf >1:
        ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    if mf >1 and incl_ave==False:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    if mf >1 and incl_ave==True:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5,texts100,texts200,texts300]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    elif mf ==1:
        adjust_text(flatten([texts1,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    if IS_EXAMPLE==True and mf >1:
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
        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300
    elif IS_EXAMPLE==False and mf >1:
        labels_ave=['M','S','R']
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

        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300
    elif IS_EXAMPLE==True and mf ==1:
        quadr[labels[0].lower()]=[el.get_text() for el in texts1]

        quadr[labels[3].lower()]=[el.get_text() for el in texts4]
        quadr[labels[4].lower()]=[el.get_text() for el in texts5]
        quadr[labels[0].lower()+'_x']=[el.get_position()[0] for el in texts1]

        quadr[labels[3].lower()+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4].lower()+'_x']=[el.get_position()[0] for el in texts5]
        quadr[labels[0].lower()+'_y']=[el.get_position()[1] for el in texts1]

        quadr[labels[3].lower()+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4].lower()+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista
        quadr[labels[0].lower()+'_pval_x']=p_val_x1
        quadr[labels[0].lower()+'_pval_y']=p_val_y1

        quadr[labels[3].lower()+'_pval_x']=p_val_x4
        quadr[labels[3].lower()+'_pval_y']=p_val_y4
        quadr[labels[4].lower()+'_pval_x']=p_val_x5
        quadr[labels[4].lower()+'_pval_y']=p_val_y5
        
    elif IS_EXAMPLE==False and mf ==1:
        quadr['a']=[el.get_text() for el in texts1]

        quadr['d']=[el.get_text() for el in texts4]
        quadr['e']=[el.get_text() for el in texts5]
        quadr['a_x']=[el.get_position()[0] for el in texts1]

        quadr['d_x']=[el.get_position()[0] for el in texts4]
        quadr['e_x']=[el.get_position()[0] for el in texts5]
        quadr['a_y']=[el.get_position()[1] for el in texts1]

        quadr['d_y']=[el.get_position()[1] for el in texts4]
        quadr['e_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista
        quadr['a_pval_x']=p_val_x1
        quadr['a_pval_y']=p_val_y1

        quadr['d_pval_x']=p_val_x4
        quadr['d_pval_y']=p_val_y4
        quadr['e_pval_x']=p_val_x5
        quadr['e_pval_y']=p_val_y5
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
    if mf>1:
        ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change*mf), (xmax-th_fold_change*mf), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((th_fold_change, th_fold_change*mf), (th_fold_change*mf-th_fold_change), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change*mf, th_fold_change), (xmax-th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change*mf, ymin), (xmax-th_fold_change*mf), (th_fold_change-ymin),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmin, th_fold_change*mf), (th_fold_change-xmin), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    elif mf==1:
        ax.add_patch(Rectangle((th_fold_change, th_fold_change), (xmax-th_fold_change), (ymax-th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))        
        ax.add_patch(Rectangle((th_fold_change, ymin), (xmax-th_fold_change), (th_fold_change-ymin),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmin, th_fold_change), (th_fold_change-xmin), (ymax-th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
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
    texts100,texts200,texts300=[],[],[]
    p_val_x100,p_val_x200,p_val_x300=[],[],[]
    p_val_y100,p_val_y200,p_val_y300=[],[],[]
    common_down=the_df[(the_df['FC cond x'] < 0) & (the_df['FC cond y'] < 0) ]
    print('Third quadrant will host ',common_down.shape[0],' entries')
    info_dict2['Quadrant entries']=common_down.shape[0]
    lista2=common_down[gene_name].tolist() 
    for my_gene in lista2:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y']
        if mf>1:
            if fch_x<-th_fold_change*mf and fch_y<-th_fold_change*mf:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)
            elif (fch_x<-th_fold_change and fch_x>=-th_fold_change*mf) and (fch_y<-th_fold_change and fch_y>=-th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts100.append( ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[0]  ))        
                    p_val_x100.append(pval_x)
                    p_val_y100.append(pval_y)
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
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts200.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[1] ))
                    p_val_x200.append(pval_x)
                    p_val_y200.append(pval_y)
            elif (fch_x>=-th_fold_change ) and (fch_y<-th_fold_change and fch_y>-th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts300.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[2]  ))
                    p_val_x300.append(pval_x)
                    p_val_y300.append(pval_y)
            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        
        elif mf==1:
            if fch_x<-th_fold_change and fch_y<-th_fold_change:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)

            elif (fch_x<-th_fold_change) and (fch_y>=-th_fold_change) :
                ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                    p_val_x4.append(pval_x)
                    p_val_y4.append(pval_y)
            elif (fch_x>=-th_fold_change ) and (fch_y<-th_fold_change):
                ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                    p_val_x5.append(pval_x)
                    p_val_y5.append(pval_y)

            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])            
    if mf >1:
        ax.axhline(y=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    #adjust_text(texts1, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts2, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts3, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts4, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts5, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    if mf >1 and incl_ave==False:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    if mf >1 and incl_ave==True:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5,texts100,texts200,texts300]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    elif mf ==1:
        adjust_text(flatten([texts1,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    if IS_EXAMPLE==True and mf >1:
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
        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300        
    elif IS_EXAMPLE==False and mf >1:
        labels_ave=['M','S','R']
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
        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300
    elif IS_EXAMPLE==True and mf ==1:
        quadr[labels[0].lower()]=[el.get_text() for el in texts1]

        quadr[labels[3].lower()]=[el.get_text() for el in texts4]
        quadr[labels[4].lower()]=[el.get_text() for el in texts5]
        quadr[labels[0].lower()+'_x']=[el.get_position()[0] for el in texts1]
       
        quadr[labels[3].lower()+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4].lower()+'_x']=[el.get_position()[0] for el in texts5]
        quadr[labels[0].lower()+'_y']=[el.get_position()[1] for el in texts1]
        
        quadr[labels[3].lower()+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4].lower()+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista2
        quadr[labels[0].lower()+'_pval_x']=p_val_x1
        quadr[labels[0].lower()+'_pval_y']=p_val_y1
     
        quadr[labels[3].lower()+'_pval_x']=p_val_x4
        quadr[labels[3].lower()+'_pval_y']=p_val_y4
        quadr[labels[4].lower()+'_pval_x']=p_val_x5
        quadr[labels[4].lower()+'_pval_y']=p_val_y5
        
    elif IS_EXAMPLE==False and mf ==1:
        quadr['a']=[el.get_text() for el in texts1]
        
        quadr['d']=[el.get_text() for el in texts4]
        quadr['e']=[el.get_text() for el in texts5]
        quadr['a_x']=[el.get_position()[0] for el in texts1]
      
        quadr['d_x']=[el.get_position()[0] for el in texts4]
        quadr['e_x']=[el.get_position()[0] for el in texts5]
        quadr['a_y']=[el.get_position()[1] for el in texts1]
     
        quadr['d_y']=[el.get_position()[1] for el in texts4]
        quadr['e_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista2
        quadr['a_pval_x']=p_val_x1
        quadr['a_pval_y']=p_val_y1
      
        quadr['d_pval_x']=p_val_x4
        quadr['d_pval_y']=p_val_y4
        quadr['e_pval_x']=p_val_x5
        quadr['e_pval_y']=p_val_y5
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
    if mf>1:
        ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change*mf), (xmin+th_fold_change*mf), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change, -th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change*mf, -th_fold_change), (xmin+th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change*mf, ymax), (xmin+th_fold_change*mf), -(ymax+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmax, -th_fold_change*mf), -(xmax+th_fold_change), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    elif mf==1:
        ax.add_patch(Rectangle((-th_fold_change, -th_fold_change), (xmin+th_fold_change), (ymin+th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change, ymax), (xmin+th_fold_change), -(ymax+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmax, -th_fold_change), -(xmax+th_fold_change), (ymin+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))        

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
    texts100,texts200,texts300=[],[],[]
    p_val_x100,p_val_x200,p_val_x300=[],[],[]
    p_val_y100,p_val_y200,p_val_y300=[],[],[]
    common_1=the_df[(the_df['FC cond x'] < 0) & (the_df['FC cond y'] > 0) ] # Quadrant II
    print('Second quadrant will host ',common_1.shape[0],' entries')
    info_dict2['Quadrant entries']=common_1.shape[0]
    lista3=common_1[gene_name].tolist()
    for my_gene in lista3:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y']
        if mf >1:
            if fch_x<-th_fold_change*mf and fch_y>th_fold_change*mf:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)
            elif (fch_x<-th_fold_change and fch_x>=-th_fold_change*mf) and (fch_y>th_fold_change and fch_y<=th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts100.append( ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[0]  ))        
                    p_val_x100.append(pval_x)
                    p_val_y100.append(pval_y)
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
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts200.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[1] ))
                    p_val_x200.append(pval_x)
                    p_val_y200.append(pval_y)
            elif (fch_x>=-th_fold_change ) and (fch_y>th_fold_change and fch_y<th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts300.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[2]  ))
                    p_val_x300.append(pval_x)
                    p_val_y300.append(pval_y)
            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        
        elif mf==1:
            if fch_x<-th_fold_change and fch_y>th_fold_change:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)
            elif (fch_x<-th_fold_change) and (fch_y<=th_fold_change) :
                ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                    p_val_x4.append(pval_x)
                    p_val_y4.append(pval_y)
            elif (fch_x>=-th_fold_change ) and (fch_y>th_fold_change):
                ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                    p_val_x5.append(pval_x)
                    p_val_y5.append(pval_y)
            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])              
    if mf>1:
        ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    #adjust_text(texts1, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts2, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts3, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts4, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts5, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    if mf >1 and incl_ave==False:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    if mf >1 and incl_ave==True:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5,texts100,texts200,texts300]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    elif mf ==1:
        adjust_text(flatten([texts1,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    
    if IS_EXAMPLE==True and mf >1:
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
        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300       
    elif IS_EXAMPLE==False and mf >1:
        labels_ave=['M','S','R']
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
        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300
    elif IS_EXAMPLE==True and mf ==1:
        quadr[labels[0].lower()]=[el.get_text() for el in texts1]

        quadr[labels[3].lower()]=[el.get_text() for el in texts4]
        quadr[labels[4].lower()]=[el.get_text() for el in texts5]
        quadr[labels[0].lower()+'_x']=[el.get_position()[0] for el in texts1]
 
        quadr[labels[3].lower()+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4].lower()+'_x']=[el.get_position()[0] for el in texts5]
        quadr[labels[0].lower()+'_y']=[el.get_position()[1] for el in texts1]
      
        quadr[labels[3].lower()+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4].lower()+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista3
        quadr[labels[0].lower()+'_pval_x']=p_val_x1
        quadr[labels[0].lower()+'_pval_y']=p_val_y1
  
        quadr[labels[3].lower()+'_pval_x']=p_val_x4
        quadr[labels[3].lower()+'_pval_y']=p_val_y4
        quadr[labels[4].lower()+'_pval_x']=p_val_x5
        quadr[labels[4].lower()+'_pval_y']=p_val_y5
        
    elif IS_EXAMPLE==False and mf ==1:
        quadr['a']=[el.get_text() for el in texts1]

        quadr['d']=[el.get_text() for el in texts4]
        quadr['e']=[el.get_text() for el in texts5]
        quadr['a_x']=[el.get_position()[0] for el in texts1]
      
        quadr['d_x']=[el.get_position()[0] for el in texts4]
        quadr['e_x']=[el.get_position()[0] for el in texts5]
        quadr['a_y']=[el.get_position()[1] for el in texts1]
       
        quadr['d_y']=[el.get_position()[1] for el in texts4]
        quadr['e_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista3
        quadr['a_pval_x']=p_val_x1
        quadr['a_pval_y']=p_val_y1
       
        quadr['d_pval_x']=p_val_x4
        quadr['d_pval_y']=p_val_y4
        quadr['e_pval_x']=p_val_x5
        quadr['e_pval_y']=p_val_y5


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
    if mf>1:
        ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change*mf), (xmin+th_fold_change*mf), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change, th_fold_change*mf), (-th_fold_change*mf+th_fold_change), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change*mf, th_fold_change), (xmin+th_fold_change*mf), (th_fold_change*mf-th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((-th_fold_change*mf, ymin), (xmin+th_fold_change*mf), (th_fold_change-ymin),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmax, th_fold_change*mf), -(xmax+th_fold_change), (ymax-th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    elif mf==1:
        ax.add_patch(Rectangle((-th_fold_change, th_fold_change), (xmin+th_fold_change), (ymax-th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((-th_fold_change, ymin), (xmin+th_fold_change), (th_fold_change-ymin),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmax, th_fold_change), -(xmax+th_fold_change), (ymax-th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))        
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
    texts100,texts200,texts300=[],[],[]
    p_val_x100,p_val_x200,p_val_x300=[],[],[]
    p_val_y100,p_val_y200,p_val_y300=[],[],[]
    common_2=the_df[(the_df['FC cond x'] > 0) & (the_df['FC cond y'] < 0) ] # Quadrant IV
    print('Forth quadrant will host ',common_2.shape[0],' entries')
    info_dict2['Quadrant entries']=common_2.shape[0]
    lista4=common_2[gene_name].tolist()
    for my_gene in lista4:
        fch_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond x']
        fch_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['FC cond y']
        pval_x=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond x']
        pval_y=  the_df.loc[the_df[gene_name] == my_gene].iloc[0]['padj cond y']
        if mf>1:
            if fch_x>th_fold_change*mf and fch_y<-th_fold_change*mf:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)
            elif (fch_x>th_fold_change and fch_x<=th_fold_change*mf) and (fch_y<-th_fold_change and fch_y>=-th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[1],marker=MarkerStyle(markers[1], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts100.append( ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[0]  ))        
                    p_val_x100.append(pval_x)
                    p_val_y100.append(pval_y)
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
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts200.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[1] ))
                    p_val_x200.append(pval_x)
                    p_val_y200.append(pval_y)
            elif (fch_x<=th_fold_change ) and (fch_y<-th_fold_change and fch_y>-th_fold_change*mf):
                ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[2],marker=MarkerStyle(markers[4], fillstyle='full'),s=sizes[2])
                if pval_x<th_significance and pval_y<th_significance and incl_ave==True:
                    texts300.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[2]  ))
                    p_val_x300.append(pval_x)
                    p_val_y300.append(pval_y)
            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])                        
        elif mf==1:
            if fch_x>th_fold_change and fch_y<-th_fold_change:
                ax.scatter( fch_x,fch_y, facecolors = colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[0], fillstyle='full'),s=sizes[0])        
                if pval_x<th_significance and pval_y<th_significance:
                    texts1.append( ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[0]  ))        
                    p_val_x1.append(pval_x)
                    p_val_y1.append(pval_y)
            elif (fch_x>th_fold_change) and (fch_y>=-th_fold_change) :
                ax.scatter( fch_x,fch_y, facecolors = colori[3], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts4.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[3]  ))
                    p_val_x4.append(pval_x)
                    p_val_y4.append(pval_y)
            elif (fch_x<=th_fold_change ) and (fch_y<-th_fold_change):
                ax.scatter( fch_x,fch_y, facecolors = colori[4], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[1])
                if pval_x<th_significance and pval_y<th_significance:
                    texts5.append(ax.text(fch_x,fch_y, my_gene,size=font_size1, ha='center', va='center',color=colori[4]  ))
                    p_val_x5.append(pval_x)
                    p_val_y5.append(pval_y)
            else:
                ax.scatter( fch_x,fch_y, facecolors = other_colori[-1], edgecolors = "k", linewidths = 0.1, alpha = trasp[-1],marker=MarkerStyle(markers[-1], fillstyle='full'),s=sizes[-1])              
    if mf>1:
        ax.axhline(y=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
    ax.axhline(y=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
    ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

    #adjust_text(texts1, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts2, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts3, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts4, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    #adjust_text(texts5, ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    if mf >1 and incl_ave==False:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    if mf >1 and incl_ave==True:
        adjust_text(flatten([texts1,texts2,texts3,texts4,texts5,texts100,texts200,texts300]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    elif mf ==1:
        adjust_text(flatten([texts1,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    quadr={}
    if IS_EXAMPLE==True and mf>1:
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
        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300        
    elif IS_EXAMPLE==False and mf>1:
        labels_ave=['M','S','R']
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
        if incl_ave==True:
            quadr[labels_ave[0]]=[el.get_text() for el in texts100]
            quadr[labels_ave[1]]=[el.get_text() for el in texts200]
            quadr[labels_ave[2]]=[el.get_text() for el in texts300]
            quadr[labels_ave[0]+'_x']=[el.get_position()[0] for el in texts100]
            quadr[labels_ave[1]+'_x']=[el.get_position()[0] for el in texts200]
            quadr[labels_ave[2]+'_x']=[el.get_position()[0] for el in texts300]
            quadr[labels_ave[0]+'_y']=[el.get_position()[1] for el in texts100]
            quadr[labels_ave[1]+'_y']=[el.get_position()[1] for el in texts200]
            quadr[labels_ave[2]+'_y']=[el.get_position()[1] for el in texts300]
            quadr[labels_ave[0]+'_pval_x']=p_val_x100
            quadr[labels_ave[0]+'_pval_y']=p_val_y100
            quadr[labels_ave[1]+'_pval_x']=p_val_x200
            quadr[labels_ave[1]+'_pval_y']=p_val_y200
            quadr[labels_ave[2]+'_pval_x']=p_val_x300
            quadr[labels_ave[2]+'_pval_y']=p_val_y300
    elif IS_EXAMPLE==True and mf==1:
        quadr[labels[0].lower()]=[el.get_text() for el in texts1]    
        quadr[labels[3].lower()]=[el.get_text() for el in texts4]
        quadr[labels[4].lower()]=[el.get_text() for el in texts5]
        
        quadr[labels[0].lower()+'_x']=[el.get_position()[0] for el in texts1]       
        quadr[labels[3].lower()+'_x']=[el.get_position()[0] for el in texts4]
        quadr[labels[4].lower()+'_x']=[el.get_position()[0] for el in texts5]
        
        quadr[labels[0].lower()+'_y']=[el.get_position()[1] for el in texts1]      
        quadr[labels[3].lower()+'_y']=[el.get_position()[1] for el in texts4]
        quadr[labels[4].lower()+'_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista4
        quadr[labels[0].lower()+'_pval_x']=p_val_x1
        quadr[labels[0].lower()+'_pval_y']=p_val_y1
     
        quadr[labels[3].lower()+'_pval_x']=p_val_x4
        quadr[labels[3].lower()+'_pval_y']=p_val_y4
        quadr[labels[4].lower()+'_pval_x']=p_val_x5
        quadr[labels[4].lower()+'_pval_y']=p_val_y5
        
    elif IS_EXAMPLE==False and mf==1:
        quadr['a']=[el.get_text() for el in texts1]
   
        quadr['d']=[el.get_text() for el in texts4]
        quadr['e']=[el.get_text() for el in texts5]
        quadr['a_x']=[el.get_position()[0] for el in texts1]
       
        quadr['d_x']=[el.get_position()[0] for el in texts4]
        quadr['e_x']=[el.get_position()[0] for el in texts5]
        quadr['a_y']=[el.get_position()[1] for el in texts1]
      
        quadr['d_y']=[el.get_position()[1] for el in texts4]
        quadr['e_y']=[el.get_position()[1] for el in texts5]
        #quadr['ALL ENTRIES']=lista4
        quadr['a_pval_x']=p_val_x1
        quadr['a_pval_y']=p_val_y1
      
        quadr['d_pval_x']=p_val_x4
        quadr['d_pval_y']=p_val_y4
        quadr['e_pval_x']=p_val_x5
        quadr['e_pval_y']=p_val_y5
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
    if mf>1:
        ax.add_patch(Rectangle((th_fold_change*mf, -th_fold_change*mf), (xmax-th_fold_change*mf), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((th_fold_change, -th_fold_change*mf), (th_fold_change*mf-th_fold_change), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change*mf,-th_fold_change), (xmax-th_fold_change*mf), (-th_fold_change*mf+th_fold_change),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
        ax.add_patch(Rectangle((th_fold_change*mf, ymax), (xmax-th_fold_change*mf), -(ymax+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmin, -th_fold_change*mf), (th_fold_change-xmin), (ymin+th_fold_change*mf),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
    elif mf==1:
        ax.add_patch(Rectangle((th_fold_change, -th_fold_change), (xmax-th_fold_change), (ymin+th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
        ax.add_patch(Rectangle((th_fold_change, ymax), (xmax-th_fold_change), -(ymax+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        ax.add_patch(Rectangle((xmin, -th_fold_change), (th_fold_change-xmin), (ymin+th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))        
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
        other_colori=my_data[quadrante]['params']['other_colors']
        mf=my_data[quadrante]['params']['multiplication factor']
        if IS_EXAMPLE:
            if mf>1:
                etichette=[xc.upper() for xc in colori]
                etichette2=[xc.upper() for xc in other_colori]
            elif mf==1:
                etichette=[xc.lower() for xc in colori]
                etichette2=[xc.lower() for xc in other_colori]                
        else:
            etichette=['A','B','C','D','E']
            etichette2=['M','S','R']
            if mf==1:
                etichette=[i_v_s.lower() for i_v_s in etichette]
                etichette2=[i_v_s.lower() for i_v_s in etichette2]
        incl_ave=my_data[quadrante]['params']['incl aver']
        if incl_ave:
            etichette=etichette+etichette2
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
        
        if mf==1:
            print('Reconstruction of a Four-Way plot. Scorecard was not created!')
        incl_ave=my_data[quadrante]['params']['incl aver']    
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        minimo,massimo=-fig_size,fig_size    
        ax.set_xlim(left=minimo,right=massimo)
        ax.set_ylim(bottom=minimo,top=massimo)
        if mf>1:
            ax.axhline(y=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
        ax.axhline(y=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
        ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
        if mf>1:
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
        texts100,texts200,texts300=[],[],[]
        if mf>1:
            
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
        elif mf==1:

            ax.add_patch(Rectangle((th_fold_change, th_fold_change), (massimo-th_fold_change), (massimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
            ax.add_patch(Rectangle((-th_fold_change, -th_fold_change), (-massimo+th_fold_change), (-massimo+th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
            ax.add_patch(Rectangle((-th_fold_change, th_fold_change), (-massimo+th_fold_change), (massimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
            ax.add_patch(Rectangle((th_fold_change, -th_fold_change), (massimo-th_fold_change), (-massimo+th_fold_change),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))

            ax.add_patch(Rectangle((th_fold_change, -th_fold_change), (massimo-th_fold_change), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
            ax.add_patch(Rectangle((-th_fold_change, th_fold_change), (th_fold_change*2), (massimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
            ax.add_patch(Rectangle((minimo,-th_fold_change), (-minimo-th_fold_change), (th_fold_change*2),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
            ax.add_patch(Rectangle((-th_fold_change,minimo), (th_fold_change*2), (-minimo-th_fold_change),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
        for quadrante in quadr_list:
            for eti in etichette:
                if eti in list(my_data[quadrante].keys()):
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
                        elif incl_ave==True and eti==etichette[5] :
                            ax.scatter( fch_x,fch_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[2])
                            texts100.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[0]  ))
                        elif incl_ave==True and eti==etichette[6]:
                            ax.scatter( fch_x,fch_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[2])
                            texts200.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[1]  ))
                        elif incl_ave==True and eti==etichette[7]:
                            ax.scatter( fch_x,fch_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[2])
                            texts300.append(ax.text(fch_x,fch_y, my_gene,size=font_size1-1, ha='center', va='center',color=other_colori[2]  ))                
        if incl_ave==False:
            adjust_text(flatten([texts1,texts2,texts3,texts4,texts5]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        else:
            adjust_text(flatten([texts1,texts2,texts3,texts4,texts5,texts100,texts200,texts300]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        if use_notation:   
            ax.set_xlabel("$log_2$ Fold Change ("+trt1+" vs "+ctrl+")")
            ax.set_ylabel("$log_2$ Fold Change ("+trt2+" vs "+ctrl+")")
            ax.set_title('Scorecard ('+titolo+')')
        else:
            ax.set_xlabel("$log_2$ Fold Change ("+trt1+")")
            ax.set_ylabel("$log_2$ Fold Change ("+trt2+")")   
            ax.set_title('Scorecard and regions of interest')
        if mf>1:
            plt.savefig(save_folder+'Scorecard.png',dpi=300,bbox_inches='tight')
        elif mf==1:
            plt.savefig(save_folder+'FourWayPlot.png',dpi=300,bbox_inches='tight')
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

    Only works for the scorecard (two thresholds) not for the four-way plot.
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
            if tmp2['params']['multiplication factor']==1:
                print('Sorry only working for the scorecard; run using multiplication factor >1')
                return
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
        mf=my_data[quadrante]['params']['multiplication factor']
        if IS_EXAMPLE:
            if mf>1:
                etichette=[xc.upper() for xc in colori]
                etichette2=[xc.upper() for xc in other_colori]
            elif mf==1:
                etichette=[xc.lower() for xc in colori]
                etichette2=[xc.lower() for xc in other_colori]                
        else:
            etichette=['A','B','C','D','E']
            etichette2=['M','S','R']
            if mf==1:
                etichette=[i_v_s.lower() for i_v_s in etichette]
                etichette2=[i_v_s.lower() for i_v_s in etichette2]
        incl_ave=my_data[quadrante]['params']['incl aver']
        if incl_ave:
            etichette=etichette+etichette2
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
        use_notation=my_data[quadrante]['params']['use_notation']            
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(18, 8),sharey=True)
        texts1x,texts2x,texts3x,texts4x,texts5x=[],[],[],[],[]
        texts1y,texts2y,texts3y,texts4y,texts5y=[],[],[],[],[]

        texts100x,texts200x,texts300x=[],[],[]
        texts100y,texts200y,texts300y=[],[],[]
        
        for quadrante in quadr_list:
            for eti in etichette:
                if eti in list(my_data[quadrante].keys()):
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

                            elif incl_ave==True and eti==etichette[5]:
                                ax1.scatter( fch_x,pval_x, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[2])
                                texts3x.append(ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=other_colori[0]  ))
                                ax2.scatter( fch_y,pval_y, facecolors = other_colori[0], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[2], fillstyle='full'),s=sizes[2])
                                texts3y.append(ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=other_colori[0]  ))
                            elif incl_ave==True and eti==etichette[6]:
                                ax1.scatter( fch_x,pval_x, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[2])
                                texts4x.append(ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=other_colori[1]  ))
                                ax2.scatter( fch_y,pval_y, facecolors = other_colori[1], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[2])
                                texts4y.append(ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=other_colori[1]  ))
                            elif incl_ave==True and eti==etichette[7]:
                                ax1.scatter( fch_x,pval_x, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[2])
                                texts5x.append(ax1.text(fch_x,pval_x, my_gene,size=font_size1, ha='center', va='center',color=other_colori[2]  ))
                                ax2.scatter( fch_y,pval_y, facecolors = other_colori[2], edgecolors = "k", linewidths = 0.1, alpha = trasp[0],marker=MarkerStyle(markers[3], fillstyle='full'),s=sizes[2])
                                texts5y.append(ax2.text(fch_y,pval_y, my_gene,size=font_size1, ha='center', va='center',color=other_colori[2]  ))
        if incl_ave==False :
            out_x=flatten([texts1x,texts2x,texts3x,texts4x,texts5x])
            out_y=flatten([texts1y,texts2y,texts3y,texts4y,texts5y])
        else:
            out_x=flatten([texts1x,texts2x,texts3x,texts4x,texts5x,texts100x,texts200x,texts300x])
            out_y=flatten([texts1y,texts2y,texts3y,texts4y,texts5y,texts100y,texts200y,texts300y])            
 

        for ax,tmp_str in zip([ax1,ax2],[trt1,trt2]):    
            if mf>1:
                ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
            ax.axvline(x=0,color='grey',linestyle='dotted',lw=0.5)
            if mf>1:
                ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)

            ax.axhline(y=-np.log10(th_significance),color='grey',linestyle='dashdot',lw=1.5)
            #ax.axhline(y=-np.log10(th_significance/mf),color='grey',linestyle='dashdot',lw=1.0)

            ax.set_xlabel("log2 Fold Change ("+tmp_str+" vs "+ctrl+")")
            ax.set_ylabel("-log10 Adjusted p-value")
            ax.set_title('Data: '+tmp_str+' versus '+ctrl)

            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            if mf>1:
                ax.add_patch(Rectangle((xmin, -np.log10(th_significance)), -(xmin+th_fold_change*mf), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
                ax.add_patch(Rectangle((th_fold_change*mf, -np.log10(th_significance)), (xmax-th_fold_change*mf), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2])                         )
                ax.add_patch(Rectangle((-th_fold_change*mf, -np.log10(th_significance)), -(-th_fold_change*mf+th_fold_change), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
                ax.add_patch(Rectangle((th_fold_change, -np.log10(th_significance)), (th_fold_change*mf-th_fold_change), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[1],alpha=trasp_rect[1]))
                ax.add_patch(Rectangle((-th_fold_change, -np.log10(th_significance)), (th_fold_change*2), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[0],alpha=trasp_rect[0]))
            elif mf==1:
                ax.add_patch(Rectangle((xmin, -np.log10(th_significance)), -(xmin+th_fold_change), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2]))
                ax.add_patch(Rectangle((th_fold_change, -np.log10(th_significance)), (xmax-th_fold_change), (ymax-np.log10(th_significance)),edgecolor='none' ,facecolor =col_rect[2],alpha=trasp_rect[2])                         )
        if len(out_x)>0 :
            adjust_text(out_x, ax=ax1,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        if len(out_y)>0:
            adjust_text(out_y, ax=ax2,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
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
        incl_ave=my_data[quadrante]['params']['incl aver']
        mf=my_data[quadrante]['params']['multiplication factor']
        if IS_EXAMPLE:
            if mf>1:
                etichette=[xc.upper() for xc in colori]
                etichette2=[xc.upper() for xc in other_colori]
            elif mf==1:
                etichette=[xc.lower() for xc in colori]
                etichette2=[xc.lower() for xc in other_colori]                
        else:
            etichette=['A','B','C','D','E']
            etichette2=['M','S','R']
            if mf==1:
                etichette=[i_v_s.lower() for i_v_s in etichette]
                etichette2=[i_v_s.lower() for i_v_s in etichette2]
        if incl_ave:
            etichette=etichette+etichette2
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
        use_notation=my_data[quadrante]['params']['use_notation']        
        Position=0
        fig = plt.figure()
        ax = fig.add_subplot(111)
        all_x,all_y,all_genes,all_colors,labels_x,labels_y=[],[],[],[],[],[]
        for quadrante in quadr_list:
            for eti in etichette:
                if eti in list(my_data[quadrante].keys()):
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

                            elif incl_ave==True and eti==etichette[5]:
                                all_x.append( fch_x)
                                all_y.append( fch_y)
                                all_genes.append( my_gene)
                                all_colors.append(other_colori[0])
                                Position = Position+1
                                labels_x.append(trt1)
                                labels_y.append(trt2)
                            elif incl_ave==True and eti==etichette[6]:
                                all_x.append( fch_x)
                                all_y.append( fch_y)
                                all_genes.append( my_gene)
                                all_colors.append(other_colori[1])
                                Position = Position+1
                                labels_x.append(trt1)
                                labels_y.append(trt2)
                            elif incl_ave==True and eti==etichette[7]:
                                all_x.append( fch_x)
                                all_y.append( fch_y)
                                all_genes.append( my_gene)
                                all_colors.append(other_colori[2])
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
            if mf>1:
                ax.axvline(x=th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
            ax.axvline(x=0,color='grey',linestyle='dotted',lw=0.5)
            if mf>1:
                ax.axvline(x=-th_fold_change*mf,color='grey',linestyle='dashdot',lw=1.5)
            ax.axvline(x=-th_fold_change,color='grey',linestyle='dashdot',lw=1.0)
            ax.tick_params(axis='y', labelsize=5)
            if the_folder[-1]!="/":
                the_folder=the_folder+"/"
            ax.spines[['right', 'top']].set_visible(False)
            if try_adj_test:
                adjust_text(flatten([text1,text2]), ax=ax,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
            plt.savefig(the_folder+'Bars.png',dpi=300,bbox_inches='tight')
            plt.close()
def save_to_file(thedir,*text):    
    with open(thedir+'log_info.txt', mode='wt', encoding='utf-8') as myfile:
        for lines in text:
            myfile.write('\n'.join(str(line) for line in lines))
            myfile.write('\n')
def count_frequencies(my_directory):
    '''
    Summary of the entries/genes identified through the scorecard.
    Reports all comparisons perfomed.
    Input the main_folder as string.
    '''
    if my_directory[-1]!="/":
        my_directory=my_directory+"/"
    my_log=[]   
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
    str1='Number of comparisons performed :'+str(len(all_dir))
    print(str1)
    print('\n')
    my_log.append(str1)
    for i,k in enumerate(VARIABLES):
        str0='-------------- Comparison '+str(i+1)+' of '+str(len(all_dir))
        str01='Comparison : '+str(k)
        print(str0)
        print(str01)
        my_log.append(str0)
        my_log.append(str01)
        conti=0
        conti2=0
        initi=True
        list_q=[]
        all_counts=[]
        list_q2=[]
        all_counts2=[]
        for qi,qr in enumerate(quadr_list):
            tmp=all_data[k][qr]
            n_entries=tmp['params']['Quadrant entries']
            tot_entries=tmp['params']['Total entries']
            ctrl=tmp['params']['Control name']
            IS_EXAMPLE=tmp['params']['is_example']
            colori=tmp['params']['colors']
            other_colori=tmp['params']['other_colors']
            incl_ave=tmp['params']['incl aver']
            mf=tmp['params']['multiplication factor']
            if IS_EXAMPLE:
                if mf>1:
                    etichette=[xc.upper() for xc in colori]
                    etichette2=[xc.upper() for xc in other_colori]
                elif mf==1:
                    etichette=[xc.lower() for xc in colori]
                    etichette2=[xc.lower() for xc in other_colori]                
            else:
                etichette=['A','B','C','D','E']
                etichette2=['M','S','R']
                if mf==1:
                    etichette=[i_v_s.lower() for i_v_s in etichette]
                    etichette2=[i_v_s.lower() for i_v_s in etichette2]

            str8=qr+' total entries were '+str(n_entries)
            print(str8)
            titolo=tmp['params']['Scorecard title']
            fig_size=tmp['params']['fig_size']
            save_folder=tmp['params']['save_dir']
            trt1=tmp['params']['Treatment1 name']
            trt2=tmp['params']['Treatment2 name']
            descr={'Quadrant1':'both '+trt1+' and '+trt2+' vs '+ctrl+' up-reg.',
                   'Quadrant2':trt1+' down.reg. and '+trt2+' up-reg. (vs '+ctrl+')',
                   'Quadrant3':'both '+trt1+' and '+trt2+' vs '+ctrl+' down-reg.',
                   'Quadrant4':trt1+' up.reg. and '+trt2+' down-reg. (vs '+ctrl+')'}
            if save_folder[-1]!="/":
                save_folder=save_folder+"/"+trt1+" "+trt2+"/"
            else:
                save_folder=save_folder+trt1+" "+trt2+"/"
            if not isdir(save_folder):
                print('Run scorecard quadrants creation before attempting reconstruction')
            th_fold_change=tmp['params']['th_fold_change']
            th_significance=tmp['params']['th_significance']
            font_size1=tmp['params']['font_size_quadrants']
            trasp=tmp['params']['marker_trasp']
            trasp_rect=tmp['params']['rect_trasp']
            col_rect=tmp['params']['rect_colors']
            markers=tmp['params']['markers']
            sizes= tmp['params']['markers_sizes']
            gene_name= tmp['params']['gene_name']            
            use_notation=tmp['params']['use_notation']
            if initi:
                if mf>1:
                    strf='Applying fold change thresholds of '+str(th_fold_change)+' and '+str(th_fold_change*mf)
                elif mf==1:
                    strf='Applying fold change threshold of '+str(th_fold_change)
                print(strf)
                print('Also, filtering the data using a statistical threshold of ',th_significance)
                print('Analyzed ',tot_entries,' entries'+'\n')
                initi=False
                my_log.append(strf )
                my_log.append('Also, filtering the data using a statistical threshold of '+str(th_significance))
                my_log.append('Analyzed '+str(tot_entries)+' entries'+'\n' )           
            for key,value in tmp.items():
                if len(value)>0 and key in etichette:
                    str1='In '+qr+' the <<'+key+'>> group contains '+str(len(value))+' entries'
                    print(str1)                    
                    conti += len(value)
                    list_q.append(qr)
                    all_counts.append(value)
                    my_log.append(str1)
                elif incl_ave==True and len(value)>0 and key in etichette2:
                    str1='In '+qr+' the '+key+' group contains '+str(len(value))+' entries'
                    print(str1)                    
                    conti2 += len(value)
                    list_q2.append(qr)
                    all_counts2.append(value)
                    my_log.append(str1)
            my_log.append(str8+'\n')
        if incl_ave==False:
            for sel_q in set(list_q):
                str4=sel_q+' means '+descr[sel_q]            
                print(str4)
                #my_log.append(str4+'\n')
        else:
            for sel_q in set(list_q+list_q2):
                str4=sel_q+' means '+descr[sel_q]            
                print(str4)
                #my_log.append(str4+'\n')            
        str1='In total the scorecard identified '+str(conti)+' entries inside the shaded regions of interest'
        str11='The scorecard also identified '+str(conti2)+' entries between F.C. thresholds'
        print(str1)
        if incl_ave:
            print(str11)
        print('\n')
        
        my_log.append(str1+'\n')
        if incl_ave:
            my_log.append(str11+'\n')
    if incl_ave==False:
        calc=Counter(flatten(all_counts))
    else:
        calc=Counter(flatten(all_counts+all_counts2))
    multi_entr=False
    for key , value in calc.items():
        if value>1:
            str1='Entry :',key,' occurred ', value,' times during the ',len(all_dir),' comparisons'
            print(str1)
            my_log.append(str1)
            multi_entr=True
            for i,k in enumerate(VARIABLES):
                for qi,qr in enumerate(quadr_list):
                    tmp=all_data[k][qr]
                    if key in tmp.keys():
                        str1='Entry found in ',k,' [',qr,']'
                        print(str1)
                        my_log.append(str1)
    if multi_entr==False:
        print('No repeated entries found among experimental conditions being compared')
        my_log.append('No repeated entries found among experimental conditions being compared')
    save_to_file(my_directory,my_log)
def heatmap(data, row_labels, col_labels, ax=None,cbar_kw=None, cbarlabel="", **kwargs):

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels,fontsize=6)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels,fontsize=5)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,           labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-60, ha="right",   rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}", textcolors=("black", "white"), threshold=None,large_im=False, **textkw):

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
    if large_im:
        fs=3
    else:
        fs=4.5
    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None),fontsize=fs, **kw)
            texts.append(text)

    return texts
def quadrants_heatmap(my_directory):
    '''
    Overview of the number of identified entries for each group/quandrant.
    It summarizes the information as a annotated heatmap, but report only the number of entries (not the expression values)
    Pass the main_folder as argument.
    '''
    if my_directory[-1]!="/":
        my_directory=my_directory+"/"
    my_log=[]   
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
    out_all=[]
    for i,k in enumerate(VARIABLES):

        conti=0
        out={}
        list_q=[]
        all_counts=[]
        for qi,qr in enumerate(quadr_list):
            tmp=all_data[k][qr]
            n_entries=tmp['params']['Quadrant entries']
            tot_entries=tmp['params']['Total entries']
            ctrl=tmp['params']['Control name']
            IS_EXAMPLE=tmp['params']['is_example']
            colori=tmp['params']['colors']
            other_colori=tmp['params']['other_colors']
            th_fold_change=tmp['params']['th_fold_change']
            th_significance=tmp['params']['th_significance']
            font_size1=tmp['params']['font_size_quadrants']
            trasp=tmp['params']['marker_trasp']
            trasp_rect=tmp['params']['rect_trasp']
            col_rect=tmp['params']['rect_colors']
            markers=tmp['params']['markers']
            sizes= tmp['params']['markers_sizes']
            gene_name= tmp['params']['gene_name']
            mf=tmp['params']['multiplication factor']
            use_notation=tmp['params']['use_notation']
            incl_ave=tmp['params']['incl aver']
            if IS_EXAMPLE==True and mf>1:
                etichette=[xc.upper() for xc in colori]
                etichette2=[xc.upper() for xc in other_colori]
            elif IS_EXAMPLE==False and mf>1:
                etichette=['A','B','C','D','E']
                etichette2=['M','S','R']
            elif IS_EXAMPLE==True and mf==1:
                property_col = [colori[i] for i in [0,3,4]]
                etichette=[xc.lower() for xc in property_col]
            elif IS_EXAMPLE==False and mf==1:
                etichette=['a','d','e']
            titolo=tmp['params']['Scorecard title']
            fig_size=tmp['params']['fig_size']
            save_folder=tmp['params']['save_dir']
            trt1=tmp['params']['Treatment1 name']
            trt2=tmp['params']['Treatment2 name']

            if save_folder[-1]!="/":
                save_folder=save_folder+"/"+trt1+" "+trt2+"/"
            else:
                save_folder=save_folder+trt1+" "+trt2+"/"
            if not isdir(save_folder):
                print('Run scorecard quadrants creation before attempting reconstruction')

            qr_str=qr.replace('uadrant', '')
            for key,value in tmp.items():
                if len(value)>0 and key in etichette:
                    str1='In '+qr+' the '+key+' group contains '+str(len(value))+' entries'    
                    conti += len(value) # Tot number of entries inside all regions of interest
                    list_q.append(qr)
                    all_counts.append(value)
                    out[qr_str+' '+key]=int(len(value))
                elif len(value)==0 and key in etichette:
                    out[qr_str+' '+key]=0
                elif incl_ave==True and len(value)>0 and key in etichette2:
                    str1='In '+qr+' the '+key+' group contains '+str(len(value))+' entries'    
                    conti += len(value) # Tot number of entries inside all regions of interest
                    list_q.append(qr)
                    all_counts.append(value)
                    out[qr_str+' '+key]=int(len(value))
                elif incl_ave==True and len(value)==0 and key in etichette2:
                    out[qr_str+' '+key]=0
        out_all.append(out)                    
    df = pd.DataFrame.from_dict(out_all)
    fig, ax = plt.subplots()
    
    im, cbar = heatmap(df.to_numpy(dtype=np.int32).T,df.columns.tolist(), VARIABLES,  ax=ax, cmap="YlGn", cbarlabel="Counts")
    if len(VARIABLES)>6:
        texts = annotate_heatmap(im, valfmt="{x:.0f}",large_im=True)
    else:
        texts = annotate_heatmap(im, valfmt="{x:.0f}",large_im=False)

    fig.tight_layout()

    plt.savefig(my_directory+'HeatmapCounts.png',dpi=300,bbox_inches='tight')
    plt.close()
