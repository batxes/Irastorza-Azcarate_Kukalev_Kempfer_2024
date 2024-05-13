#%%

"""
Script that plots unphased, cast and s129 normalized compartments for the region given as input.
This plot was used in Figure 1 and SI Figure 2.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

def normalize (score,max_value,min_value):
    if score > 0:
        return score/max_value
    elif score < 0:
        return score/abs(min_value)
    else:
        return 0

def plot_compartments(chr_,start,end,text):
    cast_df = pd.read_csv("{}compartments/GAM_F123_ALL_as_3NP_CAST.{}.AB.txt".format(root,chr_),sep="\t",usecols=["chrom","start","end","score0"])    
    s129_df = pd.read_csv("{}compartments/GAM_F123_ALL_as_3NP_S129.{}.AB.txt".format(root,chr_),sep="\t",usecols=["chrom","start","end","score0"])    
    nonphased_df = pd.read_csv("{}compartments/GAM_F123_ALL_as_3NP_no_cleaning.{}.AB.txt".format(root,chr_),sep="\t",usecols=["chrom","start","end","score0"])    
    fig,axs = plt.subplots(3,1,figsize=(12,4))
    fig.tight_layout()
    max_ = 0
    min_ = 0
    for i,cast_df in enumerate([nonphased_df,cast_df,s129_df]):
        ax1 = axs[i]
        #max_value = max(abs(cast_df.score0.values))
        max_value = max(cast_df.score0.values)
        min_value = min(cast_df.score0.values)
        cast_df = cast_df.query('start >= {} and end <= {}'.format(start,end))
        new_scores = cast_df.apply(lambda x: normalize(x.score0,max_value,min_value),axis=1)
        cast_df.loc[:,"score0"] = new_scores

        x1 = cast_df.start.values
        x2 = cast_df.end.values
        y = cast_df.score0.values

        a_ind = y>=0
        ax1.bar(x1[a_ind], width=x1[a_ind]-x2[a_ind], height=y[a_ind], color="green", align="center")
        ax1.bar(x1[~a_ind], width=x1[~a_ind]-x2[~a_ind], height=y[~a_ind], color="red", align="center")
        max2 = max(y)
        min2 = min(y)
        if max_ < max2:
            max_ = max2
        if min_ > min2:
            min_ = min2

    axs[2].set_ylim(min_,max_)
    axs[1].set_ylim(min_,max_)
    axs[0].set_ylim(min_,max_)
    axs[0].set_xlim(start,end)
    axs[1].set_xlim(start,end)
    axs[2].set_xlim(start,end)
    fig.savefig("{}comps_{}.pdf".format(save_path,text))
    return True

#%%

# Get the plot
plot_compartments("chr10",0,130690000,"chr10")
