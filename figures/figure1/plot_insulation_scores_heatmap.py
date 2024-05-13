#%%

"""
Script that plots unphased, cast and s129 insulation score heatmap, at different resolutions, for the region given as input.
This plot was used in Figure 1, SI figure 2.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

def plot_IS(chr_,from_,to_,add_nonphased=True):
    if add_nonphased:
        data3 = "{}F123.as3NPs.IS.scores.csv".format(root)
        df3 = pd.read_csv(data3,sep="\t")
        df3 = df3.drop(columns=["as3NPs_ins200K","as3NPs_ins300K"])
    data2 = "{}210910.F123.as.3NPs.curated.S129.insulation.scores.at50Kb.table".format(root)
    data1 = "{}210910.F123.as.3NPs.curated.CAST.insulation.scores.at50Kb.table".format(root)
    df1 = pd.read_csv(data1,sep="\t")
    df1 = df1.drop(columns=["ins100K","ins200K","ins900K","ins1Mb"])
    df2 = pd.read_csv(data2,sep="\t")
    df2 = df2.drop(columns=["ins100K","ins200K","ins900K","ins1Mb"])

    vmin,vmax = 0,0

    if add_nonphased:
        fig,axs = plt.subplots(3,1,sharex=True,figsize=(18,4))
        ll = [df3,df1,df2]
    else:
        fig,axs = plt.subplots(2,1,sharex=True,figsize=(18,4))
        ll = [df1,df2]
    for n,df in enumerate(ll):
        df = df.query('chrom == "{}" and start >= {} and stop <= {}'.format(chr_,from_,to_))
        df = df.iloc[:,3:].T
        df = df.iloc[::-1]

        #df_aux = df.dropna()
        df_aux = df.copy()
        all_values = []
        for sub_list in df_aux.values.tolist():
            all_values += sub_list
        vmin_aux = min(all_values)
        vmax_aux = max(all_values)
        if vmin_aux < vmin:
            vmin = vmin_aux
        if vmax_aux > vmax:
            vmax = vmax_aux

    print (vmax,vmin)
    vmin = -0.3
    vmax = 0.3
    for n,df in enumerate(ll):
        df = df.query('chrom == "{}" and start >= {} and stop <= {}'.format(chr_,from_,to_))
        df = df.iloc[:,3:].T
        df = df.iloc[::-1]
        im = axs[n].imshow(df,aspect='auto',vmin=vmin,vmax=vmax,cmap="RdBu_r")

    plt.savefig("{}IS_heatmap_{}:{}-{}.pdf".format(save_path,chr_,from_,to_),dpi=300)
    plt.show()

    fig,(ax,ax2) = plt.subplots(2,1,figsize=(40,4))
    df = df1.query('chrom == "{}" and start >= {} and stop <= {}'.format(chr_,from_,to_))
    df = df.replace(np.nan,0)
    df = df.iloc[:,3:].T
    df = df.iloc[::-1]
    df_2 = df2.query('chrom == "{}" and start >= {} and stop <= {}'.format(chr_,from_,to_))
    df_2 = df_2.replace(np.nan,0)
    df_2 = df_2.iloc[:,3:].T
    df_2 = df_2.iloc[::-1]

    pat = df.iloc[0].values
    mat = df_2.iloc[0].values
    diff = [x-y for x,y in zip(pat,mat)]
    ax.plot(diff)
    ax.axhline(0,ls="--",color="black")
    ax.set_ylabel("w = 800")
    ax.set_ylim(-0.2,0.2)
    pat = df.iloc[4].values
    mat = df_2.iloc[4].values
    diff = [x-y for x,y in zip(pat,mat)]
    ax2.plot(diff)
    ax2.axhline(0,ls="--",color="black")
    ax2.set_ylabel("w = 400")
    ax2.set_ylim(-0.2,0.2)

#%%

# Get the plot.
plot_IS("chr10",23000000,31000000,True)


# %%
