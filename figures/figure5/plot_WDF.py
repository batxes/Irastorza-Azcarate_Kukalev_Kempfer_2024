#%%
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import numpy as np
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

cast = root+"wdf_cast_50k.bed"
s129 = root+"wdf_s129_50k.bed"

cast_color = "#611163"
s129_color = "#EA8D1F"
both_color = "#7A7A7A"
non_phased_color = "#09AA9E"
plot_bars = False
plot_lines = True

colors = [cast_color,s129_color]
def plot_wdfs(chr_,start,end,title,linewidth=1.0,fill=True):
    end = end + 50000
    cast_df = pd.read_csv(cast,sep="\t",names=["chrom","start","end","wdf"])
    s129_df = pd.read_csv(s129,sep="\t",names=["chrom","start","end","wdf"])
    fig,ax = plt.subplots(1,1,figsize=(18,2))
    fig.tight_layout()
    ys = []
    for i,cast_df in enumerate([cast_df,s129_df]):
        cast_df = cast_df.query('chrom == "{}" and start >= {} and end <= {}'.format(chr_,start,end))
        x1 = cast_df.start.values
        x2 = cast_df.end.values
        y = cast_df.wdf.values
        ys.append(y)
    if plot_bars:
        purple = []
        orange = []
        grey = []
        for (y1,y2) in zip(ys[0],ys[1]):
            if y1 > y2:
                purple.append(y1)
                orange.append(0)
                grey.append(0)
            elif y2 > y1:
                purple.append(0)
                orange.append(y2)
                grey.append(0)
            else:
                grey.append(y1)
                purple.append(0)
                orange.append(0)
                orange.append(0)

        #plot background first    
        ax.bar(x1,width=50000,height=purple,color=cast_color,align="center")
        ax.bar(x1,width=50000,height=orange,color=s129_color,align="center")
        ax.bar(x1,width=50000,height=grey,color=both_color,align="center")
        #now foreground
        grey = []
        for (y1,y2) in zip(ys[0],ys[1]):
            if y1 > y2:
                grey.append(y2)
            elif y2 > y1:
                grey.append(y1)
            else:
                grey.append(y1)
        ax.bar(x1,width=50000,height=grey,color=both_color,align="center")
    if plot_lines:
        ax.plot(x1,ys[0],linestyle='-',linewidth=linewidth,color=cast_color,label="CAST WDF")
        ax.plot(x1,ys[1],linestyle='-',linewidth=linewidth,color=s129_color,label="S129 WDF")
        if fill:
            ax.fill_between(x1, ys[0], ys[1],
                    where=(ys[0] < ys[1]),
                    alpha=0.20, color=s129_color, interpolate=True,label="more open in S129")
            ax.fill_between(x1, ys[0], ys[1],
                    where=(ys[0] > ys[1]),
                    alpha=0.20, color=cast_color, interpolate=True,label="more open in CAST")
    ax.legend()
    ax.set_ylabel("WDF")
    ax.set_xlabel("{}:{}-{}".format(chr_,start,end))
    ax.set_xlim(start,end-50000)
    ax.set_ylim(0.025,0.14)
    ax.grid()
    fig.savefig("{}wdf_{}.pdf".format(save_path,title))
    return True

plot_wdfs("chr7",139800000,141200000,"fig5d",3.0,True) 
#plot_wdfs("chr7",28200000,30300000,"fig5e",3.0,True) 
#plot_wdfs("chr17",33700000,35500000,"fig5h",3.0,True)



# %%
