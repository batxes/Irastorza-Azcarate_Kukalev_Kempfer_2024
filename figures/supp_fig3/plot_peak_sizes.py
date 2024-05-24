#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Supp Figure 3, Supp Figure3

import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

cast_atac = "{}F123_bing_data030420_genome_CAST.bed".format(root)
s129_atac = "{}F123_bing_data030420_genome_S129.bed".format(root)
common_atac = "{}common_peaks.bed".format(root)

cast_ctcf = "{}CTCF.CAST.specific.peaks.bed".format(root)
s129_ctcf = "{}CTCF.S129.specific.peaks.bed".format(root)
common_ctcf = "{}CTCF.common.peaks.bed".format(root)

cast_rad21 = "{}Rad21.CAST.specific.peaks.bed".format(root)
s129_rad21 = "{}Rad21.S129.specific.peaks.bed".format(root)
common_rad21 = "{}Rad21.common.peaks.bed".format(root)

cast_h3k4me3 = "{}H3K4me3.CAST.specific.peaks.bed".format(root)
s129_h3k4me3 = "{}H3K4me3.S129.specific.peaks.bed".format(root)
common_h3k4me3 = "{}H3K4me3.common.peaks.bed".format(root)

cast_h3k27ac = "{}H3K27Ac.CAST.specific.peaks.bed".format(root)
s129_h3k27ac = "{}H3K27Ac.S129.specific.peaks.bed".format(root)
common_h3k27ac = "{}H3K27Ac.common.peaks.bed".format(root)

s5p = "{}F123.S5p.BCP_peaks.HM_mode.bed".format(root)
s7p = "{}F123.S7p.BCP_peaks.HM_mode.bed".format(root)
H3K27me3 = "{}F123.K27me3.BCP_peaks.HM_mode.bed".format(root)

def read_ (path):
    df = pd.read_csv(path,sep="\t",names=["chrom","start","end","length"]) 
    df["length"] = df.end - df.start
    return df

def read_2 (path):
    df = pd.read_csv(path,sep="\t",names=["chrom","start","end","length","score"]) 
    return df

def plot_hist(values1,values2,values3,color1,color2,color3,title,bins):
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(8,8))
    _, bins, _ = ax1.hist(values1,color=color1,bins=bins,alpha=0.8)
    ax1.axvline(np.median(values1),color="black")
    ax1.axvline(np.mean(values1),color="red")
    ax2.hist(values2,color=color2,bins=bins,alpha=0.8)
    ax2.axvline(np.median(values2),color="black")
    ax2.axvline(np.mean(values2),color="red")
    ax3.hist(values3,color=color3,bins=bins,alpha=0.8)
    ax3.axvline(np.median(values3),color="black")
    ax3.axvline(np.mean(values3),color="red")
    ax1.set_xlim([-100,2000])
    ax3.set_xlabel("Size in bp")
    ax1.set_title("CAST {} peaks. Median: {} Mean: {}".format(title,np.median(values1),np.mean(values1)))
    ax2.set_title("S129 {} peaks Median: {} Mean: {}".format(title,np.median(values2),np.mean(values2)))
    ax3.set_title("Common {} peaks Median: {} Mean: {}".format(title,np.median(values3),np.mean(values3)))
    fig.savefig("{}{}_peak_size.pdf".format(save_path,title),dpi=300)

def plot_hist2(values1,values2,values3,color1,color2,color3,title,bins):
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(8,8))
    _, bins, _ = ax1.hist(values1,color=color1,bins=bins,alpha=0.8)
    ax1.axvline(np.median(values1),color="black")
    ax1.axvline(np.mean(values1),color="red")
    ax2.hist(values2,color=color2,bins=bins,alpha=0.8)
    ax2.axvline(np.median(values2),color="black")
    ax2.axvline(np.mean(values2),color="red")
    ax3.hist(values3,color=color3,bins=bins,alpha=0.8)
    ax3.axvline(np.median(values3),color="black")
    ax3.axvline(np.mean(values3),color="red")
    ax1.set_xlim([-100,20000])
    ax3.set_xlabel("Size in bp")
    ax1.set_title("S5p {} peaks Median: {} Mean: {}".format(title,np.median(values1),np.mean(values1)))
    ax2.set_title("S7p {} peaks Median: {} Mean: {}".format(title,np.median(values2),np.mean(values2)))
    ax3.set_title("H3K27me3 {} peaks Median: {} Mean: {}".format(title,np.median(values3),np.mean(values3)))
    fig.savefig("{}{}_peak_size.pdf".format(save_path,title),dpi=300)


df1 = read_(cast_atac)
df2 = read_(s129_atac)
df3 = read_(common_atac)
df = pd.concat([df1,df2,df3])
print ("ATAC median: {}, mean: {}".format(np.median(df.length),np.mean(df.length)))
plot_hist(df1.length,df2.length,df3.length,cast_color,s129_color,both_color,"ATAC",1000)

df1 = read_(cast_ctcf)
df2 = read_(s129_ctcf)
df3 = read_(common_ctcf)
df = pd.concat([df1,df2,df3])
print ("CTCF median: {}, mean: {}".format(np.median(df.length),np.mean(df.length)))
plot_hist(df1.length,df2.length,df3.length,cast_color,s129_color,both_color,"CTCF",100)

df1 = read_(cast_h3k4me3)
df2 = read_(s129_h3k4me3)
df3 = read_(common_h3k4me3)
df = pd.concat([df1,df2,df3])
print ("K4m3 median: {}, mean: {}".format(np.median(df.length),np.mean(df.length)))
plot_hist(df1.length,df2.length,df3.length,cast_color,s129_color,both_color,"H3K4me3",200)

df1 = read_(cast_h3k27ac)
df2 = read_(s129_h3k27ac)
df3 = read_(common_h3k27ac)
df = pd.concat([df1,df2,df3])
print ("K27ac median: {}, mean: {}".format(np.median(df.length),np.mean(df.length)))
plot_hist(df1.length,df2.length,df3.length,cast_color,s129_color,both_color,"H3k27ac",200)

df1 = read_(cast_rad21)
df2 = read_(s129_rad21)
df3 = read_(common_rad21)
df = pd.concat([df1,df2,df3])
print ("rad21 median: {}, mean: {}".format(np.median(df.length),np.mean(df.length)))
plot_hist(df1.length,df2.length,df3.length,cast_color,s129_color,both_color,"Rad21",100)

df1 = read_2(s5p)
df2 = read_2(s7p)
df3 = read_2(H3K27me3)
df = pd.concat([df1,df2,df3])
print ("rad21 median: {}, mean: {}".format(np.median(df.length),np.mean(df.length)))
plot_hist2(df1.length,df2.length,df3.length,"lightblue","blue","red","POL2",2500)
# %%
