#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
date_ = datetime.datetime.now()
date = "{}-{}-{}".format(date_.day,date_.month,date_.year)
print (date)
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

# Figure 2g.
files_ = []

common_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_ATAC_Bing_shared_peaks.bed"
files_.append(common_atac)
nonphased_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_ATAC_Bing_non_phased_peaks.bed"
files_.append(nonphased_atac)
cast_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_bing_data030420_genome_CAST.bed"
files_.append(cast_atac)
s129_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_bing_data030420_genome_S129.bed"
files_.append(s129_atac)

common_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/CTCF.common.peaks.bed"
files_.append(common_ctcf)
nonphased_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/CTCF.non_phased.peaks.nosex.bed"
files_.append(nonphased_ctcf)
cast_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/CTCF.CAST.specific.peaks.bed"
files_.append(cast_ctcf)
s129_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/CTCF.S129.specific.peaks.bed"
files_.append(s129_ctcf)

common_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K4me3.common.peaks.bed"
files_.append(common_h3k4me3)
nonphased_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/H3K4me3.non_phased.peaks.nosex.bed"
files_.append(nonphased_h3k4me3)
cast_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K4me3.CAST.specific.peaks.bed"
files_.append(cast_h3k4me3)
s129_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K4me3.S129.specific.peaks.bed"
files_.append(s129_h3k4me3)

common_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K27Ac.common.peaks.bed"
files_.append(common_h3k27ac)
nonphased_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/H3K27Ac.non_phased.peaks.nosex.bed"
files_.append(nonphased_h3k27ac)
cast_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K27Ac.CAST.specific.peaks.bed"
files_.append(cast_h3k27ac)
s129_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K27Ac.S129.specific.peaks.bed"
files_.append(s129_h3k27ac)

common_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/Rad21.common.peaks.bed"
files_.append(common_rad21)
nonphased_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/Rad21.common.nonphased.peaks.nosex.bed"
files_.append(nonphased_rad21)
cast_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/Rad21.CAST.specific.peaks.bed"
files_.append(cast_rad21)
s129_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/Rad21.S129.specific.peaks.bed"
files_.append(s129_rad21)

sizes = []
for file_ in files_:
    df = pd.read_csv(file_,sep="\t",header=None)
    sizes.append(df.shape[0])

color=[both_color,non_phased_color,cast_color,s129_color,both_color,non_phased_color,cast_color,s129_color,both_color,non_phased_color,cast_color,s129_color,both_color,non_phased_color,cast_color,s129_color,both_color,non_phased_color,cast_color,s129_color,]

s = [[sizes[0],sizes[1],sizes[2],sizes[3]],
     [sizes[4],sizes[5],sizes[6],sizes[7]],
      [sizes[8],sizes[9],sizes[10],sizes[11]],
       [sizes[12],sizes[13],sizes[14],sizes[15]],
       [sizes[16],sizes[17],sizes[18],sizes[19]]]
for i,ex in enumerate(["ATAC","CTCF","H3K4me3","H3K27ac","RAD21"]):
    fig,ax = plt.subplots(1,1,figsize=(4,4))
    counts = pd.Series([s[i][0],s[i][1],s[i][2],s[i][3]],
                    index=[ex+"\n"+str(s[i][0]),ex+"\n"+str(s[i][1]),ex+"\n"+str(s[i][2]),ex+"\n"+str(s[i][3])] 
                                    )
    counts.plot(kind='pie', fontsize=17, colors=[both_color,non_phased_color,cast_color,s129_color], 
                explode=[0,0.05,0.1,0.15],shadow=False,frame=False,autopct='%1.1f%%'
                )
    plt.axis('equal')
    plt.ylabel('')
    plt.legend(labels=counts.index, loc="best")
    plt.savefig("{}pie_plots_{}.pdf".format(save_path,ex),dpi=300)
    plt.show()


# %%



# %%
