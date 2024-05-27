#%%

# coded in 26 April, 2023

"""
the idea os to calculate the density of peaks for each feature by 100kbs (this is variable)
later, we can show some regions, or chromosomes with the density colorcoded.
Finally, we can do correlations of each feature against each other and plot that in a heatmap
"""

from functools import reduce
import operator
import seaborn as sns
from pandas.plotting import parallel_coordinates
import pybedtools as pb
import operator
from functools import reduce
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

#first load the features

ase_path = "{}300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed".format(root)
ase_df = pd.read_csv(ase_path,sep="\t")
ase_df = ase_df[ase_df.chrom.isin(chrs_)]

use_bins = False

if use_bins: 
    chrom_size_path = root+"mm10.chrom.sizes"
    bin_size = 200000
    a = pb.example_bedtool('a.bed')
    master_pb = pb.BedTool.window_maker(a,g=chrom_size_path,w=bin_size)

# instead of bins use TADs
else:
    tad_cast_path = root+"TADs.F123.as3NPs.CAST.bed"
    tad_s129_path = root+"TADs.F123.as3NPs.S129.bed"
    tad_path = root+"TADs.F123.as3NPs.bed"
    tads_cast_df = pd.read_csv(tad_cast_path,sep="\t",names=["chrom","start","end"])
    tads_s129_df = pd.read_csv(tad_s129_path,sep="\t",names=["chrom","start","end"])
    tads_df = pd.read_csv(tad_path,sep="\t",names=["chrom","start","end"])

    # since there is 50K space between TADs and we want to get genes, I will add -+25 to TADs
    tads_cast_df["start"] = tads_cast_df["start"]-25000
    tads_cast_df["end"] = tads_cast_df["end"]+25000
    tads_s129_df["start"] = tads_s129_df["start"]-25000
    tads_s129_df["end"] = tads_s129_df["end"]+25000
    tads_df["start"] = tads_df["start"]-25000
    tads_df["end"] = tads_df["end"]+25000

    master_pb = pb.BedTool.from_dataframe(tads_cast_df)


cast_pb = pb.BedTool.from_dataframe(ase_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum < 0.05')[["chrom","TSS_start","TSS_end"]])
s129_pb = pb.BedTool.from_dataframe(ase_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum < 0.05')[["chrom","TSS_start","TSS_end"]])




cast_rad21 = root+"Rad21.CAST.specific.peaks.bed"
s129_rad21 = root+"Rad21.S129.specific.peaks.bed"
common_rad21 = root+"Rad21.common.peaks.bed"

cast_atac = root+"F123_bing_data030420_genome_CAST.bed"
s129_atac = root+"F123_bing_data030420_genome_S129.bed"
common_atac = root+"common_peaks.bed"

cast_ctcf = root+"CTCF.CAST.specific.peaks.bed"
s129_ctcf = root+"CTCF.S129.specific.peaks.bed"
common_ctcf = root+"CTCF.common.peaks.bed"

cast_h3k4me3 = root+"H3K4me3.CAST.specific.peaks.bed"
s129_h3k4me3 = root+"H3K4me3.S129.specific.peaks.bed"
common_h3k4me3 = root+"H3K4me3.common.peaks.bed"

cast_h3k27ac = root+"H3K27Ac.CAST.specific.peaks.bed"
s129_h3k27ac = root+"H3K27Ac.S129.specific.peaks.bed"
common_h3k27ac = root+"H3K27Ac.common.peaks.bed"

H3K27me3 = root+"F123.K27me3.BCP_peaks.HM_mode.bed"
s5p = root+"F123.S5p.BCP_peaks.HM_mode.bed"
s7p = root+"F123.S7p.BCP_peaks.HM_mode.bed"

cast_atac_pb = pb.BedTool(cast_atac).sort()
s129_atac_pb = pb.BedTool(s129_atac).sort()
common_atac_pb = pb.BedTool(common_atac).sort()
cast_ctcf_pb = pb.BedTool(cast_ctcf).sort()
s129_ctcf_pb = pb.BedTool(s129_ctcf).sort()
common_ctcf_pb = pb.BedTool(common_ctcf).sort()
cast_rad21_pb = pb.BedTool(cast_rad21).sort()
s129_rad21_pb = pb.BedTool(s129_rad21).sort()
common_rad21_pb = pb.BedTool(common_rad21).sort()
cast_h3k4me3_pb = pb.BedTool(cast_h3k4me3).sort()
s129_h3k4me3_pb = pb.BedTool(s129_h3k4me3).sort()
common_h3k4me3_pb = pb.BedTool(common_h3k4me3).sort()
cast_h3k27ac_pb = pb.BedTool(cast_h3k27ac).sort()
s129_h3k27ac_pb = pb.BedTool(s129_h3k27ac).sort()
common_h3k27ac_pb = pb.BedTool(common_h3k27ac).sort()

s5_pb = pb.BedTool(s5p).sort()
s7_pb = pb.BedTool(s7p).sort()
h3K27me3_pb = pb.BedTool(H3K27me3).sort()


master_pb = master_pb.intersect(cast_pb,c=True)
master_pb = master_pb.intersect(s129_pb,c=True)

#master_pb = master_pb.intersect(snp_pb,c=True)
master_pb = master_pb.intersect(cast_atac_pb,c=True)
master_pb = master_pb.intersect(s129_atac_pb,c=True)
master_pb = master_pb.intersect(common_atac_pb,c=True)
master_pb = master_pb.intersect(cast_ctcf_pb,c=True)
master_pb = master_pb.intersect(s129_ctcf_pb,c=True)
master_pb = master_pb.intersect(common_ctcf_pb,c=True)
master_pb = master_pb.intersect(cast_rad21_pb,c=True)
master_pb = master_pb.intersect(s129_rad21_pb,c=True)
master_pb = master_pb.intersect(common_rad21_pb,c=True)
master_pb = master_pb.intersect(cast_h3k4me3_pb,c=True)
master_pb = master_pb.intersect(s129_h3k4me3_pb,c=True)
master_pb = master_pb.intersect(common_h3k4me3_pb,c=True)
master_pb = master_pb.intersect(cast_h3k27ac_pb,c=True)
master_pb = master_pb.intersect(s129_h3k27ac_pb,c=True)
master_pb = master_pb.intersect(common_h3k27ac_pb,c=True)

master_pb = master_pb.intersect(s5_pb,c=True)
master_pb = master_pb.intersect(s7_pb,c=True)
master_pb = master_pb.intersect(h3K27me3_pb,c=True)


cols = ["chrom","start","end","n_cast","n_s129","n_cast_atac","n_s129_atac","n_common_atac","n_cast_ctcf","n_s129_ctcf","n_common_ctcf","n_cast_rad21","n_s129_rad21","n_common_rad21",
        "n_cast_h3k4me3","n_s129_h3k4me3","n_common_h3k4me3","n_cast_h3k27ac","n_s129_h3k27ac","n_common_h3k27ac","n_s5","n_s7","n_h3k27me3"]
master_df = master_pb.to_dataframe(names=cols)

master_df = master_df[master_df.chrom.isin(chrs_)]
print (master_df)

#normalize by total peaks/numbers
def normalize_by_total(df,col):
    df[col] = df[col]/sum(df[col])
    return df
#for col in cols[3:]:
#    master_df = normalize_by_total(master_df,col)

#reorder
reordered_cols = ["chrom","start","end","n_cast","n_cast_atac","n_cast_h3k4me3","n_cast_h3k27ac","n_cast_ctcf","n_cast_rad21",
                  "n_s129","n_s129_atac","n_s129_h3k4me3","n_s129_h3k27ac","n_s129_ctcf","n_s129_rad21",
                  "n_common_atac","n_common_h3k4me3","n_common_h3k27ac","n_common_ctcf","n_common_rad21","n_s5","n_s7","n_h3k27me3"]
master_df = master_df[reordered_cols]



# %%

## plot figure
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

example = "chr8"

cmap = LinearSegmentedColormap.from_list('mycmap', [(0, '#000080'),(0.1,'#008080'),(0.4,'orange'),(1,'red')])

fig,ax = plt.subplots(figsize=(30,5))
plot_df = master_df.query('chrom == "{}"'.format(example))
plot_df = plot_df[reordered_cols[3:]]
plot_df = plot_df.T
plot_df = plot_df.replace(0,np.nan)
#ax.imshow(plot_df)
sns.heatmap(plot_df,cmap="viridis",vmax=7,vmin=1)
#sns.heatmap(plot_df,cmap=cmap,vmax=10,vmin=1)
#sns.heatmap(plot_df,cmap="Purples")
ax.tick_params('x', top=False, labeltop=False,bottom=False,labelbottom=False)
if use_bins:
    fig.savefig("{}{}_density_of_features_by_bins.pdf".format(save_path,example),dpi=300)
else:
    fig.savefig("{}{}_density_of_features_by_TADs.pdf".format(save_path,example),dpi=300)

# %%

#correlations and plot
from scipy.stats.stats import pearsonr,spearmanr 
spearman = False
size = len(reordered_cols[3:])
m = np.ones((size,size))
for i,x in enumerate(reordered_cols[3:]):
    for j,y in enumerate(reordered_cols[3:]):
        if j>i:
            if spearman:
                m[j][i] = spearmanr(master_df[x].values,master_df[y].values)[0]
            else:
                m[j][i] = pearsonr(master_df[x].values,master_df[y].values)[0]
        else:
            if j == i:
                pass
            else:
                m[j][i] = np.nan
        
fig,ax=plt.subplots(figsize=(8,8))
#im = ax.imshow(m,cmap="viridis")
im = ax.imshow(m,cmap="Blues")
#im = ax.imshow(m,cmap="GnBu")
# Minor ticks
ax.set_xticks(np.arange(-.5, size, 1), minor=True)
ax.set_yticks(np.arange(-.5, size, 1), minor=True)
#ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
#im = sns.clustermap(m)
fig.colorbar(im,ax=ax)
#ax.tick_params('x', top=True, labeltop=True,bottom=False,labelbottom=False)
ax.set_yticks(range(size))
ax.set_yticklabels(reordered_cols[3:])
ax.set_xticks(range(size))
ax.set_xticklabels(reordered_cols[3:],rotation=45)
#for i,x in enumerate(reordered_cols[3:]):
#    for j,y in enumerate(reordered_cols[3:]):
#        if j>i:
#            text = ax.text(i, j,("%.2f" % m[j,i]),
#                        ha="center", va="center", color="w")
if spearman:
    fig.savefig("{}feature_density_correlation_spearmanr.pdf".format(save_path,root),dpi=300)
else:
    fig.savefig("{}feature_density_correlation_pearsonr.pdf".format(save_path,root),dpi=300)


        

# %%
