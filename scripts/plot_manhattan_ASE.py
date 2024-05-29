#%%

# There are 2 kind of plots here.
# First plots are number of CAST and S129 genes in TADs, X axis being genomic coordinates. These plots were not used.
# Second plots are position of genes in the chromosome, in a scatter plot. These are Figures 3a and SI 4a
import re
from adjustText import adjust_text
import pandas as pd
import pybedtools as pb
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *


ase_path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"
#ase_df = pd.read_csv(ase_path,sep="\t").iloc[:, 1:]
ase_df = pd.read_csv(ase_path,sep="\t")
ase_df = ase_df.drop_duplicates()
ase_df = ase_df[ase_df["chrom"].isin(chrs_)]


cast_TADs = root+"TADs.F123.as3NPs.CAST.bed"
s129_TADs = root+"TADs.F123.as3NPs.S129.bed"
type_ = "CAST_TADs" #CAST_TADs/S129_TADs

def get_ASE_in_TADs(tad_path,ase_df):
    TADs = pd.read_csv(tad_path,sep="\t",names=["chrom","start","end"])
    TADs["start"] = TADs.start - 25000
    TADs["end"] = TADs.end + 24999

    cutoff = 0.5
    plot_ = True

    cast = ase_df.query("TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05")
    s129 = ase_df.query("TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05")
    cast_pb = pb.BedTool.from_dataframe(cast)
    s129_pb = pb.BedTool.from_dataframe(s129)
    cast_TADs_pb = pb.BedTool.from_dataframe(TADs)

    aux = cast_TADs_pb.intersect(cast_pb,c=True)
    aux = aux.intersect(s129_pb,c=True).to_dataframe(names=["chrom","start","end","plus","minus"])
    aux = aux.query('plus != 0 or minus != 0')
    aux["TADsize"] = aux.end - aux.start
    if plot_:
        for chr_ in chrs_:
            aux_ = aux.query('chrom == "{}"'.format(chr_))
            fig,ax = plt.subplots(figsize=(40,4))
            ax.bar(aux_.start+(aux_.TADsize/2),aux_.plus,aux_.TADsize,color=cast_color)
            ax.bar(aux_.start+(aux_.TADsize/2),-1*(aux_.minus),aux_.TADsize,color=s129_color)
            aux_ = TADs.query('chrom == "{}"'.format(chr_))
            for a in aux_.start:
                ax.axvline(a,color="black",lw=0.5)
            for a in np.arange(-10,10):
                ax.axhline(a,color="grey",lw=0.5)
            ax.set_title(chr_)
            fig.savefig("{}figure2_manhattan_plot_ASE_in_TADs_{}_{}.pdf".format(save_path,chr_,type_))

    mix_ = aux.query('plus != 0 and minus != 0')
    minus_ = aux.query('plus == 0 and minus != 0')
    plus_ = aux.query('plus != 0 and minus == 0')
    print (len(mix_),len(plus_),len(minus_))

get_ASE_in_TADs(cast_TADs,ase_df)
#get_ASE_in_TADs(s129_TADs,ase_df)

# %%

#Bing examples
ex = "chr12"
aux_df = ase_df[ase_df.chrom == ex]
cast = aux_df.query("TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05")
s129 = aux_df.query("TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05")
exp = aux_df.query("TPM_transcript >= 1 and (log2foldchange > -1 and log2foldchange < 1) and number_SNPs >= 1")
exp_no_SNP = aux_df.query("TPM_transcript >= 1 and number_SNPs == 0")
print (len(cast),len(s129),len(exp),len( exp_no_SNP))
cast_y = cast.log2foldchange.values
cast_x = cast.TSS_start.values
s129_y = s129.log2foldchange.values
s129_x = s129.TSS_start.values
exp_y = exp.log2foldchange.values
exp_x = exp.TSS_start.values
exp_no_x = exp_no_SNP.TSS_start.values
exp_no_y = [0]*len(exp_no_x)


fig,ax = plt.subplots(figsize=(20,2))
ax.scatter(exp_x,exp_y,s=20, color="darkgreen",zorder=10)
ax.scatter(exp_no_x,exp_no_y,s=20, color="limegreen",zorder=11)
ax.scatter(cast_x,cast_y,s=20, color=cast_color,zorder=12)
ax.scatter(s129_x,s129_y,s=20, color=s129_color,zorder=13)

#First Hist
r = re.compile('^Hist.*')
r = re.compile('^fake.*')
highlight = list(filter(r.match, aux_df.gene_name))
highlight_df =aux_df[aux_df.gene_name.isin(highlight)]
y = highlight_df.log2foldchange.values.tolist() 
x = highlight_df.TSS_start.values.tolist()
for l in x:
    ax.axvline(l,color="crimson",lw=10)

#now ribos
r = re.compile('^Rps.*|^Rpl.*')
r = re.compile('^fake.*')
highlight = list(filter(r.match, aux_df.gene_name))
highlight_df =aux_df[aux_df.gene_name.isin(highlight)]
highlight_df.drop_duplicates(subset=["gene_name"],inplace=True)
y = highlight_df.log2foldchange.values.tolist() 
x = highlight_df.TSS_start.values.tolist()
names = highlight_df.gene_name.values.tolist()
texts = [ax.text(x[i], y[i], names[i], ha='center', va='center',zorder=14) for i in range(len(x))]
adjust_text(texts,autoalign=False, arrowprops=dict(arrowstyle='-', color='black'))
for l in x:
    ax.axvline(l,color="cornflowerblue",lw=10)
plt.savefig("{}manhattan_{}.pdf".format(save_path,ex),dpi=300)

# %%
fig,axs = plt.subplots(19,1,figsize=(20,10),sharex=True,sharey=True)
for i,chr_ in enumerate(chrs_):
    ax = axs[i]
    aux_df = ase_df[ase_df.chrom == chr_]
    cast = aux_df.query("TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05")
    s129 = aux_df.query("TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05")
    exp = aux_df.query("TPM_transcript >= 1 and (log2foldchange > -1 and log2foldchange < 1) and number_SNPs >= 1")
    exp_no_SNP = aux_df.query("TPM_transcript >= 1 and number_SNPs == 0")

    cast_y = cast.log2foldchange.values
    cast_x = cast.TSS_start.values
    s129_y = s129.log2foldchange.values
    s129_x = s129.TSS_start.values
    exp_y = exp.log2foldchange.values
    exp_x = exp.TSS_start.values
    exp_no_x = exp_no_SNP.TSS_start.values
    exp_no_y = [0]*len(exp_no_x)

    ax.scatter(exp_x,exp_y,s=5, color="darkgreen",zorder=10)
    ax.scatter(exp_no_x,exp_no_y,s=5, color="limegreen",zorder=11)
    ax.scatter(cast_x,cast_y,s=5, color=cast_color,zorder=12)
    ax.scatter(s129_x,s129_y,s=5, color=s129_color,zorder=13)

    #First Hist
    r = re.compile('^Hist.*')
    r = re.compile('^fake.*')
    highlight = list(filter(r.match, aux_df.gene_name))
    highlight_df =aux_df[aux_df.gene_name.isin(highlight)]
    y = highlight_df.log2foldchange.values.tolist() 
    x = highlight_df.TSS_start.values.tolist()
    for l in x:
        ax.axvline(l,color="crimson",lw=4)
    #names = highlight_df.gene_name.values.tolist()
    #now ribos
    r = re.compile('^Rps.*|^Rpl.*')
    r = re.compile('^fake.*')
    highlight = list(filter(r.match, aux_df.gene_name))
    highlight_df =aux_df[aux_df.gene_name.isin(highlight)]
    highlight_df.drop_duplicates(subset="gene_name",inplace=True)


    y = highlight_df.log2foldchange.values.tolist() 
    x = highlight_df.TSS_start.values.tolist()
    for l in x:
        ax.axvline(l,color="cornflowerblue",lw=4)
    #names = highlight_df.gene_name.values.tolist()
    #texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    #adjust_text(texts,autoalign=False, arrowprops=dict(arrowstyle='-', color='black'))
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("{}manhattan_all.pdf".format(save_path),dpi=300)
# %%
