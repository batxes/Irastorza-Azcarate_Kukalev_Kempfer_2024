#%%
import pybedtools as pb
import re
from os import listdir
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

# This script plots different figures related to compartments and their overlap with features and genes

def sort_df(df, column_idx, key):
    '''Takes a dataframe, a column index and a custom function for sorting, 
    returns a dataframe sorted by that column using that function'''

    col = df.loc[:,column_idx]
    df = df.loc[[i[1] for i in sorted(zip(col,range(len(col))), key=key)]]
    return df

# function that merges two dataframes using intersections
# both dataframes must contain chrom start end in the beginning
# overlap is for master_df. Values in master_df will overlap anything from compartments. From all the intersections
# the promoter will be taken, depending on strand
# we return the dataframe without index

def merge_dfs_with_intersect(master_df,new_df,slop=0):
    print (master_df.columns)
    print (new_df.columns)
    print ("number of genes before: {}".format(len(master_df)))
    master_df = master_df.sort_values(by=['chrom', 'start'])
    if slop != 0:
        genome = root+"mm10.chrom.sizes"
        master_pb = pb.BedTool.from_dataframe(master_df).slop(g=genome,b=slop)
    else:
        master_pb = pb.BedTool.from_dataframe(master_df)
    columns_new = list(new_df.columns.values)
    columns_new[0] = "chrom_new"
    columns_new[1] = "start_new"
    columns_new[2] = "end_new"
    columns_master = list(master_df.columns.values)
    columns = columns_new+columns_master
    new_pb = pb.BedTool.from_dataframe(new_df)
    aux_pb = new_pb.intersect(master_pb,wb=True) 
    # now, remove chrom_new, start_new and end_new, and put the rest of the joined dataframe to the end of the master df
    aux_df = aux_pb.to_dataframe(names=columns,usecols=columns_master+columns_new[3:])
    aux_df = aux_df[columns_master+columns_new[3:]]
    if slop != 0:
        aux_df["start"] = aux_df["start"]+slop
        aux_df["end"] = aux_df["end"]-slop
    # take the promoter always, depending on strand
    aux_df_pos_strand = aux_df.query('strand == "+"')
    aux_df_neg_strand = aux_df.query('strand == "-"')
    aux_df_pos_strand = aux_df_pos_strand.groupby(['chrom','start','end'],as_index=False).first()
    aux_df_neg_strand = aux_df_neg_strand.groupby(['chrom','start','end'],as_index=False).last()
    aux_df_pos_strand = aux_df_pos_strand.set_index(['chrom','start','end'])
    aux_df_neg_strand = aux_df_neg_strand.set_index(['chrom','start','end'])
    aux_df_pos_strand = aux_df_pos_strand[columns_new[3:]]
    aux_df_neg_strand = aux_df_neg_strand[columns_new[3:]]
    aux_df = pd.concat([aux_df_pos_strand,aux_df_neg_strand])
    master_df = master_df.set_index(['chrom','start','end'])

    master_df = master_df.merge(aux_df,on=['chrom','start','end'],how="outer")
    print ("number of genes after: {}".format(len(master_df)))
    return master_df.reset_index()

if __name__ == '__main__':

    #### set up ASE file
    ase_path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"
    ase_df = pd.read_csv(ase_path, sep="\t")
    master_df = ase_df.rename(columns={"TSS_start":"start","TSS_end":"end"})

    #set up compartment DFs and merge them  
    print ("Adding compartments, 3DF")
    np_ = "All_as_3NP"
    path = root+"compartments/"
    ab_files_all = [f for f in listdir(path) if "no_cleaning" in f]
    ab_files_cast = [f for f in listdir(path) if "CAST" in f]
    ab_files_s129 = [f for f in listdir(path) if "S129" in f]
    names = ["PCA_score_nonphased","PCA_score_cast","PCA_score_s129"]
    namesAB = ["AB_nonphased","AB_cast","AB_s129"]
    ab_files = [ab_files_all,ab_files_cast,ab_files_s129]
    final_comp_dfs = []
    for n,ab_file in enumerate(ab_files):
        comp_dfs = []
        for chr_ in chrs_:
            r = re.compile(".*{}\.AB\.txt".format(chr_))
            try:
                abfile = list(filter(r.match, ab_file))[0]
            except:
                continue
            comp_df = pd.read_csv(path+abfile,sep="\t",usecols=["chrom","start","end","ABchosen","chosen_score"])
            comp_df = comp_df.rename(columns={"ABchosen": namesAB[n],"chosen_score": names[n]})
            comp_dfs.append(comp_df)
        #concatenate all chromosomes and then merge
        master_df = merge_dfs_with_intersect(master_df,pd.concat(comp_dfs))

    master_df = master_df.replace('.', np.NaN)
    master_df.to_csv("{}Master_table_genes_{}_test.csv".format(save_path,np_))
    print ("Master table written!")


# %%

# create a scatter plot with both PCA scores in the axes

font = {'family': 'sansserif',
        'color':  'black',
        'weight': 'normal',
        'size': 25,
        }

font1 = {'family': 'sansserif',
        'color':  cast_color,
        'weight': 'normal',
        'size': 20,
        }
font2 = {'family': 'sansserif',
        'color':  non_phased_color,
        'weight': 'normal',
        'size': 20,
        }
font3 = {'family': 'sansserif',
        'color':  s129_color,
        'weight': 'normal',
        'size': 20,
        }


cast_genes_df = master_df.query('TPM_transcript >= 1 & log2foldchange >= 1 & p_adj_sum <= 0.05')
n_AA_cast_genes = len(cast_genes_df.query('AB_cast == "A" and AB_s129 == "A"'))
n_BB_cast_genes = len(cast_genes_df.query('AB_cast == "B" and AB_s129 == "B"'))
n_AB_cast_genes = len(cast_genes_df.query('AB_cast == "A" and AB_s129 == "B"'))
n_BA_cast_genes = len(cast_genes_df.query('AB_cast == "B" and AB_s129 == "A"'))
cast_genes_paternal_PCA = cast_genes_df["PCA_score_cast"].values.tolist()
cast_genes_maternal_PCA = cast_genes_df["PCA_score_s129"].values.tolist()

s129_genes_df = master_df.query('TPM_transcript >= 1 & log2foldchange <= -1 & p_adj_sum <= 0.05')
n_AA_s129_genes = len(s129_genes_df.query('AB_cast == "A" and AB_s129 == "A"'))
n_BB_s129_genes = len(s129_genes_df.query('AB_cast == "B" and AB_s129 == "B"'))
n_AB_s129_genes = len(s129_genes_df.query('AB_cast == "A" and AB_s129 == "B"'))
n_BA_s129_genes = len(s129_genes_df.query('AB_cast == "B" and AB_s129 == "A"'))
s129_genes_paternal_PCA = s129_genes_df["PCA_score_cast"].values.tolist()
s129_genes_maternal_PCA = s129_genes_df["PCA_score_s129"].values.tolist()

mono_cast_genes_df = master_df.query('TPM_transcript >= 1 & J129count_sum == 0 & p_adj_sum <= 0.05')
n_AA_mono_cast_genes = len(mono_cast_genes_df.query('AB_cast == "A" and AB_s129 == "A"'))
n_BB_mono_cast_genes = len(mono_cast_genes_df.query('AB_cast == "B" and AB_s129 == "B"'))
n_AB_mono_cast_genes = len(mono_cast_genes_df.query('AB_cast == "A" and AB_s129 == "B"'))
n_BA_mono_cast_genes = len(mono_cast_genes_df.query('AB_cast == "B" and AB_s129 == "A"'))
mono_cast_genes_paternal_PCA = mono_cast_genes_df["PCA_score_cast"].values.tolist()
mono_cast_genes_maternal_PCA = mono_cast_genes_df["PCA_score_s129"].values.tolist()

mono_s129_genes_df = master_df.query('TPM_transcript >= 1 & CASTcount_sum == 0 & p_adj_sum <= 0.05')
n_AA_mono_s129_genes = len(mono_s129_genes_df.query('AB_cast == "A" and AB_s129 == "A"'))
n_BB_mono_s129_genes = len(mono_s129_genes_df.query('AB_cast == "B" and AB_s129 == "B"'))
n_AB_mono_s129_genes = len(mono_s129_genes_df.query('AB_cast == "A" and AB_s129 == "B"'))
n_BA_mono_s129_genes = len(mono_s129_genes_df.query('AB_cast == "B" and AB_s129 == "A"'))
mono_s129_genes_paternal_PCA = mono_s129_genes_df["PCA_score_cast"].values.tolist()
mono_s129_genes_maternal_PCA = mono_s129_genes_df["PCA_score_s129"].values.tolist()

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.scatter(cast_genes_paternal_PCA,cast_genes_maternal_PCA,color=cast_color,label="CAST genes")
ax.scatter(mono_cast_genes_paternal_PCA,mono_cast_genes_maternal_PCA,color=non_phased_color,label="Monoallelic")
ax.text(4,5,n_AA_cast_genes,fontdict=font1)
ax.text(4,-5,n_AB_cast_genes,fontdict=font1)
ax.text(-4,5,n_BA_cast_genes,fontdict=font1)
ax.text(-4,-5,n_BB_cast_genes,fontdict=font1)
ax.text(5,5,n_AA_mono_cast_genes,fontdict=font2)
ax.text(5,-5,n_AB_mono_cast_genes,fontdict=font2)
ax.text(-5,5,n_BA_mono_cast_genes,fontdict=font2)
ax.text(-5,-5,n_BB_mono_cast_genes,fontdict=font2)
ax.set_xlim([-6,6])
ax.set_ylim([-6,6])
ax.axvline(0,color="grey")
ax.axhline(0,color="grey")
ax.set_ylabel("compartment scores (S129)", fontdict=font)
ax.set_xlabel("compartment scores (CAST)", fontdict=font)
plt.legend()
fig.savefig("{}CAST_genes_in_compartments_{}.pdf".format(save_path,np_),dpi=300,bbox_inches='tight')

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.scatter(s129_genes_paternal_PCA,s129_genes_maternal_PCA,color=s129_color,label="S129 genes")
ax.scatter(mono_s129_genes_paternal_PCA,mono_s129_genes_maternal_PCA,color=non_phased_color,label="Monoallelic")
ax.text(4,5,n_AA_s129_genes,fontdict=font3)
ax.text(4,-5,n_AB_s129_genes,fontdict=font3)
ax.text(-4,5,n_BA_s129_genes,fontdict=font3)
ax.text(-4,-5,n_BB_s129_genes,fontdict=font3)
ax.text(5,5,n_AA_mono_s129_genes,fontdict=font2)
ax.text(5,-5,n_AB_mono_s129_genes,fontdict=font2)
ax.text(-5,5,n_BA_mono_s129_genes,fontdict=font2)
ax.text(-5,-5,n_BB_mono_s129_genes,fontdict=font2)
ax.set_xlim([-6,6])
ax.set_ylim([-6,6])
ax.axvline(0,color="grey")
ax.axhline(0,color="grey")
ax.set_ylabel("compartment scores (S129)", fontdict=font)
ax.set_xlabel("compartment scores (CAST)", fontdict=font)
plt.legend()
fig.savefig("{}S129_genes_in_compartments_{}.pdf".format(save_path,np_),dpi=300,bbox_inches='tight')

#%%

# stacked bar plot with genes in compartments
exp_no_ASE_SNP_df = master_df.query('TPM_transcript >= 1  & (log2foldchange > -1 & log2foldchange < 1 ) & number_SNPs >= 1')
n_AA_exp_no_ase_SNP_genes = len(exp_no_ASE_SNP_df.query('AB_cast == "A" and AB_s129 == "A"'))
n_BB_exp_no_ase_SNP_genes = len(exp_no_ASE_SNP_df.query('AB_cast == "B" and AB_s129 == "B"'))
n_AB_exp_no_ase_SNP_genes = len(exp_no_ASE_SNP_df.query('AB_cast == "A" and AB_s129 == "B"'))
n_BA_exp_no_ase_SNP_genes = len(exp_no_ASE_SNP_df.query('AB_cast == "B" and AB_s129 == "A"'))

exp_no_ASE_no_SNP_df = master_df.query('TPM_transcript >= 1 & number_SNPs == 0')
n_AA_exp_no_ase_no_SNP_genes = len(exp_no_ASE_no_SNP_df.query('AB_cast == "A" and AB_s129 == "A"'))
n_BB_exp_no_ase_no_SNP_genes = len(exp_no_ASE_no_SNP_df.query('AB_cast == "B" and AB_s129 == "B"'))
n_AB_exp_no_ase_no_SNP_genes = len(exp_no_ASE_no_SNP_df.query('AB_cast == "A" and AB_s129 == "B"'))
n_BA_exp_no_ase_no_SNP_genes = len(exp_no_ASE_no_SNP_df.query('AB_cast == "B" and AB_s129 == "A"'))

no_exp_df = master_df.query('TPM_transcript < 1')
n_AA_no_exp_genes = len(no_exp_df.query('AB_cast == "A" and AB_s129 == "A"'))
n_BB_no_exp_genes = len(no_exp_df.query('AB_cast == "B" and AB_s129 == "B"'))
n_AB_no_exp_genes = len(no_exp_df.query('AB_cast == "A" and AB_s129 == "B"'))
n_BA_no_exp_genes = len(no_exp_df.query('AB_cast == "B" and AB_s129 == "A"'))

fig,ax = plt.subplots(figsize=(2,4))
labels = ['A/A', 'B/B', 'A/B', 'B/A']
width=0.8
cast_values = [n_AA_cast_genes,n_BB_cast_genes,n_AB_cast_genes,n_BA_cast_genes]
s129_values = [n_AA_s129_genes,n_BB_s129_genes,n_AB_s129_genes,n_BA_s129_genes]
exp_no_ASE_SNP_values = [n_AA_exp_no_ase_SNP_genes,n_BB_exp_no_ase_SNP_genes,n_AB_exp_no_ase_SNP_genes,n_BA_exp_no_ase_SNP_genes]
exp_no_ASE_no_SNP_values = [n_AA_exp_no_ase_no_SNP_genes,n_BB_exp_no_ase_no_SNP_genes,n_AB_exp_no_ase_no_SNP_genes,n_BA_exp_no_ase_no_SNP_genes]
no_exp_values = [n_AA_no_exp_genes,n_BB_no_exp_genes,n_AB_no_exp_genes,n_BA_no_exp_genes]

print (cast_values)
print (s129_values)
print(exp_no_ASE_SNP_values)
print(exp_no_ASE_no_SNP_values)
print(no_exp_values)

ax.bar(labels,cast_values,width=width,label="CAST",color=cast_color)
ax.bar(labels,s129_values,width=width,label="S129",bottom=cast_values,color=s129_color)
summup = [x+y for x,y in zip(cast_values, s129_values)]
ax.bar(labels,exp_no_ASE_SNP_values,width=width,label="exp_no_ASE_w_SNP",bottom=summup,color="darkgreen")
summup = [x+y for x,y in zip(summup, exp_no_ASE_SNP_values)]
ax.bar(labels,exp_no_ASE_no_SNP_values,width=width,label="exp_no_ASE_wo_SNP",bottom=summup,color="lightgreen")
summup = [x+y for x,y in zip(summup, exp_no_ASE_no_SNP_values)]
ax.bar(labels,no_exp_values,width=width,label="not_exp",bottom=summup,color="grey")
ax.legend()
ax.set_ylabel("Number of genes")
plt.savefig(save_path+"stacked_bar_comps.eps",dpi=300)

total_n_genes = [x+y for x,y in zip(summup, no_exp_values)]
fig,ax = plt.subplots(figsize=(2,4))
cast_per = [x/y for x,y in zip(cast_values, total_n_genes)]
ax.bar(labels,cast_per,width=width,label="CAST",color=cast_color)

s129_per = [x/y for x,y in zip(s129_values, total_n_genes)]
ax.bar(labels,s129_per,width=width,label="S129",bottom=cast_per,color=s129_color)

exp_w_SNP_per = [x/y for x,y in zip(exp_no_ASE_SNP_values, total_n_genes)]
summup = [x+y for x,y in zip(cast_per, s129_per)]
ax.bar(labels,exp_w_SNP_per,width=width,label="exp_no_ASE_w_SNP",bottom=summup,color="darkgreen")

exp_wo_SNP_per = [x/y for x,y in zip(exp_no_ASE_no_SNP_values, total_n_genes)]
summup = [x+y for x,y in zip(summup, exp_w_SNP_per)]
ax.bar(labels,exp_wo_SNP_per,width=width,label="exp_no_ASE_wo_SNP",bottom=summup,color="lightgreen")

not_exp_per = [x/y for x,y in zip(no_exp_values, total_n_genes)]
summup = [x+y for x,y in zip(summup, exp_wo_SNP_per)]
ax.bar(labels,not_exp_per,width=width,label="not_exp",bottom=summup,color="grey")
ax.legend()
ax.set_ylabel("Percentage of genes")
plt.savefig(save_path+"stacked_bar_comps2.eps",dpi=300)

# %%

# Figure2 plots

font = {'family': 'sansserif',
        'color':  'black',
        'weight': 'normal',
        'size': 8,
        }
# A compartment in cast, take all genes and get TPM
cast_a_exp = master_df.query('AB_cast == "A"')
cast_ab_exp = master_df.query('AB_cast == "A" & AB_s129 == "B"')
cast_b_exp = master_df.query('AB_cast == "B"')

fig,ax = plt.subplots(figsize=(2,8))
#bp1 = ax.boxplot([cast_a_exp["TPM_transcript"].values.tolist(),cast_ab_exp["TPM_transcript"].values.tolist(),cast_b_exp["TPM_transcript"].values.tolist()],showfliers=False,patch_artist=True,labels=["A","specific A","B"])
bp1 = ax.boxplot([cast_a_exp["TPM_transcript"].values.tolist(),cast_ab_exp["TPM_transcript"].values.tolist(),cast_b_exp["TPM_transcript"].values.tolist()],showfliers=False,patch_artist=True,labels=["A","specific A","B"])
colors = [cast_color, 'mediumpurple','grey']
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp1["medians"],color="black")
ax.set_xlabel("CAST compartments")
ax.set_ylabel("TPM")
ax.text(0.5, 18,"n: {}".format(len(cast_a_exp)) , fontdict=font)
ax.text(1.5, 18,"n: {}".format(len(cast_ab_exp)) , fontdict=font)
ax.text(2.5, 18,"n: {}".format(len(cast_b_exp)) , fontdict=font)
fig.savefig("{}CAST_compartments_TPM_{}.png".format(save_path,np_),dpi=300,bbox_inches='tight')


s129_a_exp = master_df.query('AB_s129 == "A"')
s129_ab_exp = master_df.query('AB_s129 == "A" & AB_cast == "B"')
s129_b_exp = master_df.query('AB_s129 == "B"')

fig,ax = plt.subplots(figsize=(2,8))
bp1 = ax.boxplot([s129_a_exp["TPM_transcript"].values.tolist(),s129_ab_exp["TPM_transcript"].values.tolist(),s129_b_exp["TPM_transcript"].values.tolist()],showfliers=False,patch_artist=True,labels=["A","specific A","B"])
colors = [s129_color, 'wheat','grey']
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp1["medians"],color="black")
ax.set_xlabel("S129 compartments")
ax.set_ylabel("TPM")
ax.text(0.5, 18,"n: {}".format(len(s129_a_exp)) , fontdict=font)
ax.text(1.5, 18,"n: {}".format(len(s129_ab_exp)) , fontdict=font)
ax.text(2.5, 18,"n: {}".format(len(s129_b_exp)) , fontdict=font)
fig.savefig("{}S129_compartments_TPM_{}.png".format(save_path,np_),dpi=300,bbox_inches='tight')

#%%

labels = ["non_phased","CAST","S129"]
non_phased_comp = []
cast_comp = []
s129_comp = []

for n,v in enumerate(labels):
    for chr_ in chrs_:
        if n == 0:
            df = pd.read_csv("{}GAM_F123_ALL_as_3NP_no_cleaning.{}.AB.txt".format(path,chr_),sep="\t")
            non_phased_comp.append(df)
        elif n == 1:
            df = pd.read_csv("{}GAM_F123_ALL_as_3NP_CAST.{}.AB.txt".format(path,chr_),sep="\t")
            cast_comp.append(df)
        else:
            df = pd.read_csv("{}GAM_F123_ALL_as_3NP_S129.{}.AB.txt".format(path,chr_),sep="\t")
            s129_comp.append(df)

non_phased_AB = pd.concat(non_phased_comp)
cast_AB = pd.concat(cast_comp)
s129_AB = pd.concat(s129_comp)

non_phased_AB_pb = pb.BedTool.from_dataframe(non_phased_AB)
cast_AB_pb = pb.BedTool.from_dataframe(cast_AB)
s129_AB_pb = pb.BedTool.from_dataframe(s129_AB)
import numpy
import matplotlib.pyplot as plt
cast_AB = cast_AB.rename(columns={"ABchosen": "AB_cast"})
s129_AB = s129_AB.rename(columns={"ABchosen": "AB_s129"})
try:
    cast_AB.set_index(['chrom','start','end'], inplace = True)
    s129_AB.set_index(['chrom','start','end'], inplace = True)
except:
    pass
cast_s129_AB = pd.concat([cast_AB,s129_AB],axis=1, join='outer')
cast_s129_AB.reset_index(inplace=True)

cast_A_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_cast == "A"')[['chrom','start','end']])
cast_B_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_cast == "B"')[['chrom','start','end']])
cast_A_s129_B_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_cast == "A" & AB_s129 == "B"')[['chrom','start','end']])
cast_B_s129_A_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_cast == "B" & AB_s129 == "A"')[['chrom','start','end']])
cast_A_s129_A_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_cast == "A" & AB_s129 == "A"')[['chrom','start','end']])
cast_B_s129_B_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_cast == "B" & AB_s129 == "B"')[['chrom','start','end']])
s129_A_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_s129 == "A"')[['chrom','start','end']])
s129_B_pb = pb.BedTool.from_dataframe(cast_s129_AB.query('AB_s129 == "B"')[['chrom','start','end']])

#%%
fig,ax = plt.subplots(figsize=(6,8))
ind = numpy.arange(4)    # the x locations for the groups

n_comps = len(cast_s129_AB)
ax.bar(ind,[len(cast_A_s129_A_pb)/n_comps,len(cast_B_s129_B_pb)/n_comps,len(cast_A_s129_B_pb)/n_comps,len(cast_B_s129_A_pb)/n_comps],color=["green","red","blue","blue"])
print([len(cast_A_s129_A_pb)/n_comps,len(cast_B_s129_B_pb)/n_comps,len(cast_A_s129_B_pb)/n_comps,len(cast_B_s129_A_pb)/n_comps])
plt.xticks(ind, ('A/A', 'B/B', 'A/B', 'B/A'))
ax.set_xlabel("CAST/S129 compartments")
ax.set_ylabel("Percentage")
plt.title("{} dataset".format(np))
fig.savefig("{}Percentage_of_compartments_{}.png".format(save_path,np_),dpi=300,bbox_inches='tight')
#%%
# now we will overlap other features with AB

cast_ctcf = "{}CTCF.CAST.specific.peaks.bed".format(root)
s129_ctcf = "{}CTCF.S129.specific.peaks.bed".format(root)
common_ctcf = "{}CTCF.common.peaks.bed".format(root)
ctcf_data = [common_ctcf,cast_ctcf,s129_ctcf]
cast_atac = "{}F123_bing_data030420_genome_CAST.bed".format(root)
s129_atac = "{}F123_bing_data030420_genome_S129.bed".format(root)
common_atac = "{}common_peaks.bed".format(root)
atac_data = [common_atac,cast_atac,s129_atac]
cast_k4me3 = "{}H3K4me3.CAST.specific.peaks.bed".format(root)
s129_k4me3 = "{}H3K4me3.S129.specific.peaks.bed".format(root)
common_k4me3 = "{}H3K4me3.common.peaks.bed".format(root)
k4me3_data = [common_k4me3,cast_k4me3,s129_k4me3]
cast_k27ac = "{}H3K27Ac.CAST.specific.peaks.bed".format(root)
s129_k27ac = "{}H3K27Ac.S129.specific.peaks.bed".format(root)
common_k27ac = "{}H3K27Ac.common.peaks.bed".format(root)
k27ac_data = [common_k27ac,cast_k27ac,s129_k27ac]
s5 = "{}F123.S5p.BCP_peaks.HM_mode.bed".format(root)
s7 = "{}F123.S7p.BCP_peaks.HM_mode.bed".format(root)
k27me3 = "{}F123.K27me3.BCP_peaks.HM_mode.bed".format(root)
s5_data = [s5,s5,s5]
s7_data = [s7,s7,s7]
k27me3_data = [k27me3,k27me3,k27me3]

datasets = [ctcf_data,atac_data,k4me3_data,k27ac_data,s5_data,s7_data,k27me3_data]
colors = ["grey","skyblue","royalblue","palegreen","slateblue","steelblue","crimson"]

#allele specific features in compartments

fig,axs = plt.subplots(1,4,figsize=(8,4),sharey=True)
features = ["cast CTCF","cast ATAC","cast K4me3","cast K27ac"]
for n,feature in enumerate(features):
    cast_pb = pb.BedTool(datasets[n][1]).sort()

    cast_a_cov = cast_A_pb.coverage(cast_pb)
    cast_b_cov = cast_B_pb.coverage(cast_pb)
    cast_ab_cov = cast_A_s129_B_pb.coverage(cast_pb)
    cast_a_cov_df = cast_a_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    cast_b_cov_df = cast_b_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    cast_ab_cov_df = cast_ab_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    c_A = sum(i > 0 for i in cast_a_cov_df["n_peaks"])/len(cast_a_cov_df)
    c_B = sum(i > 0 for i in cast_b_cov_df["n_peaks"])/len(cast_b_cov_df)
    c_AB = sum(i > 0 for i in cast_ab_cov_df["n_peaks"])/len(cast_ab_cov_df)
    print ("Number of non 0 {} peaks in CAST A comps, B, and AB (ratio): {:.3f} {:.3f} {:.3f}".format(feature, c_A,c_B,c_AB))

    #if ratio of peaks in comps
    axs[n].bar(numpy.arange(3),[c_A,c_AB,c_B],color=colors[n])
    axs[0].set_ylabel("% of compartments with 1 peak or more")
    axs[n].set_title(features[n])
    axs[n].set_xticks([0,1,2])
    axs[n].set_xticklabels(["A","specific A","B"])
    axs[n].set_xlabel("CAST compartments")
    axs[n].set_ylim(0,0.15)


fig.savefig("{}CAST_compartments_n_peaks_CAST_features_{}.png".format(save_path,np_),dpi=300,bbox_inches='tight')

fig,axs = plt.subplots(1,4,figsize=(8,4),sharey=True)
features = ["s129 CTCF","s129 ATAC","s129 K4me3","s129 K27ac"]
for n,feature in enumerate(features):
    s129_pb = pb.BedTool(datasets[n][2]).sort()

    s129_a_cov = s129_A_pb.coverage(s129_pb)
    s129_b_cov = s129_B_pb.coverage(s129_pb)
    s129_ab_cov = cast_B_s129_A_pb.coverage(s129_pb)
    s129_a_cov_df = s129_a_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    s129_b_cov_df = s129_b_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    s129_ab_cov_df = s129_ab_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])

    c_A = sum(i > 0 for i in s129_a_cov_df["n_peaks"])/len(s129_a_cov_df)
    c_B = sum(i > 0 for i in s129_b_cov_df["n_peaks"])/len(s129_b_cov_df)
    c_AB = sum(i > 0 for i in s129_ab_cov_df["n_peaks"])/len(s129_ab_cov_df)
    print ("Number of non 0 {} peaks in S129 A comps, B, and AB (ratio): {:.3f} {:.3f} {:.3f} ".format(feature, c_A,c_B,c_AB))

    #if ratio of peaks in comps
    axs[n].bar(numpy.arange(3),[c_A,c_AB,c_B],color=colors[n])
    axs[0].set_ylabel("% of compartmenbts with 1 peak or more")
    axs[n].set_title(features[n])
    axs[n].set_xticks([0,1,2])
    axs[n].set_xticklabels(["A","specific A","B"])
    axs[n].set_xlabel("S129 compartments")
    axs[n].set_ylim(0,0.15)
fig.savefig("{}S129_compartments_n_peaks_S129_features_{}.png".format(save_path,np_),dpi=300,bbox_inches='tight')

#%%

fig,axs = plt.subplots(1,4,figsize=(18,12))
features = ["cast CTCF","cast ATAC","cast K4me3","cast K27ac"]
for n,feature in enumerate(features):
    cast_pb = pb.BedTool(datasets[n][1]).sort()

    cast_a_cov = cast_A_pb.coverage(cast_pb)
    cast_b_cov = cast_B_pb.coverage(cast_pb)
    cast_ab_cov = cast_A_s129_B_pb.coverage(cast_pb)
    cast_a_cov_df = cast_a_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    cast_b_cov_df = cast_b_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    cast_ab_cov_df = cast_ab_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])

    bp1 = axs[n].boxplot([cast_a_cov_df["coverage"].values.tolist(),cast_ab_cov_df["coverage"].values.tolist(),cast_b_cov_df["coverage"].values.tolist()],showfliers=True,patch_artist=True,labels=["A","specific A","B"])
    plt.setp(bp1["boxes"],color=colors[n])
    axs[0].set_ylabel("coverage")
    axs[n].set_title(features[n])
    axs[n].set_xlabel("CAST compartments {}".format(np))
fig.savefig("{}CAST_compartments_coverage_CAST_features_{}.png".format(save_path,np_),dpi=300,bbox_inches='tight')

fig,axs = plt.subplots(1,4,figsize=(18,12))
features = ["s129 CTCF","s129 ATAC","s129 K4me3","s129 K27ac"]
for n,feature in enumerate(features):
    s129_pb = pb.BedTool(datasets[n][2]).sort()

    s129_a_cov = s129_A_pb.coverage(s129_pb)
    s129_b_cov = s129_B_pb.coverage(s129_pb)
    s129_ab_cov = cast_B_s129_A_pb.coverage(s129_pb)
    s129_a_cov_df = s129_a_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    s129_b_cov_df = s129_b_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    s129_ab_cov_df = s129_ab_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])

    bp1 = axs[n].boxplot([s129_a_cov_df["coverage"].values.tolist(),s129_ab_cov_df["coverage"].values.tolist(),s129_b_cov_df["coverage"].values.tolist()],showfliers=True,patch_artist=True,labels=["A","specific A","B"])
    plt.setp(bp1["boxes"],color=colors[n])
    axs[0].set_ylabel("coverage")
    axs[n].set_title(features[n])
    axs[n].set_xlabel("S129 compartments {}".format(np))
fig.savefig("{}S129_compartments_coverage_S129_features_{}.png".format(save_path,np_),dpi=300,bbox_inches='tight')

fig,axs = plt.subplots(1,7,figsize=(12,5),sharey=False)
features = ["common CTCF","common ATAC","common K4me3","common K27ac", "S5p", "S7p", "K27m3"]
for n,feature in enumerate(features):
    common_pb = pb.BedTool(datasets[n][0]).sort()

    cast_a_cov = cast_A_pb.coverage(common_pb)
    cast_b_cov = cast_B_pb.coverage(common_pb)
    cast_ab_cov = cast_A_s129_B_pb.coverage(common_pb)
    cast_a_cov_df = cast_a_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    cast_b_cov_df = cast_b_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    cast_ab_cov_df = cast_ab_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    s129_a_cov = s129_A_pb.coverage(common_pb)
    s129_b_cov = s129_B_pb.coverage(common_pb)
    s129_ab_cov = cast_B_s129_A_pb.coverage(common_pb)
    s129_a_cov_df = s129_a_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    s129_b_cov_df = s129_b_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])
    s129_ab_cov_df = s129_ab_cov.to_dataframe(names=["chrom","start","end","n_peaks","covered","window","coverage"])

    bp1 = axs[n].boxplot([cast_a_cov_df["coverage"].values.tolist(),cast_ab_cov_df["coverage"].values.tolist(),cast_b_cov_df["coverage"].values.tolist(),s129_a_cov_df["coverage"].values.tolist(),s129_ab_cov_df["coverage"].values.tolist(),s129_b_cov_df["coverage"].values.tolist()],showfliers=False,patch_artist=True,labels=["CAST A","CAST specific A","CAST B","S129 A", "S129 specific A", "S129 B"])
    plt.setp(bp1["boxes"],color=colors[n])
    axs[0].set_ylabel("coverage")
    axs[n].set_title(features[n])
    axs[n].set_xlabel("CAST compartments {}".format(np))
    #axs[n].set_ylim([0,0.8])
fig.savefig("{}CAST_S129_compartments_coverage_common_features_{}.pdf".format(save_path,np_),dpi=300,bbox_inches='tight')


# %%