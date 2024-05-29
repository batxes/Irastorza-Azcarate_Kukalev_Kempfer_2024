#%%
import seaborn as sns
import pybedtools as pb
import operator
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *


def get_WDF_changes():
    def flatten(l):
        flattenned = []
        for x in l:
            for y in x:
                flattenned.append(y)
        return flattenned


    path_cast = root+"Curated_F123_1122_ALL_as_3NP_CAST_at50000.passed_qc_fc5_cw6_s11.table"
    path_s129 = root+"Curated_F123_1122_ALL_as_3NP_S129_at50000.passed_qc_fc5_cw6_s11.table"

    seg_table = pd.read_csv(path_cast,sep="\t")
    seg_table.set_index(['chrom','start','stop'], inplace = True)
    seg_table_cast = seg_table.loc[chrs_]
    seg_table = pd.read_csv(path_s129,sep="\t")
    seg_table.set_index(['chrom','start','stop'], inplace = True)
    seg_table_s129 = seg_table.loc[chrs_]
    n_NPs = seg_table_cast.shape[1]

    def normalize_wdf(l):
        aux_l = []
        max_v = max(l)
        for v in l:
            aux_l.append(v/max_v)
        return np.array(aux_l)
    
    def subtract_and_remove_zeros(l1,l2,starts,ends):
        aux_l = []
        wdf1 = []
        wdf2 = []
        s_l = []
        e_l = []
        for v1,v2,s,e in zip(l1,l2,starts,ends):
            if v1 == 0 or v2 == 0:
                pass
            else:
                aux_l.append(v1-v2)
                wdf1.append(v1)
                wdf2.append(v2)
                s_l.append(s)
                e_l.append(e)
        return np.array(aux_l),s_l,e_l,wdf1,wdf2

    wdf_diff_values = []
    wdf1_list = []
    wdf2_list = []
    chroms = []
    starts = []
    ends = []
    for chr_ in chrs_:
        df = seg_table_cast.loc[chr_]
        df2 = seg_table_s129.loc[chr_]
        aux_df_sum = df.sum(axis=1)
        aux_df_sum2 = df2.sum(axis=1)
        wdf_values = aux_df_sum.values/n_NPs
        wdf_values2 = aux_df_sum2.values/n_NPs
        index = list(df.index.values)
        s_l = []
        e_l = []
        for x in index:
            s_l.append(x[0])
            e_l.append(x[1])
        wdf_diff_v,s,e,wdf1,wdf2 = subtract_and_remove_zeros(wdf_values,wdf_values2,s_l,e_l)
        wdf1_list.append(wdf1)
        wdf2_list.append(wdf2)
        starts.extend(s)
        ends.extend(e)
        wdf_diff_values.append(wdf_diff_v)
        chroms.append([chr_]*len(wdf_diff_v))
    wdf_diff_values = flatten(wdf_diff_values)
    chroms = flatten(chroms)
    out_df = pd.DataFrame({"chrom":chroms,"start":starts,"end":ends,"wdf_diff":wdf_diff_values})
    out_df['wdf_diff'] = out_df['wdf_diff'].fillna(0)

    #plot wdf pero allele and chr
    fig,ax = plt.subplots(figsize=(12,12))
    wdf1_list = flatten(wdf1_list)
    allele = ["cast"]*len(wdf1_list)
    out_df1 = pd.DataFrame({"chrom":chroms,"start":starts,"end":ends,"wdf":wdf1_list,"allele":allele})
    wdf2_list = flatten(wdf2_list)
    allele = ["s129"]*len(wdf2_list)
    out_df2 = pd.DataFrame({"chrom":chroms,"start":starts,"end":ends,"wdf":wdf2_list,"allele":allele})
    out_df1['wdf'] = out_df1['wdf'].fillna(0)
    out_df2['wdf'] = out_df2['wdf'].fillna(0)
    out_df_both = pd.concat([out_df1,out_df2])
    out_df1 = pd.DataFrame({"chrom":chroms,"start":starts,"end":ends,"wdf":wdf1_list})
    out_df2 = pd.DataFrame({"chrom":chroms,"start":starts,"end":ends,"wdf":wdf2_list})
    ax = sns.violinplot(x="chrom",y="wdf",hue="allele",ax=ax,data=out_df_both,palette="PuOr_r",split=True,figsize=(18,12),inner="quartile")

    fig, ax = plt.subplots(figsize = (12,12))
    ax = sns.barplot(x="chrom", y="wdf", hue="allele", data=out_df_both,palette="PuOr_r",estimator=np.median)

    return out_df,out_df1,out_df2

path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"
df = pd.read_csv(path,sep="\t")
ase_df = df[df["chrom"].isin(chrs_)]
cols = ase_df.columns 
ase_pb = pb.BedTool.from_dataframe(ase_df)

cast_genes = ase_df.query('log2foldchange >= 1 and p_adj_sum < 0.05 and TPM_transcript >= 1')[["chrom","TSS_start","TSS_end"]]
s129_genes = ase_df.query('log2foldchange <= -1 and p_adj_sum < 0.05 and TPM_transcript >= 1')[["chrom","TSS_start","TSS_end"]]
exp_genes_w_SNP = ase_df.query('TPM_transcript >= 1 and number_SNPs > 0')[["chrom","TSS_start","TSS_end"]]
cast_pb = pb.BedTool.from_dataframe(cast_genes)
s129_pb = pb.BedTool.from_dataframe(s129_genes)
exp_w_SNP_pb = pb.BedTool.from_dataframe(exp_genes_w_SNP)

diff_wdf_df,wdf_cast, wdf_s129 = get_WDF_changes()
diff_wdf_pb = pb.BedTool.from_dataframe(diff_wdf_df)
cast_wdf_pb = pb.BedTool.from_dataframe(wdf_cast)
s129_wdf_pb = pb.BedTool.from_dataframe(wdf_s129)
# %%
names=["chrom_bin","bin_start","bin_end","diff_wdf"] 
names.extend(cols)
ase_wdf_df = diff_wdf_pb.intersect(ase_pb,F=0.5,wa=True,wb=True).to_dataframe(names=names)
ase_wdf_df = ase_wdf_df.replace(".",0)
ase_wdf_df.to_csv(save_path+"test.tsv",sep="\t")
ase_wdf_df['log2foldchange'] = ase_wdf_df['log2foldchange'].astype(float)
ase_wdf_df['p_adj_sum'] = ase_wdf_df['p_adj_sum'].astype(float)
aux_ase = ase_wdf_df.copy()
aux_ase = aux_ase.rename(columns={"bin_start":"start","bin_end":"end"})

master_df = wdf_cast.copy()
master_pb = pb.BedTool.from_dataframe(master_df)

master_df = master_df.rename(columns={"wdf":"wdf_cast"})
master_df = master_df.merge(wdf_s129,on=["chrom","start","end"])
master_df = master_df.rename(columns={"wdf":"wdf_s129"})
master_df = master_df.merge(diff_wdf_df,on=["chrom","start","end"])
master_with_genes_df = master_df.merge(aux_ase,on=["chrom","start","end"],how="outer")

#%%

# WDF vs TPM

fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type != 'no' and type != 'exp_wo_SNP' and type != 'exp_w_SNP'")

sns.kdeplot(
    data=df_to_plot, x="log2foldchange", y="wdf_diff", hue="type", fill=True,levels=10,
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")

fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type == 'exp_w_SNP'")

sns.kdeplot(
    data=df_to_plot, x="log2foldchange", y="wdf_diff", color="green", fill=True,levels=10,
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
#%%
fig, ax = plt.subplots(figsize=(8,8))
df_to_plot = master_with_genes_df.query("type == 'cast' and type != 's129' and Promoter_state != 'Overlapping_TSSs' and Promoter_state != 'Internal_gene'")
sns.kdeplot(
    data=df_to_plot, x="log2foldchange", y="wdf_diff", color=cast_color, fill=False,levels=15,
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
df_to_plot = master_with_genes_df.query("type == 's129'")
sns.kdeplot(
    data=df_to_plot, x="log2foldchange", y="wdf_diff", color=s129_color, fill=False,levels=15,
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
fig.savefig(save_path+"diff_wdf_vs_log2.pdf",dpi=300)

fig, (ax2,ax1) = plt.subplots(1,2,figsize=(12,8))
df_to_plot = master_with_genes_df.query("type != 'cast' and type == 's129' and Promoter_state != 'Overlapping_TSSs' and Promoter_state != 'Internal_gene'")
sns.kdeplot(
    data=df_to_plot, x="log2foldchange", y="wdf_s129", color=s129_color, fill=False,levels=15,ax=ax2
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
ax2.set_ylim(0,0.15)
df_to_plot = master_with_genes_df.query("type == 'cast'")
sns.kdeplot(
    data=df_to_plot, x="log2foldchange", y="wdf_cast", color=cast_color, fill=False,levels=15,ax=ax1
)
ax1.set_ylim(0,0.15)
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
fig.tight_layout()
fig.savefig(save_path+"wdf_vs_log2_density.pdf",dpi=300)


fig, (ax2,ax1) = plt.subplots(1,2,figsize=(12,8))
#df_to_plot = master_with_genes_df.query("type == 's129'")
df_to_plot = master_with_genes_df.query("type != 'cast' and type == 's129' and Promoter_state != 'Overlapping_TSSs' and Promoter_state != 'Internal_gene'")
sns.scatterplot(
    data=df_to_plot, x="log2foldchange", y="wdf_s129", color=s129_color,ax=ax2
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
ax2.set_ylim(0,0.15)
df_to_plot = master_with_genes_df.query("type == 'cast'")
sns.scatterplot(
    data=df_to_plot, x="log2foldchange", y="wdf_cast", color=cast_color,ax=ax1
)
ax1.set_ylim(0,0.15)
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
fig.tight_layout()
fig.savefig(save_path+"wdf_vs_log2.pdf",dpi=300)

#%%
#get reads per bin for each allele
# first, remove overlapping genes
df_to_plot = master_with_genes_df.query("type == 'cast' or type == 's129' or type == 'exp_w_SNP'")
df_to_plot["CASTcount_sum"] = df_to_plot["CASTcount_sum"].astype(float)
df_to_plot["J129count_sum"] = df_to_plot["J129count_sum"].astype(float)
master_with_non_overlapping_genes_df = df_to_plot.query("Promoter_state != 'Overlapping_TSSs' and Promoter_state != 'Internal_gene'")
reads_per_bin = master_with_non_overlapping_genes_df.groupby(by=["chrom","start","end","wdf_cast","wdf_s129"]).sum()
reads_per_bin = reads_per_bin.reset_index()
print (reads_per_bin[["chrom","start","end","CASTcount_sum","J129count_sum","wdf_cast","wdf_s129"]])
fig, ax = plt.subplots(figsize=(12,12))

ax.scatter(reads_per_bin.CASTcount_sum,reads_per_bin.wdf_cast, color=cast_color,alpha=0.2)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
#ax.set_xlim(-1,30)
ax.set_xlim(-1000,20000)
ax.set_xlabel("Number of CAST RNA-seq reads")
ax.set_ylabel("CAST WDF")
ax.set_title("CAST S129 and exp_w_SNP genes with non overlapping genes")
fig.savefig(save_path+"wdf_vs_cast_reads.pdf",dpi=300)
fig, ax = plt.subplots(figsize=(12,12))
ax.scatter(reads_per_bin.J129count_sum,reads_per_bin.wdf_s129, color=s129_color,alpha=0.2)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
ax.set_xlim(-1000,20000)
ax.set_xlabel("Number of S129 RNA-seq reads")
ax.set_ylabel("S129 WDF")
ax.set_title("CAST S129 and exp_w_SNP genes with non overlapping genes")
fig.savefig(save_path+"wdf_vs_s129_reads.pdf",dpi=300)


#%%
fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type == 'cast' and type != 's129' and Promoter_state != 'Overlapping_TSSs' and Promoter_state != 'Internal_gene'")
df_to_plot["CASTcount_sum"] = df_to_plot["CASTcount_sum"].astype(float)
df_to_plot["J129count_sum"] = df_to_plot["J129count_sum"].astype(float)
ax.scatter(df_to_plot.CASTcount_sum,df_to_plot.wdf_cast, color=cast_color,alpha=0.3)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
ax.set_xlabel("Number of CAST RNA-seq reads")
ax.set_ylabel("CAST WDF")
ax.set_title("CAST genes with non overlapping genes")
ax.set_xlim(-100,5000)
fig.savefig(save_path+"wdf_vs_cast_reads_cast_genes.pdf",dpi=300)
fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type == 's129' and type != 'cast' and Promoter_state != 'Overlapping_TSSs' and Promoter_state != 'Internal_gene'")
df_to_plot["CASTcount_sum"] = df_to_plot["CASTcount_sum"].astype(float)
df_to_plot["J129count_sum"] = df_to_plot["J129count_sum"].astype(float)
ax.scatter(df_to_plot.J129count_sum,df_to_plot.wdf_s129, color=s129_color,alpha=0.3)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
ax.set_xlabel("Number of S129 RNA-seq reads")
ax.set_ylabel("S129 WDF")
ax.set_title("S129 genes with non overlapping genes")
ax.set_xlim(-100,5000)
#ax.set_xlim(-1,30)
#ax.set_xlim(-10,300)
fig.savefig(save_path+"wdf_vs_s129_reads_s129_genes.pdf",dpi=300)

#%%
fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type == 'exp_w_SNP'")
sns.kdeplot(
    data=df_to_plot, x="TPM_transcript", y="wdf_cast", color="green", fill=True,levels=20,
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
ax.set_xlim(-100,500)
fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type == 'exp_w_SNP'")
sns.kdeplot(
    data=df_to_plot, x="TPM_transcript", y="wdf_s129", color="green", fill=True,levels=20,
)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
ax.set_xlim(-100,500)
#%%
fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type == 'exp_w_SNP'")
ax.scatter(df_to_plot.TPM_transcript,df_to_plot.wdf_cast)
ax.axhline(0,color="black")
ax.axvline(0,color="black")
fig, ax = plt.subplots(figsize=(12,12))
df_to_plot = master_with_genes_df.query("type == 'exp_w_SNP'")
ax.scatter(df_to_plot.TPM_transcript,df_to_plot.wdf_s129)
ax.axhline(0,color="black")
ax.axvline(0,color="black")


## %%
