#%%

#this master table generator is meant to:

# First, populate the ASE table with peaks and more info
# second, use this table to check Promoters and see different classes. Upset plots would help
# third, the same as second but starting from PRC2 peaks
# All of this is due to the fact that I see very interesting combinations of prc2 with CpG islands and more allelic k27me3 reads in Allelic genes (of the same type)

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
from upsetplot import UpSet

ase_path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype.bed"

recalculate = False
if recalculate:
    ase_df = pd.read_csv(ase_path,sep="\t")
    ase_df = ase_df.drop_duplicates()
    ase_df = ase_df[ase_df["chrom"].isin(chrs_)]

    path_dnamet = "{}DNA_met_candidates.bed".format(root) # from DNAmet.py
    dnamet_df = pd.read_csv(path_dnamet,sep="\t")
    merge_on = ase_df.columns.values.tolist()

    ase_df = pd.merge(ase_df,dnamet_df,on=merge_on,how="left")
    master_pb = pb.BedTool.from_dataframe(ase_df)

    cast_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_bing_data030420_genome_CAST.bed"
    s129_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_bing_data030420_genome_S129.bed"
    common_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_ATAC_Bing_shared_peaks.bed"
    nonphased_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_ATAC_Bing_non_phased_peaks.bed"
    all_atac = "/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/common_peaks.bed"

    cast_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/CTCF.CAST.specific.peaks.bed"
    s129_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/CTCF.S129.specific.peaks.bed"
    common_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/CTCF.common.peaks.bed"
    nonphased_ctcf = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/CTCF.non_phased.peaks.nosex.bed"

    cast_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/Rad21.CAST.specific.peaks.bed"
    s129_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/Rad21.S129.specific.peaks.bed"
    common_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/Rad21.common.peaks.bed"
    nonphased_rad21 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/Rad21.common.nonphased.peaks.nosex.bed"

    cast_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K4me3.CAST.specific.peaks.bed"
    s129_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K4me3.S129.specific.peaks.bed"
    common_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K4me3.common.peaks.bed"
    nonphased_h3k4me3 = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/H3K4me3.non_phased.peaks.nosex.bed"

    cast_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K27Ac.CAST.specific.peaks.bed"
    s129_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K27Ac.S129.specific.peaks.bed"
    common_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks/H3K27Ac.common.peaks.bed"
    nonphased_h3k27ac = "/home/ibai/pombo_johnny5/F123/tracks/ChIPseq_allele_specific_peaks_ibai/H3K27Ac.non_phased.peaks.nosex.bed"

    s5p = "/home/ibai/pombo_johnny5/F123/tracks/S5p_S7p_K27/F123.S5p.BCP_peaks.HM_mode.bed"
    s7p = "/home/ibai/pombo_johnny5/F123/tracks/S5p_S7p_K27/F123.S7p.BCP_peaks.HM_mode.bed"
    H3K27me3 = "/home/ibai/pombo_johnny5/F123/tracks/S5p_S7p_K27/F123.K27me3.only_peaks.bed"

    path = "/home/ibai/pombo_johnny5/Sasha/Projects/F123/DNAMethylation/"
    cast_met1 = path+"Bismark_mapping_R1/Final_processed/Coverage/F123_WGBS_Repl1.bismark.deduplicated.genome1.bismark.cov"
    cast_met2 = path+"Bismark_mapping_R2/Final_processed/Coverage/F123_WGBS_Repl2.bismark.deduplicated.genome1.bismark.cov"
    s129_met1 = path+"Bismark_mapping_R1/Final_processed/Coverage/F123_WGBS_Repl1.bismark.deduplicated.genome2.bismark.cov"
    s129_met2 = path+"Bismark_mapping_R2/Final_processed/Coverage/F123_WGBS_Repl2.bismark.deduplicated.genome2.bismark.cov"

    cpg_path = "/home/ibai/pombo_johnny5/F123/tracks/CpG_mm10.table"
    cpg_df = pd.read_csv(cpg_path, sep="\t")
    cpg_df = cpg_df[cpg_df["chrom"].isin(chrs_)]
    cpg_pb = pb.BedTool.from_dataframe(cpg_df[["chrom","start","end"]]).sort()

    cast_atac_pb = pb.BedTool(cast_atac).sort()
    s129_atac_pb = pb.BedTool(s129_atac).sort()
    common_atac_pb = pb.BedTool(common_atac).sort()
    nonphased_atac_pb = pb.BedTool(nonphased_atac).sort()
    cast_ctcf_pb = pb.BedTool(cast_ctcf).sort()
    s129_ctcf_pb = pb.BedTool(s129_ctcf).sort()
    common_ctcf_pb = pb.BedTool(common_ctcf).sort()
    nonphased_ctcf_pb = pb.BedTool(nonphased_ctcf).sort()
    cast_rad21_pb = pb.BedTool(cast_rad21).sort()
    s129_rad21_pb = pb.BedTool(s129_rad21).sort()
    common_rad21_pb = pb.BedTool(common_rad21).sort()
    nonphased_rad21_pb = pb.BedTool(nonphased_rad21).sort()
    cast_h3k4me3_pb = pb.BedTool(cast_h3k4me3).sort()
    s129_h3k4me3_pb = pb.BedTool(s129_h3k4me3).sort()
    common_h3k4me3_pb = pb.BedTool(common_h3k4me3).sort()
    nonphased_h3k4me3_pb = pb.BedTool(nonphased_h3k4me3).sort()
    cast_h3k27ac_pb = pb.BedTool(cast_h3k27ac).sort()
    s129_h3k27ac_pb = pb.BedTool(s129_h3k27ac).sort()
    common_h3k27ac_pb = pb.BedTool(common_h3k27ac).sort()
    nonphased_h3k27ac_pb = pb.BedTool(nonphased_h3k27ac).sort()
    s5_pb = pb.BedTool(s5p).sort()
    s7_pb = pb.BedTool(s7p).sort()
    h3K27me3_pb = pb.BedTool(H3K27me3).sort()
    print ("ALL files read")


    master_pb = master_pb.intersect(cast_atac_pb,c=True)
    master_pb = master_pb.intersect(s129_atac_pb,c=True)
    master_pb = master_pb.intersect(common_atac_pb,c=True)
    master_pb = master_pb.intersect(nonphased_atac_pb,c=True)
    master_pb = master_pb.intersect(cast_ctcf_pb,c=True)
    master_pb = master_pb.intersect(s129_ctcf_pb,c=True)
    master_pb = master_pb.intersect(common_ctcf_pb,c=True)
    master_pb = master_pb.intersect(nonphased_ctcf_pb,c=True)
    master_pb = master_pb.intersect(cast_rad21_pb,c=True)
    master_pb = master_pb.intersect(s129_rad21_pb,c=True)
    master_pb = master_pb.intersect(common_rad21_pb,c=True)
    master_pb = master_pb.intersect(nonphased_rad21_pb,c=True)
    master_pb = master_pb.intersect(cast_h3k4me3_pb,c=True)
    master_pb = master_pb.intersect(s129_h3k4me3_pb,c=True)
    master_pb = master_pb.intersect(common_h3k4me3_pb,c=True)
    master_pb = master_pb.intersect(nonphased_h3k4me3_pb,c=True)
    master_pb = master_pb.intersect(cast_h3k27ac_pb,c=True)
    master_pb = master_pb.intersect(s129_h3k27ac_pb,c=True)
    master_pb = master_pb.intersect(common_h3k27ac_pb,c=True)
    master_pb = master_pb.intersect(nonphased_h3k27ac_pb,c=True)
    master_pb = master_pb.intersect(s5_pb,c=True)
    master_pb = master_pb.intersect(s7_pb,c=True)
    master_pb = master_pb.intersect(h3K27me3_pb,c=True)
    master_pb = master_pb.intersect(cpg_pb,c=True)

    columns = ase_df.columns.values.tolist()
    columns2 = ["n_cast_ATAC","n_s129_ATAC","n_common_ATAC","n_nonphased_ATAC","n_cast_CTCF","n_s129_CTCF","n_common_CTCF","n_nonphased_CTCF","n_cast_Rad21","n_s129_Rad21","n_common_Rad21","n_nonphased_Rad21",
                    "n_cast_h3k4me3","n_s129_h3k4me3","n_common_h3k4me3","n_nonphased_h3k4me3","n_cast_h3k27ac","n_s129_h3k27ac","n_common_h3k27ac","n_nonphased_h3k27ac",
                    "S5p","S7p","H3K27me3","n_cpg"]


    master_df = master_pb.to_dataframe(names=columns+columns2)
    print (master_df)
    master_df = master_df.replace(".",np.nan)

    master_df.log2foldchange = master_df.log2foldchange.astype(float)
    master_df.p_adj_sum = master_df.p_adj_sum.astype(float)
    master_df.by_DNAmet = master_df.by_DNAmet.astype(float)

    master_df['type'] = 'no'
    master_df.loc[(master_df['TPM_transcript'] >= 1) & (master_df['log2foldchange'] <= -1) & (master_df['p_adj_sum'] <= 0.05),'type'] = 's129'
    master_df.loc[(master_df['TPM_transcript'] >= 1) & (master_df['log2foldchange'] >=  1) & (master_df['p_adj_sum'] <= 0.05),'type'] = 'cast'
    master_df.loc[(master_df['TPM_transcript'] >= 1) & (((master_df['log2foldchange'] <  1) & (master_df['log2foldchange'] > -1)) | (master_df['p_adj_sum'] > 0.05)),'type'] = 'exp_w_SNP'
    master_df.loc[(master_df['TPM_transcript'] >= 1) & (master_df['number_SNPs'] == 0),'type'] = 'exp_wo_SNP'
    master_df["by_DNAmet"] = master_df["by_DNAmet"].replace(np.nan,0)
    master_df.to_csv(root+"/300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed",sep="\t",index=False)
else:
    master_df = pd.read_csv(root+"/300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed",sep="\t",index_col=None)

def plot_boxplot(values,labels,colors,ylabel,save_path,scatter=True):
    fig,ax = plt.subplots()
    bp1 = ax.boxplot(values,patch_artist=True ,showfliers=False,widths=0.5)
    ax.set_xticklabels(labels)
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set(alpha=0.5)

    if scatter:
        for pos,v in enumerate(values): 
            xv = np.random.normal(pos+1, 0.05, len(v))
            ax.scatter(xv, v,s=1,alpha=1,color="black")
    plt.setp(bp1["medians"],color="white")
    ax.set_ylabel(ylabel)
    plt.savefig(save_path,dpi=300)


def merge_cols (df,l,new_col):
    subset_df = df[l]
    df[new_col] = subset_df.sum(axis=1)
    df = df.drop(columns=l)
    cols = df.columns.values.tolist()
    cols.remove("type")
    df[cols] = df[cols].apply(pd.to_numeric)
    a = np.array(df[new_col].values.tolist())
    df[new_col] = np.where(a > 1, 1, a).tolist()
    df = df[cols+["type"]]
    return df

#%%

# Filter information depending on which genes we want
# This plot is not used in final manuscript

#subset_df = master_df.query('TPM_transcript >= 1 and number_SNPs > 0 and log2foldchange > -1 and log2foldchange < 1')
#label = "exp_not_ASE"
subset_df = master_df.query('TPM_transcript >= 1 and (log2foldchange <= -1 or log2foldchange >= 1) and p_adj_sum <= 0.05')
label = "ASE"
#subset_df = master_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05') #S129
#label = "s129"
#subset_df = master_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05') #CAST
#label="cast"
subset_df = subset_df.query("Promoter_state != 'Internal_gene' and Promoter_state != 'Overlapping_TSSs'")

upset_df = subset_df[["n_cast_ATAC","n_s129_ATAC","n_common_ATAC",
                    "n_cast_h3k4me3","n_s129_h3k4me3","n_common_h3k4me3",
                    "n_cast_h3k27ac","n_s129_h3k27ac","n_common_h3k27ac",
                    "S5p","S7p","H3K27me3","type"]]
categories_df = upset_df[upset_df.columns[:-1]]
print (categories_df)

binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))
facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=5)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
upset.style_subsets(present=["H3K27me3_", "n_common_h3k4me3_"],
                    facecolor="#4eb3d3",
                    label="bivalent")
upset.style_subsets(present=["H3K27me3_", "n_cast_h3k4me3_"],
                    facecolor="#4eb3d3",
                    label="bivalent")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=9)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=9)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}all_features_upset_plot_{}.pdf".format(save_path,label),dpi=300)


#%%

# If we add nonphased h3k4me4
# This plot is not used in final manuscript
subset_df["H3K4me3"] = subset_df.n_cast_h3k4me3 +subset_df.n_s129_h3k4me3 +subset_df.n_common_h3k4me3 +subset_df.n_nonphased_h3k4me3 

def xor(x):
    if x >= 1:
        return 1
    else:
        return 0

subset_df['H3K4me3'] = subset_df['H3K4me3'].apply(xor)
upset_df = subset_df[["S7p","S5p","H3K4me3","type","TPM_transcript"]]
categories_df = upset_df[upset_df.columns[:-2]]

binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))
facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=5)

if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=3)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=3)

upset.add_catplot(value='TPM_transcript', kind='box',color="lightgrey",elements=9,showfliers=False)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
ax = plt.gca()
ax.set_ylim(0,50)

plt.title(label)
fig.savefig("{}upset_plot_k4me3_and_s7{}.pdf".format(save_path,label),dpi=300)


# If we add nonphased h3k27m3
# Figure 2e and 2f

upset_df = subset_df[["S7p","S5p","H3K27me3","type","TPM_transcript"]]
categories_df = upset_df[upset_df.columns[:-2]]

binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))
facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=5)

if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=3)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=3)

upset.add_catplot(value='TPM_transcript', kind='box',color="lightgrey",elements=9,showfliers=False)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
ax = plt.gca()
ax.set_ylim(0,50)

plt.title(label)
fig.savefig("{}upset_plot_k27me3_and_s7{}.pdf".format(save_path,label),dpi=300)

# Expression boxplots for S5 S7 and K27me3
plot_boxplot([  subset_df.query('S7p == 1 and H3K27me3 == 0 and S5p == 0 and (type == "cast" or type == "s129")').TPM_transcript,
                subset_df.query('S7p == 1 and H3K27me3 == 1 and S5p == 0 and (type == "cast" or type == "s129")').TPM_transcript,
                subset_df.query('S7p == 0 and H3K27me3 == 0 and S5p == 0 and (type == "cast" or type == "s129")').TPM_transcript,
                subset_df.query('S7p == 1 and H3K27me3 == 0 and S5p == 1 and (type == "cast" or type == "s129")').TPM_transcript,
                subset_df.query('S7p == 0 and H3K27me3 == 1 and S5p == 0 and (type == "cast" or type == "s129")').TPM_transcript,
                subset_df.query('S7p == 0 and H3K27me3 == 0 and S5p == 1 and (type == "cast" or type == "s129")').TPM_transcript,
                subset_df.query('S7p == 1 and H3K27me3 == 1 and S5p == 1 and (type == "cast" or type == "s129")').TPM_transcript,
                subset_df.query('S7p == 0 and H3K27me3 == 1 and S5p == 1 and (type == "cast" or type == "s129")').TPM_transcript],
                ["lol"]*8,
                ["grey"]*8,
                "Gene expression",
                "/home/ibai/test.pdf",
                scatter=False
              )



#%%

# These are other plots not used in final manuscript

upset_df = subset_df[["n_cast_ATAC","n_s129_ATAC",
                    "n_cast_h3k4me3","n_s129_h3k4me3",
                    "n_cast_h3k27ac","n_s129_h3k27ac",
                    "H3K27me3","type"]]
categories_df = upset_df[upset_df.columns[:-1]]
binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))
facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=1)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=9)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=9)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}upset_plot_{}_only_specific.pdf".format(save_path,label),dpi=300)


###################################################################################
upset_df = subset_df[["n_cast_ATAC","n_s129_ATAC",
                    "n_cast_h3k4me3","n_s129_h3k4me3",
                    "n_cast_h3k27ac","n_s129_h3k27ac",
                    "n_nonphased_ATAC","n_nonphased_h3k4me3","n_nonphased_h3k27ac",
                    "H3K27me3","type"]]
categories_df = upset_df[upset_df.columns[:-1]]
binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))

facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=1)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=9)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=9)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}upset_plot_{}_only_specific_and_nonphased.pdf".format(save_path,label),dpi=300)
###
###################################################################################
upset_df = subset_df[["n_common_ATAC",
                    "n_common_h3k4me3",
                    "n_common_h3k27ac",
                    "H3K27me3","type"]]
categories_df = upset_df[upset_df.columns[:-1]]
binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))

facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=1)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
upset.style_subsets(present=["H3K27me3_", "n_common_h3k4me3_"],
                    facecolor="#4eb3d3",
                    label="bivalent")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=3)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=3)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}upset_plot_{}_only_common.pdf".format(save_path,label),dpi=300)

###################################################################################
upset_df = subset_df[["n_cast_ATAC",
                    "n_cast_h3k4me3",
                    "n_cast_h3k27ac",
                    "H3K27me3","S7p","by_DNAmet","type"]]
upset_df = merge_cols(upset_df,["n_cast_ATAC","n_cast_h3k4me3","n_cast_h3k27ac"],"cast_k4m3_k27ac_atac")
categories_df = upset_df[upset_df.columns[:-1]]
binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))

facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=1)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=3)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=3)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}upset_plot_{}_only_cast.pdf".format(save_path,label),dpi=300)

###################################################################################

upset_df = subset_df[["n_s129_ATAC",
                    "n_s129_h3k4me3",
                    "n_s129_h3k27ac",
                    "H3K27me3","S7p","by_DNAmet","type"]]
print (upset_df)
upset_df = merge_cols(upset_df,["n_s129_ATAC","n_s129_h3k4me3","n_s129_h3k27ac"],"s129_k4m3_k27ac_atac")
print (upset_df)

categories_df = upset_df[upset_df.columns[:-1]]
binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))

facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=1)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=3)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=3)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}upset_plot_{}_only_s129.pdf".format(save_path,label),dpi=300)

###################################################################################

#now only with K4m3
upset_df = subset_df[["n_cast_h3k4me3","H3K27me3","S7p","type"]]
categories_df = upset_df[upset_df.columns[:-1]]
binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))

facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=1)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=3)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=3)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}upset_plot_{}_only_cast_and_k4me3.pdf".format(save_path,label),dpi=300)


###################################################################################


upset_df = subset_df[["n_s129_h3k4me3","H3K27me3","S7p","type"]]

categories_df = upset_df[upset_df.columns[:-1]]
binary_categories_df = categories_df > 0 #make categories True or False
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '_')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))

facecolor="black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True, min_subset_size=1)
upset.style_subsets(present=["H3K27me3_"],
                    facecolor="red",
                    label="PRC2")
if label == "cast":
    upset.add_stacked_bars(by="type", colors=[cast_color],
                        title="Number of ASE genes", elements=3)
else:
    upset.add_stacked_bars(by="type", colors=[cast_color,s129_color],
                        title="Number of ASE genes", elements=3)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(label)
fig.savefig("{}upset_plot_{}_only_s129_and_k4me3.pdf".format(save_path,label),dpi=300)
