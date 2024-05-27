#%%

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

non_phased_rep12 = root+"F123.WGBS.Merged.R1R2.deduplicated.bismark.cov.gz"
cast_rep12 = root+"F123.WGBS.Merged.R1R2.deduplicated.genome1.bismark.cov.gz"
s129_rep12 = root+"F123.WGBS.Merged.R1R2.deduplicated.genome2.bismark.cov.gz"

path_candidates = root+"DNA_met_candidates.bed"
cast_fasta = root+"CAST.mm10.fa" 
s129_fasta = root+"S129_SvJae.mm10.fa" 

#set promoters

calculate = False

#I set promoter as Sasha did, TSS -1000 +1000
promoter_up = 200 #100
promoter_down = 400 #200
promoter = [promoter_up,promoter_down]
dnamet_cutoff = 15 #Percentile
dnamet_cutoff = 5 #Percentile
cast_df = pd.read_csv(cast_rep12,sep="\t",compression="gzip",names=["chrom","start","end","per","met","unmet"])
s129_df = pd.read_csv(s129_rep12,sep="\t",compression="gzip",names=["chrom","start","end","per","met","unmet"])
cast_df = cast_df[cast_df["chrom"].isin(chrs_)]
s129_df = s129_df[s129_df["chrom"].isin(chrs_)]
fig,ax = plt.subplots(figsize=(8,6))
ax.hist(cast_df["per"])
fig,ax = plt.subplots(figsize=(8,6))
ax.hist(s129_df["per"])

cast_met_df = cast_df.query('per > 50')
s129_met_df = s129_df.query('per > 50')
cast_met_pb = pb.BedTool.from_dataframe(cast_met_df[["chrom","start","end"]])
s129_met_pb = pb.BedTool.from_dataframe(s129_met_df[["chrom","start","end"]])
ase_path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"

ase_df = pd.read_csv(ase_path, sep="\t")
ase_df.drop_duplicates(inplace=True)
ase_df.rename(columns={"start_transcript": "start","end_transcript":"end"},inplace=True)
ase_df = ase_df[ase_df.chrom.isin(chrs_)]
ase_df["size"] = abs(ase_df.end - ase_df.start)

#%% 
#we will get possible CG for each gene/promoter...
def findall(p, s):
    '''Yields all the positions of
    the pattern p in the string s.'''
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)

def get_CG (sequence):
    return [(i, sequence[i:i+2]) for i in findall('CG', sequence)]

def add_cg_to_df(df,fasta,start,end):
    df_pb = pb.BedTool.from_dataframe(df[["chrom",start,end]])
    seq = df_pb.sequence(fi=fasta)

    CG_numbers = []
    with open(seq.seqfn) as stdin:
        for line in stdin:
            if line[0] != ">":
                cg_list = get_CG(line)
                n_cg = len(cg_list)*2 ## actually both strands can be methylated wit CpG
                CG_numbers.append(n_cg)
    return CG_numbers
def annotate_master_table(df,fasta,normalize="CG"):
    df = df[df["size"] > 2000]
    promoter_df = df[["chrom","TSS_start","TSS_end","gene_name"]]

    #promoter_df = promoter_df.rename(columns={"TSS_start":"start","TSS_end":"end"})
    # or the whole gene (for LncRNA better) 
    # dont take into account the TSS
    body_negative_df = df.query("strand == '-'")
    body_negative_df = body_negative_df[["chrom","start","TSS_start","gene_name","gene_id","transcript_id"]]
    body_negative_df = body_negative_df.rename(columns={"TSS_start":"end"})
    body_positive_df = df.query("strand == '+'")
    body_positive_df = body_positive_df[["chrom","TSS_end","end","gene_name","gene_id","transcript_id"]]
    body_positive_df = body_positive_df.rename(columns={"TSS_end":"start"})
    body_df = pd.concat([body_positive_df,body_negative_df])
    body_df["body_size"] = body_df["end"].values - body_df["start"].values

    sizes = promoter_df["TSS_end"].values - promoter_df["TSS_start"].values
    sizes = [abs(x) for x in sizes]

    df_pb = pb.BedTool.from_dataframe(promoter_df)
    df_cast_met = df_pb.intersect(cast_met_pb,c=True).to_dataframe(names=["chrom","TSS_start","TSS_end","gene_name","DNAmet_CAST_promoter"])
    df = pd.merge(df,df_cast_met,on=["chrom","TSS_start","TSS_end","gene_name"])
    df_s129_met = df_pb.intersect(s129_met_pb,c=True).to_dataframe(names=["chrom","TSS_start","TSS_end","gene_name","DNAmet_S129_promoter"])
    df = pd.merge(df,df_s129_met,on=["chrom","TSS_start","TSS_end","gene_name"])
    promoter_df["all_cg_promoter"] = add_cg_to_df(promoter_df,fasta,"TSS_start","TSS_end")
    df = pd.merge(df,promoter_df,on=["chrom","TSS_start","TSS_end","gene_name"])

    df_body_pb = pb.BedTool.from_dataframe(body_df)
    df_cast_met_body = df_body_pb.intersect(cast_met_pb,c=True).to_dataframe(names=["chrom","start","end","gene_name","gene_id","transcript_id","body_size","DNAmet_CAST_body"])
    df = pd.merge(df,df_cast_met_body,on=["chrom","gene_name","gene_id","transcript_id"])
    df = df.rename(columns={"start_x":"start","end_x":"end"})
    df.drop(["start_y","end_y"],axis=1,inplace=True)
    df_s129_met_body = df_body_pb.intersect(s129_met_pb,c=True).to_dataframe(names=["chrom","start","end","gene_name","gene_id","transcript_id","body_size","DNAmet_S129_body"])
    df = pd.merge(df,df_s129_met_body,on=["chrom","gene_name","gene_id","transcript_id"])
    df = df.rename(columns={"start_x":"start","end_x":"end","body_size_x":"body_size"})
    df.drop(["start_y","end_y","body_size_y"],axis=1,inplace=True)
    body_df["all_cg_body"] = add_cg_to_df(body_df,fasta,"start","end")
    df = pd.merge(df,body_df,on=["chrom","gene_name","gene_id","transcript_id"])
    df = df.rename(columns={"start_x":"start","end_x":"end","body_size_x":"body_size"})
    df.drop(["start_y","end_y","body_size_y"],axis=1,inplace=True)

    print (df)

    if normalize == "CG":
        #normalized by possible CGs in the region
        ax.set_ylabel("Differential (CAST-S129) Percentage of possible methylated CpGs")

        df["DNAmet_promoter_CAST_ratio"] = (df.DNAmet_CAST_promoter+1)/(df.all_cg_promoter+1)
        df["DNAmet_promoter_S129_ratio"] = (df.DNAmet_S129_promoter+1)/(df.all_cg_promoter+1)
        df["DNAmet_promoter_ratio_difference"] = df.DNAmet_promoter_CAST_ratio - df.DNAmet_promoter_S129_ratio

        df["DNAmet_body_CAST_ratio"] = (df.DNAmet_CAST_body+1)/(df.all_cg_body+1)
        df["DNAmet_body_S129_ratio"] = (df.DNAmet_S129_body+1)/(df.all_cg_body+1)
        df["DNAmet_body_ratio_difference"] = df.DNAmet_body_CAST_ratio - df.DNAmet_body_S129_ratio

    # call the genes that are methylated

    cutoff_neg = (np.percentile(df.DNAmet_promoter_ratio_difference.values.tolist(), dnamet_cutoff))
    cutoff_pos = (np.percentile(df.DNAmet_promoter_ratio_difference.values.tolist(), 100-dnamet_cutoff))
    print ("cutoff pos = {}".format(cutoff_pos))
    print ("cutoff neg = {}".format(cutoff_neg))

    def call_DNAmet_genes(row):  #1 means that passes the cutoff, 2 means that it a potential ASE candidate
        if row["size"] > 2000:
            if row['TPM_transcript'] >= 1:
                if row['number_SNPs'] > 0:
                    if row["p_adj_sum"] <= 0.05:
                        if row["log2foldchange"] <= -1:
                            if row["DNAmet_promoter_ratio_difference"] >= cutoff_pos:
                                val = 2
                            else:
                                val = 0
                        elif row["log2foldchange"] >= 1:
                            if row["DNAmet_promoter_ratio_difference"] <= cutoff_neg:
                                val = 2
                            else:
                                val = 0
                        else:
                            if row["DNAmet_promoter_ratio_difference"] <= cutoff_neg or row["DNAmet_promoter_ratio_difference"] >= cutoff_pos:
                                val = 1
                            else:
                                val = 0
                    else:
                        if row["DNAmet_promoter_ratio_difference"] <= cutoff_neg or row["DNAmet_promoter_ratio_difference"] >= cutoff_pos:
                            val = 1
                        else:
                            val = 0
                else:
                    if row["DNAmet_promoter_ratio_difference"] <= cutoff_neg or row["DNAmet_promoter_ratio_difference"] >= cutoff_pos:
                        val = 1
                    else:
                        val = 0
            else:
                if row["DNAmet_promoter_ratio_difference"] <= cutoff_neg or row["DNAmet_promoter_ratio_difference"] >= cutoff_pos:
                    val = 1
                else:
                    val = 0
        else: 
            val = 0
        return val

    df['by_DNAmet'] = df.apply(call_DNAmet_genes, axis=1)


    
    df.rename(columns={"start": "start_transcript","end":"end_transcript"},inplace=True)
    df = df.drop_duplicates("gene_name_type",keep="last")
    print (len(df))
    df.to_csv(path_candidates, index=False,sep="\t")
    df.query('by_DNAmet == 1 or by_DNAmet == 2').to_csv("{}DNA_met_candidates_{}.bed".format(save_path,dnamet_cutoff))
    df.query('by_DNAmet == 2').to_csv("{}DNA_met_ASE_final_candidates_{}.bed".format(save_path,dnamet_cutoff))
    return df

def get_overlaps(regions):
    exp_n = len(regions.intersect(exp_no_ASE_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    cast_n = len(regions.intersect(cast_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    s129_n = len(regions.intersect(s129_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    cast_atac_n = len(regions.intersect(cast_atac_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    s129_atac_n = len(regions.intersect(s129_atac_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    common_atac_n = len(regions.intersect(common_atac_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    cast_ctcf_n = len(regions.intersect(cast_ctcf_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    s129_ctcf_n = len(regions.intersect(s129_ctcf_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    common_ctcf_n = len(regions.intersect(common_ctcf_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    cast_k4me3_n = len(regions.intersect(cast_h3k4me3_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    s129_k4me3_n = len(regions.intersect(s129_h3k4me3_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    common_k4me3_n = len(regions.intersect(common_h3k4me3_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    cast_k27ac_n = len(regions.intersect(cast_h3k27ac_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    s129_k27ac_n = len(regions.intersect(s129_h3k27ac_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    common_k27ac_n = len(regions.intersect(common_h3k27ac_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    k27me3_n = len(regions.intersect(k27me3_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    s5p_n = len(regions.intersect(s5p_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    s7p_n = len(regions.intersect(s7p_pb,u=True).to_dataframe(names=["chrom","start","end"]))
    cpg_n = len(regions.intersect(cpg_pb,u=True).to_dataframe(names=["chrom","start","end"]))

    return ([[cast_n,cast_atac_n,cast_ctcf_n,cast_k4me3_n,cast_k27ac_n,s5p_n,0],[s129_n,s129_atac_n,s129_ctcf_n,s129_k4me3_n,s129_k27ac_n,s7p_n,0],[exp_n,common_atac_n,common_ctcf_n,common_k4me3_n,common_k27ac_n,k27me3_n,cpg_n]])

def plot_overlaps(ll,size,name):
    lcast = [int(x)/size if int(x) != 0 else 0 for x in ll[0]]
    ls129 = [int(x)/size if int(x) != 0 else 0 for x in ll[1]]
    lcommon = [int(x)/size if int(x) != 0 else 0 for x in ll[2]]
    fig,ax = plt.subplots(1,1,figsize=(10,8))
    width = 0.25  # the width of the bars
    ticks = np.arange(len(lcast))
    bp1 = ax.bar(ticks-width,lcast,width,label="CAST",color=cast_color)
    bp1 = ax.bar(ticks,ls129,width,label="S129",color=s129_color)
    bp1 = ax.bar(ticks+width,lcommon,width,label="Common",color=both_color)
    plt.legend()
    ax.set_xticks(range(len(lcast)))
    ax.set_xticklabels(["Gene promoter","ATAC","CTCF","H3K4me3","H3K27ac","S5/S7/K27me3","CpG"])
    ax.set_xlabel("Overlap with")
    ax.set_title("Percentage of {} ({})".format(name,size))
    for i, v in enumerate(lcast):
        ax.text(i - width-0.05, v + .005, str(ll[0][i]), color=cast_color, fontweight='bold')
    for i, v in enumerate(ls129):
        ax.text(i-0.05, v + .005, str(ll[1][i]), color=s129_color, fontweight='bold')
    for i, v in enumerate(lcommon):
        ax.text(i + width-0.05, v + .005, str(ll[2][i]), color=both_color, fontweight='bold')
def prepare_df_for_pybed(ase_df): 
    #CHECK COLS
    # chrom, start, end, strand, gene_name
    #set gene_name column in the 5th place and change chrom column name
    cols = list(ase_df.columns.values)
    #cols = [[cols[0]],cols[3:5],cols[1:3],cols[5:]]
    cols = [[cols[0]],cols[3:5],[cols[5]],[cols[7]],cols[1:3],[cols[6]],cols[8:]]
    cols = reduce(operator.concat, cols)
    ase_df = ase_df[cols]
    ase_df.rename(columns={"chr": "chrom"},inplace=True)
    ase_df = ase_df.sort_values(by=['chrom', 'start'])
    return ase_df
if calculate:
    df = annotate_master_table(ase_df,s129_fasta,"CG")
else:
    df = pd.read_csv(path_candidates,sep="\t")

cutoff_neg = (np.percentile(df.DNAmet_promoter_ratio_difference.values.tolist(), dnamet_cutoff))
cutoff_pos = (np.percentile(df.DNAmet_promoter_ratio_difference.values.tolist(), 100-dnamet_cutoff))
print ("cutoff pos = {}".format(cutoff_pos))
print ("cutoff neg = {}".format(cutoff_neg))

print (df)
#%%

ase_df = df.copy() #df is just ase_df with more annotations
ase_df.rename(columns={"start_transcript": "start","end_transcript":"end"},inplace=True)

#cast_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & log2foldchange >= 1 & (p_adj_sum <= 0.05 | p_adj_mean <= 0.05)')
cast_df = ase_df.query('TPM_transcript >= 1 & log2foldchange >= 1 & p_adj_sum < 0.05')
cast_df.reset_index(inplace=True)
cast_df = prepare_df_for_pybed(cast_df)
cast_df = cast_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
cast_df = cast_df[list(cast_df.columns.values)[1:]]

cast_not_exp_df = ase_df.query('TPM_transcript < 1 & log2foldchange >= 1 & p_adj_sum < 0.05')
cast_not_exp_df.reset_index(inplace=True)
cast_not_exp_df = prepare_df_for_pybed(cast_not_exp_df)
cast_not_exp_df = cast_not_exp_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
cast_not_exp_df = cast_not_exp_df[list(cast_not_exp_df.columns.values)[1:]]

#cast_lnc_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & log2foldchange >= 1 & (p_adj_sum <= 0.05 | p_adj_mean <= 0.05) & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
cast_lnc_df = ase_df.query('TPM_transcript >= 1 & log2foldchange >= 1 & p_adj_sum < 0.05 & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
cast_lnc_df.reset_index(inplace=True)
cast_lnc_df = prepare_df_for_pybed(cast_lnc_df)
cast_lnc_df = cast_lnc_df.sort_values(by=["chrom","start"])
cast_lnc_df = cast_lnc_df[list(cast_lnc_df.columns.values)[1:]]

#cast_prot_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & log2foldchange >= 1 & (p_adj_sum <= 0.05 | p_adj_mean <= 0.05) & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
cast_prot_df = ase_df.query('TPM_transcript >= 1 & log2foldchange >= 1 & p_adj_sum < 0.05  & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
cast_prot_df.reset_index(inplace=True)
cast_prot_df = prepare_df_for_pybed(cast_prot_df)
cast_prot_df = cast_prot_df.sort_values(by=["chrom","start"])
cast_prot_df = cast_prot_df[list(cast_prot_df.columns.values)[1:]]

#s129_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & log2foldchange <= -1 & (p_adj_sum <= 0.05 | p_adj_mean <= 0.05)')
s129_df = ase_df.query('TPM_transcript >= 1  & log2foldchange <= -1 & p_adj_sum < 0.05 ')
s129_df.reset_index(inplace=True)
s129_df = prepare_df_for_pybed(s129_df)
s129_df = s129_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
s129_df = s129_df[list(s129_df.columns.values)[1:]]

s129_not_exp_df = ase_df.query('TPM_transcript <1 & log2foldchange <= -1 & p_adj_sum < 0.05 ')
s129_not_exp_df.reset_index(inplace=True)
s129_not_exp_df = prepare_df_for_pybed(s129_not_exp_df)
s129_not_exp_df = s129_not_exp_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
s129_not_exp_df = s129_not_exp_df[list(s129_not_exp_df.columns.values)[1:]]

#s129_lnc_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & log2foldchange <= -1 & (p_adj_sum <= 0.05 | p_adj_mean <= 0.05) & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
s129_lnc_df = ase_df.query('TPM_transcript >= 1  & log2foldchange <= -1 & p_adj_sum < 0.05  & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
s129_lnc_df.reset_index(inplace=True)
s129_lnc_df = prepare_df_for_pybed(s129_lnc_df)
s129_lnc_df = s129_lnc_df.sort_values(by=["chrom","start"])
s129_lnc_df = s129_lnc_df[list(s129_lnc_df.columns.values)[1:]]

#s129_prot_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & log2foldchange <= -1 & (p_adj_sum <= 0.05 | p_adj_mean <= 0.05) & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
s129_prot_df = ase_df.query('TPM_transcript >= 1  & log2foldchange <= -1 & p_adj_sum < 0.05  & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
s129_prot_df.reset_index(inplace=True)
s129_prot_df = prepare_df_for_pybed(s129_prot_df)
s129_prot_df = s129_prot_df.sort_values(by=["chrom","start"])
s129_prot_df = s129_prot_df[list(s129_prot_df.columns.values)[1:]]

#exp_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1)')
exp_df = ase_df.query('TPM_transcript >= 1')
exp_df.reset_index(inplace=True)
exp_df = prepare_df_for_pybed(exp_df)
exp_df = exp_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
exp_df = exp_df[list(exp_df.columns.values)[1:]]

#exp_no_ASE_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & (log2foldchange > -1 & log2foldchange < 1 )')
exp_no_ASE_df = ase_df.query('TPM_transcript >= 1  & (log2foldchange > -1 & log2foldchange < 1 )')
exp_no_ASE_df.reset_index(inplace=True)
exp_no_ASE_df = prepare_df_for_pybed(exp_no_ASE_df)
exp_no_ASE_df = exp_no_ASE_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
exp_no_ASE_df = exp_no_ASE_df[list(exp_no_ASE_df.columns.values)[1:]]

#exp_lnc_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
exp_lnc_df = ase_df.query('TPM_transcript >= 1  & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
exp_lnc_df.reset_index(inplace=True)
exp_lnc_df = prepare_df_for_pybed(exp_lnc_df)
exp_lnc_df = exp_lnc_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
exp_lnc_df = exp_lnc_df[list(exp_lnc_df.columns.values)[1:]]

#exp_lnc_no_ASE_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & (log2foldchange > -1 & log2foldchange < 1 ) & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
exp_lnc_no_ASE_df = ase_df.query('TPM_transcript >= 1 & (log2foldchange > -1 & log2foldchange < 1 ) & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
exp_lnc_no_ASE_df.reset_index(inplace=True)
exp_lnc_no_ASE_df = prepare_df_for_pybed(exp_lnc_no_ASE_df)
exp_lnc_no_ASE_df = exp_lnc_no_ASE_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
exp_lnc_no_ASE_df = exp_lnc_no_ASE_df[list(exp_lnc_no_ASE_df.columns.values)[1:]]

#exp_prot_no_ASE_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & (log2foldchange > -1 & log2foldchange < 1 ) & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
exp_prot_no_ASE_df = ase_df.query('TPM_transcript >= 1 & (log2foldchange > -1 & log2foldchange < 1 ) & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
exp_prot_no_ASE_df.reset_index(inplace=True)
exp_prot_no_ASE_df = prepare_df_for_pybed(exp_prot_no_ASE_df)
exp_prot_no_ASE_df = exp_prot_no_ASE_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
exp_prot_no_ASE_df = exp_prot_no_ASE_df[list(exp_prot_no_ASE_df.columns.values)[1:]]

#exp_prot_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
exp_prot_df = ase_df.query('TPM_transcript >= 1  & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
exp_prot_df.reset_index(inplace=True)
exp_prot_df = prepare_df_for_pybed(exp_prot_df)
exp_prot_df = exp_prot_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
exp_prot_df = exp_prot_df[list(exp_prot_df.columns.values)[1:]]

#not_exp_df = ase_df.query('(TPM_sumexonlength < 1 & TPM_transcriptlength < 1)')
not_exp_df = ase_df.query('TPM_transcript < 1 ')
#not_exp_df = ase_df.query('TPM_transcript < 1 and Promoter_state != "Inactive"')
not_exp_df.reset_index(inplace=True)
not_exp_df = prepare_df_for_pybed(not_exp_df)
not_exp_df = not_exp_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
not_exp_df = not_exp_df[list(not_exp_df.columns.values)[1:]]

#not_exp_lnc_df = ase_df.query('(TPM_sumexonlength < 1 & TPM_transcriptlength < 1) & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
not_exp_lnc_df = ase_df.query('TPM_transcript < 1  & transcript_biotype == "lncRNA" & gene_biotype == "lncRNA"')
not_exp_lnc_df.reset_index(inplace=True)
not_exp_lnc_df = prepare_df_for_pybed(not_exp_lnc_df)
not_exp_lnc_df = not_exp_lnc_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
not_exp_lnc_df = not_exp_lnc_df[list(not_exp_lnc_df.columns.values)[1:]]

#not_exp_prot_df = ase_df.query('(TPM_sumexonlength < 1 & TPM_transcriptlength < 1) & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
not_exp_prot_df = ase_df.query('TPM_transcript < 1 & transcript_biotype == "protein_coding" & gene_biotype == "protein_coding"')
not_exp_prot_df.reset_index(inplace=True)
not_exp_prot_df = prepare_df_for_pybed(not_exp_prot_df)
not_exp_prot_df = not_exp_prot_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
not_exp_prot_df = not_exp_prot_df[list(not_exp_prot_df.columns.values)[1:]]

#exp_snp_df = ase_df.query('(TPM_sumexonlength >= 1 | TPM_transcriptlength >= 1) & number_SNPs >= 1')
exp_snp_df = ase_df.query('TPM_transcript >= 1  & number_SNPs >= 1')
exp_snp_df.reset_index(inplace=True)
exp_snp_df = prepare_df_for_pybed(exp_snp_df)
exp_snp_df = exp_snp_df.sort_values(by=['gene_name', 'transcript_biotype']).drop_duplicates("gene_name",keep="last").sort_values(by=["chrom","start"])
exp_snp_df = exp_snp_df[list(exp_snp_df.columns.values)[1:]]

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

k27me3 = "F123.K27me3.BCP_peaks.HM_mode.bed"
s5p = "F123.S5p.BCP_peaks.HM_mode.bed"
s7p = "F123.S7p.BCP_peaks.HM_mode.bed"

cast_pb = pb.BedTool.from_dataframe(cast_df)
s129_pb = pb.BedTool.from_dataframe(s129_df)
exp_pb = pb.BedTool.from_dataframe(exp_df)
exp_no_ASE_pb = pb.BedTool.from_dataframe(exp_no_ASE_df)
exp_snp_pb = pb.BedTool.from_dataframe(exp_snp_df)
not_exp_pb = pb.BedTool.from_dataframe(not_exp_df)
cast_atac_pb = pb.BedTool(cast_atac).sort()
s129_atac_pb = pb.BedTool(s129_atac).sort()
common_atac_pb = pb.BedTool(common_atac).sort()
cast_ctcf_pb = pb.BedTool(cast_ctcf).sort()
s129_ctcf_pb = pb.BedTool(s129_ctcf).sort()
common_ctcf_pb = pb.BedTool(common_ctcf).sort()
cast_h3k4me3_pb = pb.BedTool(cast_h3k4me3).sort()
s129_h3k4me3_pb = pb.BedTool(s129_h3k4me3).sort()
common_h3k4me3_pb = pb.BedTool(common_h3k4me3).sort()
cast_h3k27ac_pb = pb.BedTool(cast_h3k27ac).sort()
s129_h3k27ac_pb = pb.BedTool(s129_h3k27ac).sort()
common_h3k27ac_pb = pb.BedTool(common_h3k27ac).sort()

k27me3_pb = pb.BedTool(k27me3).sort()
s5p_pb = pb.BedTool(s5p).sort()
s7p_pb = pb.BedTool(s7p).sort()

cpg_path = root+"CpG_mm10.table"
cpg_df = pd.read_csv(cpg_path, sep="\t")
cpg_df = cpg_df[cpg_df["chrom"].isin(chrs_)]
cpg_pb = pb.BedTool.from_dataframe(cpg_df[["chrom","start","end"]]).sort()

def set_only_promoters(aux_df,slop):
    df = aux_df.copy()
    mask_pos = (df["strand"] == "+")
    df_pos = df[mask_pos]
    df.loc[mask_pos,'start'] = df_pos["start"]-slop[0]
    df.loc[mask_pos,'end'] = df_pos["start"]+slop[1]
    mask_neg = (df["strand"] == "-")
    df_neg = df[mask_neg]
    df.loc[mask_neg,'end'] = df_neg["end"]+slop[0]
    df.loc[mask_neg,'start'] = df_neg["end"]-slop[1]
    return df

cast_promoter_df = set_only_promoters(cast_df,promoter)[["chrom","start","end"]]
s129_promoter_df = set_only_promoters(s129_df,promoter)[["chrom","start","end"]]
exp_promoter_df = set_only_promoters(exp_df,promoter)[["chrom","start","end"]]
cast_promoter_pb = pb.BedTool.from_dataframe(cast_promoter_df)
s129_promoter_pb = pb.BedTool.from_dataframe(s129_promoter_df)
exp_promoter_pb = pb.BedTool.from_dataframe(exp_promoter_df)




#%%
#scatter_plot log2fc and methylation
import itertools
from scipy.stats import norm
import scipy.stats
import os

#set promoters
promoter_up = 200 #100
promoter_down = 400 #200
promoter = [promoter_up,promoter_down]
 
## delete only if file exists ##
if os.path.exists(path_candidates):
    os.remove(path_candidates)
else:
    print("Sorry, I can not remove %s file." % path_candidates)
all_n_of_cpgs = []

import  math 
imp_file = "{}imprinted_genes_mouse_updated.csv".format(root)
imp_df = pd.read_csv(imp_file,sep="\t")
imp_genes = imp_df.query('Status == "Imprinted"')["Gene"].values.tolist()
imp_aliases = imp_df.query('Status == "Imprinted"')["Aliases"].values.tolist()
new_imp_aliases = []
for alias in imp_aliases:
    if type(alias) != float:
        aux = alias.split(",")
        for n,gene in enumerate(aux):
            if n == len(aux)-1:
                new_imp_aliases.append(gene[1:-1])
            elif n != 0:
                new_imp_aliases.append(gene[1:])
            else:
                new_imp_aliases.append(gene)
imp_genes = imp_genes+new_imp_aliases


def do_scatter_plot(df,ax,color,fasta,type_,annotate=False,normalize="CG"):
    print (len(df))

    df = df[df["size"] > 2000]
    print (len(df))

    log2fc_values = df.log2foldchange.values
    aux = []
    for x in log2fc_values:
        if x == np.inf:
            aux.append(13)
        elif x == -np.inf:
            aux.append(-13)
        else:
            aux.append(x)
    log2fc_values = aux
    #log2fc_values = df.ASE_ratio.values

    # if we want only promoters
    promoter_df = df[["chrom","TSS_start","TSS_end","gene_name"]]

    # or the whole gene (for LncRNA better) 
    # dont take into account the TSS
    body_negative_df = df.query("strand == '-'")
    body_negative_df = body_negative_df[["chrom","start","TSS_start","gene_name","gene_id","transcript_id"]]
    body_negative_df = body_negative_df.rename(columns={"TSS_start":"end"})
    body_positive_df = df.query("strand == '+'")
    body_positive_df = body_positive_df[["chrom","TSS_end","end","gene_name","gene_id","transcript_id"]]
    body_positive_df = body_positive_df.rename(columns={"TSS_end":"start"})
    body_df = pd.concat([body_positive_df,body_negative_df])
    body_df["body_size"] = body_df["end"].values - body_df["start"].values



    sizes = promoter_df["TSS_end"].values - promoter_df["TSS_start"].values
    sizes = [abs(x) for x in sizes]


    print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    print (df.columns)
    if normalize == "CG":
        #normalized by possible CGs in the region
        ax.set_ylabel("Differential (CAST-S129) Percentage of possible methylated CpGs")
        annotate_limit = 0.1666
        annotate_limit = 0.2

        df["DNAmet_promoter_CAST_ratio"] = (df.DNAmet_CAST_promoter+1)/(df.all_cg_promoter+1)
        df["DNAmet_promoter_S129_ratio"] = (df.DNAmet_S129_promoter+1)/(df.all_cg_promoter+1)
        df["DNAmet_promoter_ratio_difference"] = df.DNAmet_promoter_CAST_ratio - df.DNAmet_promoter_S129_ratio

        df["DNAmet_body_CAST_ratio"] = (df.DNAmet_CAST_body+1)/(df.all_cg_body+1)
        df["DNAmet_body_S129_ratio"] = (df.DNAmet_S129_body+1)/(df.all_cg_body+1)
        df["DNAmet_body_ratio_difference"] = df.DNAmet_body_CAST_ratio - df.DNAmet_body_S129_ratio

    elif normalize == "size":
 
        pass

    filter = False
    if filter:
        ##filter for plotting
        #log2fc_values_plot = []
        #difference_ratios_plot = []
        #n_of_cpgs_plot = []
        #names = []
        #chrs = []
        #starts = []
        #ends = []
        #strands = []
        #print ("@@@@@@@@@@")
        #print (len(log2fc_values))
        #for n,x in enumerate(zip(log2fc_values,difference_ratios,difference_ratios_body,n_of_cpgs,df.gene_name,df.chrom, df.start, df.end, df.strand)):
        #    if abs(x[1]-x[2]) > 0.1:
        #        if x[1] > 0.05 or x[1] < -0.05:
        #            log2fc_values_plot.append(x[0])
        #            difference_ratios_plot.append(x[1])
        #            n_of_cpgs_plot.append(x[3])
        #            names.append(x[4])
        #            chrs.append(x[5])
        #            starts.append(x[6])
        #            ends.append(x[7])
        #            strands.append(x[8])
        #print (len(log2fc_values_plot))
        #print ("@@@@@@@@@@")
        #log2fc_values = log2fc_values_plot
        #difference_ratios = difference_ratios_plot
        #n_of_cpgs = n_of_cpgs_plot
        pass
    else:
        #df.to_csv(path_candidates,sep="\t",mode="a")
        #names = df.gene_name.values
        #chrs = df.chrom.values
        #starts = df.start.values
        #ends = df.end.values
        #strands = df.strand.values
        pass
    print (len(df))
    df["type"] = len(df)*[type_]
    if not os.path.isfile(path_candidates):
        df.to_csv(path_candidates, index=False,sep="\t")
    else:  # else it exists so append without writing the header
        df.to_csv(path_candidates, mode='a', header=False, index=False,sep="\t")
    ax.axhline(cutoff_neg,color="black")
    ax.axhline(cutoff_pos,color="black")
    grey_df = df.query('by_DNAmet != 2')
    not_grey_df = df.query('by_DNAmet == 2')
    print ("Number of ASE genes passing cutoff: {}, {}".format(cutoff_neg, cutoff_pos))
    print (len(not_grey_df))
    sc = ax.scatter(grey_df.log2foldchange,grey_df.DNAmet_promoter_ratio_difference,color=both_color)
    sc = ax.scatter(not_grey_df.log2foldchange,not_grey_df.DNAmet_promoter_ratio_difference,color=color,edgecolors="black")
    handles, labels = sc.legend_elements(prop="sizes", alpha=0.6,num=4)
    ax.legend(handles, labels)
    ax.set_xlabel("log2FC")
    ax.yaxis.grid(True)

        
    if annotate:
        annotate_df = not_grey_df.query("((log2foldchange <= -1 or log2foldchange >= 1) and (DNAmet_promoter_ratio_difference >= {} or DNAmet_promoter_ratio_difference <= {})) or (log2foldchange <= -5 or log2foldchange >= 7)".format(annotate_limit,annotate_limit*-1))
        annotate_df = annotate_df.loc[~annotate_df.gene_name.str.startswith('Gm')]
        x = annotate_df.log2foldchange.values
        y = annotate_df.DNAmet_promoter_ratio_difference.values
        z = annotate_df.gene_name.values
        for (a,b,c) in zip(x,y,z):
            #if c in imp_genes:
            #    ax.annotate(c, xy=(a, b))
            #    print (c)
            ax.annotate(c, xy=(a, b))
        print (z)




#%%

#fig,ax = plt.subplots(figsize=(18,12))
fig = plt.figure(figsize=(8,6))
ax1 = plt.subplot2grid((1,6),(0,0), colspan=5)
ax2 = plt.subplot2grid((1,6),(0,5),sharey=ax1)

shift = False
if shift:
    exp_no_ASE_df["start"] = exp_no_ASE_df["start"]+100000
    cast_df["start"] = cast_df["start"]+100000
    s129_df["start"] = s129_df["start"]+100000
    exp_no_ASE_df["end"] = exp_no_ASE_df["end"]+100000
    cast_df["end"] = cast_df["end"]+100000
    s129_df["end"] = s129_df["end"]+100000


print (cast_df)
diff_met_cast = do_scatter_plot(cast_df,ax1,cast_color,cast_fasta,"cast",True,"CG")
print (len(diff_met_cast))
diff_met_s129 = do_scatter_plot(s129_df,ax1,s129_color,s129_fasta,"s129",True,"CG")
print (len(diff_met_s129))

bp = ax2.boxplot([diff_met_cast,diff_met_s129],
            showfliers=True,
            patch_artist=True)
colors = ["darkgreen", cast_color, s129_color,"grey"]
colors = [cast_color, s129_color,"grey"]
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
ax2.yaxis.grid(True)
ax2.legend([bp["boxes"][0], bp["boxes"][1]], ['CAST', 'S129'], loc='upper right')
ax1.set_title("Active genes (Promoters)")
plt.savefig(save_path+"DNA_met_scatter.pdf",dpi=300)
plt.savefig(save_path+"DNA_met_scatter.eps",dpi=300)
fig,ax = plt.subplots()
hist_list = np.concatenate([diff_met_cast,diff_met_s129])  
ax.hist(hist_list,bins=50)

cutoff_pos = (np.percentile(np.concatenate([diff_met_cast,diff_met_s129]), 100-dnamet_cutoff))
cutoff_neg = (np.percentile(np.concatenate([diff_met_cast,diff_met_s129]), dnamet_cutoff))
print ("cutoff pos = {}".format(cutoff_pos))
print ("cutoff neg = {}".format(cutoff_neg))


#%%

#candidates_df = pd.read_csv(path_candidates,sep="\t")
##candidates_df = candidates_df.query("size > 2000")
##df.replace([np.inf, -np.inf], np.nan,inplace=True)
#
##pseudo_count. I add 1 only if value == 0 for s5 and s7
##like df[S5p_CAST_reads"]+1
#candidates_df["S5p_CAST_reads"].replace(0,1,inplace=True)
#candidates_df["S5p_S129_reads"].replace(0,1,inplace=True)
#candidates_df["S7p_CAST_reads"].replace(0,1,inplace=True)
#candidates_df["S7p_S129_reads"].replace(0,1,inplace=True)
#candidates_df["K27_S7_CAST_ratio"] = candidates_df["K27me3_CAST_reads"]/candidates_df["S7p_CAST_reads"]
#candidates_df["K27_S7_S129_ratio"] = candidates_df["K27me3_S129_reads"]/candidates_df["S7p_S129_reads"]
##candidates_df["K27_S7_CAST_ratio"] = (candidates_df["K27me3_CAST_reads"]+1)/(candidates_df["S7p_CAST_reads"]+1)
##candidates_df["K27_S7_S129_ratio"] = (candidates_df["K27me3_S129_reads"]+1)/(candidates_df["S7p_S129_reads"]+1)
#
#candidates_df["K27_S7_CAST_ratio"].replace([np.inf, np.nan], 0,inplace=True)
#candidates_df["K27_S7_S129_ratio"].replace([np.inf, np.nan], 0,inplace=True)
#candidates_df["diff_PRC_ratio"] = candidates_df.K27_S7_CAST_ratio - candidates_df.K27_S7_S129_ratio
#df.replace([np.inf, -np.inf], np.nan,inplace=True)
#%%

# Code to add data to the ASE table

## merge with ASE_DF DNA MET DATA and prc ratios
#candidates_df.drop_duplicates(inplace=True)
#candidates_df = candidates_df.rename(columns={"start":"start_transcript","end":"end_transcript"})
##print (candidates_df.columns)
#gene_df = pd.read_csv(root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype.bed",sep="\t")
#merged_df = pd.merge(gene_df,candidates_df,on=gene_df.columns.values.tolist(),how="left")
#
#
##add ratios to all
#merged_df["S5p_CAST_reads"].replace(0,1,inplace=True)
#merged_df["S5p_S129_reads"].replace(0,1,inplace=True)
#merged_df["S7p_CAST_reads"].replace(0,1,inplace=True)
#merged_df["S7p_S129_reads"].replace(0,1,inplace=True)
#merged_df["K27_S7_CAST_ratio"] = merged_df["K27me3_CAST_reads"]/merged_df["S7p_CAST_reads"]
#merged_df["K27_S7_S129_ratio"] = merged_df["K27me3_S129_reads"]/merged_df["S7p_S129_reads"]
#merged_df["K27_S7_CAST_ratio"].replace([np.inf, np.nan], 0,inplace=True)
#merged_df["K27_S7_S129_ratio"].replace([np.inf, np.nan], 0,inplace=True)
#merged_df["diff_PRC_ratio"] = merged_df.K27_S7_CAST_ratio - merged_df.K27_S7_S129_ratio
#
#
#merged_df.to_csv(root+"/F123/tracks/300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated.bed",sep="\t")

#%%
