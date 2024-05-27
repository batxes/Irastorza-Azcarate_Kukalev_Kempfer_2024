
#%%

# Script that plots all figures related to ATAC-seq characterization and motifs related to genes
# hommer software was used in order to characterize the ATAC peaks
# First we need to run homer like this (examle):
# annotatePeaks.pl F123_ATAC_Bing_non_phased_peaks.bed mm10 > nonphased_ATAC.txt
# to get GO: annotatePeaks.pl CAST_specific_ATAC_in_promoters.bed mm10 -go CAST_specific_ATAC_in_promoters_GO

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import Counter
from  operator import itemgetter
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

ase_path = root+"/300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"
ase_df = pd.read_csv(ase_path,sep="\t")
ase_df = ase_df.drop_duplicates()
ase_df = ase_df[ase_df["chrom"].isin(chrs_)]

cast_genes = ase_df.query("type == 'cast'").gene_name.values.tolist()
s129_genes = ase_df.query("type == 's129'").gene_name.values.tolist()
exp_w_genes = ase_df.query("type == 'exp_w_SNP'").gene_name.values.tolist()
no_genes = ase_df.query("type == 'no'").gene_name.values.tolist()
exp_wo_genes = ase_df.query("type == 'exp_wo_SNP'").gene_name.values.tolist()

path = "{}homer/".format(root)

s129 = "S129_specific_ATAC.txt"
cast = "CAST_specific_ATAC.txt"
nonphased = "nonphased_ATAC.txt"
shared = "shared_ATAC.txt"

def type_(x):
    #name = x["Gene Name"]
    name = x
    value = "nan"
    if name in no_genes:
        value = "no"
    if name in exp_wo_genes:
        value = "exp_wo_SNP"
    if name in exp_w_genes:
        value = "exp_w_SNP"
    if name in s129_genes:
        value = "s129"
    if name in cast_genes:
        value = "cast"
    return value

def tpm_(x):
    #name = x["Gene Name"]
    name = x
    subset_df = ase_df[ase_df.gene_name == name]
    tpms = subset_df.TPM_transcript.values
    for tpm in tpms:
        if tpm >= 1:
            return True
        else:
            return False
    return False
 
types = ["cast","s129","exp_w_SNP","exp_wo_SNP","no","nan"]

def get_annotation(file):
    data = pd.read_csv(file,"\t")
    print ("Size of df")
    print (len(data))
    data['type'] = data["Gene Name"].apply(type_)
    data['tpm'] = data["Gene Name"].apply(tpm_)
    print ("Size of df with nan")
    print (len(data[data.type == "nan"]))
    print ("Size of df with tpm < 1")
    print (len(data[data.tpm < 1]))
    data.to_csv("{}.with_type.txt".format(file),sep="\t",index=False)

    final_annotations = []
    for t in types:
        aux = data.query("type == '{}'".format(t))
        annotations = aux.Annotation
        annotations_list = []
        for annot in annotations:
            annotations_list.append(annot.split(" ")[0])
        final_annotations.append(annotations_list)
    return final_annotations

annot = ["3'", "5'", 'Intergenic', 'TTS', 'exon', 'intron', 'non-coding', 'promoter-TSS']
dic = dict.fromkeys(annot,0)

def do_plot (annotations,title):
    c = [cast_color,s129_color,"darkgreen","limegreen","grey","black"]
    bottom = [0,0,0,0,0,0,0,0]
    totals = [0,0,0,0,0,0,0,0]
    bottom_ratio = [0,0,0,0,0,0,0,0]
    fig,ax = plt.subplots()
    for n,t in enumerate(types):
        annotations_set = Counter(annotations[n])
        keys = []
        numbers = []
        counter = 0
        for key,number in sorted(annotations_set.items()):
            while key != annot[counter]:
                dic[annot[counter]] = 0
                counter += 1
            dic[key] = number
            counter += 1
        for key,number in sorted(dic.items()):
            keys.append(key)
            numbers.append(number)
        print (keys,numbers,c[n],keys,bottom)
        ax.bar(keys,numbers,color = c[n],label=keys,bottom=bottom)
        bottom = [sum(x) for x in zip(bottom,numbers)]
        totals = bottom
    print (totals)
    ax.tick_params(rotation=45)
    ax.set_title(title)
    fig.savefig("{}ATAC_annotation_{}.pdf".format(save_path,title),dpi=300)
    
    fig,ax = plt.subplots()

    for n,t in enumerate(types):
        annotations_set = Counter(annotations[n])
        print (annotations_set)
        keys = []
        numbers = []
        counter = 0
        for key,number in sorted(annotations_set.items()):
            while key != annot[counter]:
                dic[annot[counter]] = 0/totals[counter]
                print ("{} {}: {} = {} / {}".format(t,annot[counter],dic[annot[counter]],0,totals[counter]))
                counter += 1
            dic[key] = number/totals[counter]
            print ("{} {}: {} = {} / {}".format(t,key,dic[key],number,totals[counter]))
            counter += 1
        for key,number in sorted(dic.items()):
            keys.append(key)
            numbers.append(number)
        print (keys,numbers,c[n],keys,bottom_ratio)
        ax.bar(keys,numbers,color = c[n],label=keys,bottom=bottom_ratio)
        bottom_ratio = [sum(x) for x in zip(bottom_ratio,numbers)]
    ax.tick_params(rotation=45)
    ax.set_title(title)
    fig.savefig("{}ATAC_annotation_normalized_{}.pdf".format(save_path,title),dpi=300)

    #ax.legend(title='Annotation')
#%%
do_plot(get_annotation(path+cast),"CAST_specific")
do_plot(get_annotation(path+s129),"S129_specific")
do_plot(get_annotation(path+nonphased),"non_phased")
do_plot(get_annotation(path+shared),"Shared")

#%%
#check size of peaks in different places:

aux = pd.read_csv(path+cast,sep="\t")
aux["size_"] = aux.End - aux.Start
cast_promoter = aux[aux["Annotation"].str.startswith('promoter-TSS')]
print (np.median(cast_promoter.size_))
cast_intron = aux[aux["Annotation"].str.startswith('intron')]
print (np.median(cast_intron.size_))
cast_inter = aux[aux["Annotation"].str.startswith('Intergenic')]
print (np.median(cast_inter.size_))

aux = pd.read_csv(path+s129,sep="\t")
aux["size_"] = aux.End - aux.Start
s129_promoter = aux[aux["Annotation"].str.startswith('promoter-TSS')]
print (np.median(s129_promoter.size_))
s129_intron = aux[aux["Annotation"].str.startswith('intron')]
print (np.median(s129_intron.size_))
s129_inter = aux[aux["Annotation"].str.startswith('Intergenic')]
print (np.median(s129_inter.size_))

aux = pd.read_csv(path+shared,sep="\t")
aux["size_"] = aux.End - aux.Start
shared_promoter = aux[aux["Annotation"].str.startswith('promoter-TSS')]
print (np.median(shared_promoter.size_))
shared_intron = aux[aux["Annotation"].str.startswith('intron')]
print (np.median(shared_intron.size_))
shared_inter = aux[aux["Annotation"].str.startswith('Intergenic')]
print (np.median(shared_inter.size_))

def plot_hist(values1,values2,values3,color1,color2,color3,title,bins):
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(8,8))
    _, bins, _ = ax1.hist(values1,color=color1,bins=bins,alpha=0.8)
    ax1.axvline(np.median(values1),color="black")
    ax1.axvline(np.mean(values1),color="red")
    ax2.hist(values2,color=color2,bins=bins,alpha=0.8)
    ax2.axvline(np.median(values2),color="black")
    ax2.axvline(np.mean(values2),color="red")
    ax3.hist(values3,color=both_color,bins=bins,alpha=0.8)
    ax3.axvline(np.median(values3),color="black")
    ax3.axvline(np.mean(values3),color="red")
    ax1.set_xlim([-100,2000])
    ax3.set_xlabel("Size in bp")
    ax1.set_title("CAST {} peaks. Median: {} Mean: {}".format(title,np.median(values1),np.mean(values1)))
    ax2.set_title("S129 {} peaks Median: {} Mean: {}".format(title,np.median(values2),np.mean(values2)))
    ax3.set_title("Common {} peaks Median: {} Mean: {}".format(title,np.median(values3),np.mean(values3)))
    fig.savefig("{}{}_peak_size.pdf".format(save_path,title),dpi=300)

plot_hist(cast_promoter.size_,s129_promoter.size_,shared_promoter.size_,cast_color,s129_color,both_color,"Promoters",100)
plot_hist(cast_intron.size_,s129_intron.size_,shared_intron.size_,cast_color,s129_color,both_color,"Intron",1000)
plot_hist(cast_inter.size_,s129_inter.size_,shared_inter.size_,cast_color,s129_color,both_color,"Intergenic",100)


# %%

# get specific data
aux = pd.read_csv("{}CAST_specific_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('promoter-TSS')][["Chr","Start","End"]]
aux.to_csv("{}CAST_specific_ATAC_in_Promoters.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}CAST_specific_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('promoter-TSS')]
aux = aux.query('type == "cast"')[["Chr","Start","End"]]
aux.to_csv("{}CAST_specific_cast_genes_ATAC_in_Promoters.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}S129_specific_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('promoter-TSS')][["Chr","Start","End"]]
aux.to_csv("{}S129_specific_ATAC_in_Promoters.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}S129_specific_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('promoter-TSS')]
aux = aux.query('type == "s129"')[["Chr","Start","End"]]
aux.to_csv("{}S129_specific_s129_genes_ATAC_in_Promoters.bed".format(path),sep="\t",index=False)
print (len(aux))

#intron
aux = pd.read_csv("{}CAST_specific_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('intron')][["Chr","Start","End"]]
aux.to_csv("{}CAST_specific_ATAC_in_intron.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}CAST_specific_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('intron')]
aux = aux.query('type == "cast"')[["Chr","Start","End"]]
aux.to_csv("{}CAST_specific_cast_genes_ATAC_in_intron.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}S129_specific_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('intron')][["Chr","Start","End"]]
aux.to_csv("{}S129_specific_ATAC_in_intron.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}S129_specific_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('intron')]
aux = aux.query('type == "s129"')[["Chr","Start","End"]]
aux.to_csv("{}S129_specific_s129_genes_ATAC_in_intron.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}CAST_specific_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('Intergenic')][["Chr","Start","End"]]
aux.to_csv("{}CAST_specific_ATAC_in_intergenic.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}CAST_specific_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('Intergenic')]
aux = aux.query('type == "cast"')[["Chr","Start","End"]]
aux.to_csv("{}CAST_specific_cast_genes_ATAC_in_intergenic.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}S129_specific_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('Intergenic')][["Chr","Start","End"]]
aux.to_csv("{}S129_specific_ATAC_in_intergenic.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}S129_specific_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('Intergenic')]
aux = aux.query('type == "s129"')[["Chr","Start","End"]]
aux.to_csv("{}S129_specific_s129_genes_ATAC_in_intergenic.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}shared_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('Intergenic')][["Chr","Start","End"]]
aux.to_csv("{}shared_ATAC_in_intergenic.bed".format(path),sep="\t",index=False)
print (len(aux))
aux = pd.read_csv("{}shared_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('intron')][["Chr","Start","End"]]
aux.to_csv("{}shared_ATAC_in_intron.bed".format(path),sep="\t",index=False)
print (len(aux))
aux = pd.read_csv("{}shared_ATAC.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('promoter-TSS')][["Chr","Start","End"]]
aux.to_csv("{}shared_ATAC_in_Promoters.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}shared_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('promoter-TSS')]
aux = aux.query('type == "cast"')[["Chr","Start","End"]]
aux.to_csv("{}shared_ATAC_cast_genes_in_Promoters.bed".format(path),sep="\t",index=False)
print (len(aux))

aux = pd.read_csv("{}shared_ATAC.txt.with_type.txt".format(path),sep="\t")
aux = aux[aux["Annotation"].str.startswith('promoter-TSS')]
aux = aux.query('type == "s129"')[["Chr","Start","End"]]
aux.to_csv("{}shared_ATAC_s129_genes_in_Promoters.bed".format(path),sep="\t",index=False)
print (len(aux))

"""
Then we launch this:
findMotifsGenome.pl CAST_specific_ATAC_in_Promoters.bed mm10 CAST_specific_Promoters_MotifOutput/ -size 200 -mask
findMotifsGenome.pl CAST_specific_cast_genes_ATAC_in_Promoters.bed mm10 CAST_specific_cast_genes_Promoters_MotifOutput/ -size 200 -mask
findMotifsGenome.pl CAST_specific_ATAC_in_intron.bed mm10 CAST_specific_intron_MotifOutput/ -size 200 -mask
findMotifsGenome.pl CAST_specific_cast_genes_ATAC_in_intron.bed mm10 CAST_specific_cast_genes_intron_MotifOutput/ -size 200 -mask
findMotifsGenome.pl CAST_specific_ATAC_in_intergenic.bed mm10 CAST_specific_intergenic_MotifOutput/ -size 200 -mask
findMotifsGenome.pl CAST_specific_cast_genes_ATAC_in_intergenic.bed mm10 CAST_specific_cast_genes_intergenic_MotifOutput/ -size 200 -mask

findMotifsGenome.pl S129_specific_ATAC_in_Promoters.bed mm10 S129_specific_Promoters_MotifOutput/ -size 200 -mask
findMotifsGenome.pl S129_specific_ATAC_in_intron.bed mm10 S129_specific_intron_MotifOutput/ -size 200 -mask
findMotifsGenome.pl S129_specific_ATAC_in_intergenic.bed mm10 S129_specific_intergenic_MotifOutput/ -size 200 -mask


findMotifsGenome.pl S129_specific_s129_genes_ATAC_in_Promoters.bed mm10 S129_specific_S129_genes_Promoters_MotifOutput/ -size 200 -mask
findMotifsGenome.pl S129_specific_s129_genes_ATAC_in_intron.bed mm10 S129_specific_S129_genes_intron_MotifOutput/ -size 200 -mask
findMotifsGenome.pl S129_specific_s129_genes_ATAC_in_intergenic.bed mm10 S129_specific_S129_genes_intergenic_MotifOutput/ -size 200 -mask


findMotifsGenome.pl shared_ATAC_in_Promoter.bed mm10 shared_Promoters_MotifOutput/ -size 200 -mask
findMotifsGenome.pl shared_ATAC_in_intron.bed mm10 shared_intron_MotifOutput/ -size 200 -mask
findMotifsGenome.pl shared_ATAC_in_intergenic.bed mm10 shared_intergenic_MotifOutput/ -size 200 -mask

findMotifsGenome.pl shared_ATAC_cat_genes_in_Promoters.bed mm10 shared_cast_genes_Promoters_MotifOutput/ -size 200 -mask
findMotifsGenome.pl shared_ATAC_s129_genes_in_Promoters.bed mm10 shared_s129_genes_Promoters_MotifOutput/ -size 200 -mask

"""

# %%

import scipy.cluster.hierarchy as sch


# get matrix

folders = ["CAST_specific_Promoters_MotifOutput","CAST_specific_cast_genes_Promoters_MotifOutput",
           "S129_specific_Promoters_MotifOutput","S129_specific_S129_genes_Promoters_MotifOutput",
           "shared_Promoters_MotifOutput","shared_cast_genes_Promoters_MotifOutput","shared_s129_genes_Promoters_MotifOutput"]

# curated list extracted filtering the TFs for their expression
expressed_TFs = ["n-Myc","bhLHe40","Usf2","USF1","Stat3", "Six1","STAT5","STAT1","Rfx5","Rfx1","Pknox1","Oct6","Oct4","Nanog","NRF1","NFIL3","Maz","Max","MYB","MITF","JunD","JunB","Jun-AP1","IRF8","GRHL2","Egr1","E2F4","E2F3","CTCF","CLOCK","BORIS","BMYB","Atf7","Atf3","AP-1","Atf2","Zfp57","Rfx2","YY1","NRF","KLF10","ELF3","TFE3","ETV4","ZBTB33","Sp2","Sp1","Klf9","KLF6","KLF5","KLF3","KLF1","GABPA","Elk4","Elk1","ETV1","ELF1","Ronin","Elf4","ZNF143|STAF","Etv2","Klf4"]

all_lists = []
list_of_lists = []
for f in folders:
    df = pd.read_csv("{}{}/knownResults.txt".format(path,f),sep="\t")
    list_names = df.query('`q-value (Benjamini)` <= 0.05 and `P-value` <= 0.001')["Motif Name"].values.tolist()
    final_list = []
    for l in list_names:
        final_list.append(l.split("(")[0].lower())
    ase_df["gene_name"] = ase_df["gene_name"].str.lower()
    final_expressed_list = []
    for gene in final_list:
        tpms = ase_df[ase_df.gene_name == gene].TPM_transcript.values
        for tpm in tpms:
            if tpm >= 1:
                final_expressed_list.append(gene)
    all_lists = all_lists+final_expressed_list
    list_of_lists.append(final_expressed_list)

all_names = list(set(all_lists))
all_names.sort()

# Create an empty matrix to store presence (1) or absence (0) of names in lists
matrix = np.zeros((len(all_names), len(folders)), dtype=int)

# Fill the matrix with 1 if the name is present in the respective list
for n,f in enumerate(folders):
    for i, name in enumerate(all_names):
        if name in list_of_lists[n]:
            matrix[i, n] = 1

row_linkage = sch.linkage(matrix, method='ward', metric='euclidean')
col_linkage = sch.linkage(matrix.T, method='ward', metric='euclidean')

# Reorder rows and columns based on clustering
row_order = sch.dendrogram(row_linkage, no_plot=True)['leaves']

matrix = matrix[row_order, :]

# Create the heatmap

plt.figure(figsize=(8, 15))
plt.imshow(matrix, cmap='viridis', aspect='auto', interpolation='none')
plt.colorbar(label='Presence')
plt.xticks(range(len(folders)), folders)
plt.yticks(range(len(all_names)), [all_names[i] for i in row_order])
plt.title('Clustered Presence Heatmap')
plt.xlabel('Lists')
plt.ylabel('Names')
plt.xticks(rotation=45)
plt.savefig("{}atac_all.pdf".format(save_path),dpi=300)

plt.show()


# %%

# for Figure 4

# CTCF strand

######## annotation of CTCF peaks orientation

# annotatePeaks.pl {}homer/CTCF.common.peaks.bed mm10 -m ctcf.motif > ctcf_in_common.txt 

l_pos = []
l_neg = []
for path in ["ctcf_in_cast.txt","ctcf_in_s129.txt","ctcf_in_common.txt"]:

    df = pd.read_csv("{}homer/{}".format(root,path),sep="\t")
    df = df[["Chr","Start","End","CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer Distance From Peak(sequence,strand,conservation)"]]
    #remove peaks with no CTCF match
    df = df.dropna()

    def get_strand(r):
        p = r.split("(")[1]
        return p.split(",")[1]

    df["CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer Distance From Peak(sequence,strand,conservation)"] = df["CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer Distance From Peak(sequence,strand,conservation)"].apply(lambda x:get_strand(x))
    df = df.rename(columns={"CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer Distance From Peak(sequence,strand,conservation)":"Strand"}).sort_values(["Chr","Start"])
    print (df)
    df.to_csv("{}homer/motifs_{}".format(root,path),index=False,sep="\t")
    df.query('Strand == "+"').to_csv("{}homer/positive_strand_motifs_{}".format(root,path),index=False,sep="\t")
    df.query('Strand == "-"').to_csv("{}homer/negative_strand_motifs_{}".format(root,path),index=False,sep="\t")
    l_pos.append(df.query("Strand == '+'"))
    l_neg.append(df.query("Strand == '-'"))
neg = pd.concat(l_neg).sort_values(["Chr","Start"])[["Chr","Start","End"]]
pos = pd.concat(l_pos).sort_values(["Chr","Start"])[["Chr","Start","End"]]
neg.to_csv("{}homer/negative_strand_motifs_all.txt".format(root),index=False,sep="\t")
pos.to_csv("{}homer/positive_strand_motifs_all.txt".format(root),index=False,sep="\t")



# %%

# characterize ZFP57

## annotatePeaks.pl {}homer/F123_R1R2_nonphased_chroma_all_peaks.bed mm10 -m ZFP57.motif > zfp57_in_all_ATAC_peaks.txt

# ATAC are all the peaks, regardless of phased or not

df = pd.read_csv("{}homer/zfp57_in_all_ATAC_peaks.txt".format(root),sep="\t")
df = df[["Chr","Start","End","Annotation","Gene Name","Zfp57(Zf)/H1-ZFP57.HA-ChIP-Seq(GSE115387)/Homer Distance From Peak(sequence,strand,conservation)"]]
df = df.dropna()
df = df.rename(columns={"Gene Name":"gene_name"})
merged_df = ase_df.merge(df,on="gene_name",how="inner")
#sort by type
sorter = ["cast","s129","exp_w_SNP","exp_wo_SNP","no"]
merged_df.type = merged_df.type.astype("category")
merged_df.type = merged_df.type.cat.set_categories(sorter)
merged_df = merged_df.sort_values(["type"])
print(merged_df)
zfp_df = merged_df.drop_duplicates(subset="gene_name",keep="first")
print(zfp_df)
zfp_df.to_csv("/home/ibai/test.tsv",sep="\t")

print ("Total genes: {}".format(len(zfp_df)))
print ("CAST genes: {}".format(len(zfp_df[zfp_df.type == "cast"])))
print ("% of total CAST genes: {}".format(100*len(zfp_df[zfp_df.type == "cast"])/len(ase_df[ase_df.type == "cast"])))
print ("S129 genes: {}".format(len(zfp_df[zfp_df.type == "s129"])))
print ("% of total S129 genes: {}".format(100*len(zfp_df[zfp_df.type == "s129"])/len(ase_df[ase_df.type == "s129"])))
print ("EXP w SNP genes: {}".format(len(zfp_df[zfp_df.type == "exp_w_SNP"])))
print ("% of total exp_w_SNP genes: {}".format(100*len(zfp_df[zfp_df.type == "exp_w_SNP"])/len(ase_df[ase_df.type == "exp_w_SNP"])))
print ("EXP no SNP genes: {}".format(len(zfp_df[zfp_df.type == "exp_wo_SNP"])))
print ("% of total exp_wo_SNP genes: {}".format(100*len(zfp_df[zfp_df.type == "exp_wo_SNP"])/len(ase_df[ase_df.type == "exp_wo_SNP"])))
print ("Not exp genes: {}".format(len(zfp_df[zfp_df.type == "no"])))
print ("median gene expression of all: {}".format(zfp_df.TPM_transcript.median()))
print ("median gene expression of cast: {}".format(zfp_df[zfp_df.type == "cast"].TPM_transcript.median()))
print ("median gene expression of s129: {}".format(zfp_df[zfp_df.type == "s129"].TPM_transcript.median()))
print ("median gene expression of exp_w_SNP: {}".format(zfp_df[zfp_df.type == "exp_w_SNP"].TPM_transcript.median()))
print ("are these genes methylated: {}".format(len(zfp_df[zfp_df.by_DNAmet == 1])))
print ("Total number of methylated genes in the genome: {}".format(len(ase_df[ase_df.by_DNAmet == 1])))
print ("are these genes methylated and ASE: {}".format(len(zfp_df.query('by_DNAmet == 1 and (type == "cast" or type == "s129")'))))
print ("have these genes PRC: {}".format(len(zfp_df[zfp_df.H3K27me3 == 1])))

imp_file = "{}imprinted_genes_mouse_updated.csv".format(root)
imp_df = pd.read_csv(imp_file,sep="\t")
imp_genes = imp_df.query('Status == "Imprinted"')["Gene"].values.tolist()
print ("Number of known imprinted genes: {}".format(len(imp_df)))
for x in imp_df.query('Status == "Imprinted"')["Aliases"].values:
    y = x.split(",")
    for gene in y:
        imp_genes.append(gene.strip())

zfp_and_imprinted_df = zfp_df[zfp_df.gene_name.isin(imp_genes)]
print ("Number of ZFP related genes that are imprinted: {}".format(len(zfp_and_imprinted_df)))
print ("From these, how many are differentially methylated: {}".format(len(zfp_and_imprinted_df[zfp_and_imprinted_df.by_DNAmet == 1])))
print (zfp_and_imprinted_df[zfp_and_imprinted_df.by_DNAmet == 1])



