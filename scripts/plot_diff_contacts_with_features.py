#%%

# Figure4, Figure 4
# About Differential contacts EP and CTCF loops
# the scripts to generate the dataset for pairs/arcs are in the bottom of this script

from collections import Counter
import pybedtools as pb
from adjustText import adjust_text
import ast
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

distance = 2000000
recalculate_all = False

ase_path = root+"/300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"
ase_df = pd.read_csv(ase_path,sep="\t")
ase_df = ase_df.drop_duplicates()
ase_df = ase_df[ase_df["chrom"].isin(chrs_)]

use_TSS = False
if use_TSS:
    s_ = "TSS_start"
    e_ = "TSS_end"
else:
    s_ = "start_transcript"
    e_ = "end_transcript"

cast_df = ase_df.query("(log2foldchange >= 1) and p_adj_sum < 0.05 and TPM_transcript >= 1")[["chrom",s_,e_]]
s129_df = ase_df.query("(log2foldchange <= -1) and p_adj_sum < 0.05 and TPM_transcript >= 1")[["chrom",s_,e_]]
cast_prc_df = ase_df.query("(log2foldchange >= 1) and p_adj_sum < 0.05 and TPM_transcript >= 1 and K27me3_positive_TSS == 'Yes'")[["chrom",s_,e_]]
s129_prc_df = ase_df.query("(log2foldchange <= -1) and p_adj_sum < 0.05 and TPM_transcript >= 1 and K27me3_positive_TSS == 'Yes'")[["chrom",s_,e_]]
print (len(cast_df),len(cast_prc_df))
print (len(s129_df),len(s129_prc_df))
exp_df = ase_df.query("((log2foldchange < 1 and log2foldchange > -1) or p_adj_sum > 0.05) and TPM_transcript >= 1")[["chrom",s_,e_]]
if not recalculate_all:
    contacts_with_features_df = pd.read_csv(root+"full_contacts_with_features_and_genes_and_strands.tsv",sep="\t")


#%%

def flatten_list(lst):
    flattened = []
    lst = lst.tolist()
    for item in lst:
        if isinstance(item,float):
            continue
        items = item.split(",")
        if len(items) == 1:
            flattened.append(item)
        else:
            flattened.extend(items)
    return list(set(flattened))

def remove_nan_from_list(lst):
    filtered_list = [x for x in lst if not isinstance(x, float)]
    return filtered_list

def print_genes(df,type_,only_top5=False):
    if type_ == "cast_contacts" or type_ == "cast_specific":
        l1 = flatten_list(df.bin1_cast_gene_names)
        l2 = flatten_list(df.bin2_cast_gene_names)
    else:
        l1 = flatten_list(df.bin1_s129_gene_names)
        l2 = flatten_list(df.bin2_s129_gene_names)
    g = list(set(l1+l2))
    if not only_top5:
        print ("ALL {}: {}".format(len(g),g))
    aux_df = ase_df[ase_df["gene_name_type"].isin(g)]
    top5 = aux_df.sort_values("TPM_transcript",ascending=False).head(5)["gene_name"].values.tolist()
    import datetime
    now = datetime.datetime.now()
    t = now.time()
    print ("n_random: {}".format(t))
    aux_df.to_csv("{}{}_len{}_{}_list_genes.tsv".format(save_path,t,len(aux_df),type_),sep="\t")
    print ("Top5: {}\n".format(top5))
    return (g)

def plot_universe_npmi(query,title):
    fig,ax = plt.subplots(figsize=(8,8))
    x = contacts_with_features_df.query(query)
    x["npmi_max"] = x[["value.cast", "value.s129"]].max(axis=1)
    ax.plot([-4,4],[-4,4],color="black",linestyle='dashed')
    ax.plot([-4,3],[-3,4],color="black",linestyle='dashed')
    ax.plot([-3,4],[-4,3],color="black",linestyle='dashed')
    ax.axhline(1,color="black",linestyle='dashed')
    ax.axvline(1,color="black",linestyle='dashed')
    x = x.query('`zscore.dist.cast` >= 1 or `zscore.dist.s129` >= 1')
    prc = x.query('bin1_PRC > 0 and bin2_PRC > 0')
    sc = ax.scatter(x["zscore.dist.cast"],x["zscore.dist.s129"],s=5,c=x["npmi_max"],cmap="jet",vmax=0.7,vmin=0.3)
    ax.set_xlabel("CAST zscore")
    ax.set_ylabel("S129 zscore")
    ax.set_title(title)
    ax.set_xlim([-4,4])
    ax.set_ylim([-4,4])
    plt.colorbar(sc)


def plot_universe_npmi_axis(query,title,type_,add_prc=False):
    fig,ax = plt.subplots(figsize=(8,8))
    x = contacts_with_features_df.query(query)
    ax.plot([0.1,1],[0.1,1],color="black",linestyle='dashed')
    ax.plot([0.2,1],[0.1,0.9],color="black",linestyle='dashed')
    ax.plot([0.1,0.9],[0.2,1],color="black",linestyle='dashed')
    cast_specific = x.query('`zscore.dist.cast` >= 2 and delta >= 2')
    s129_specific = x.query('delta <= -2 and `zscore.dist.s129` >= 2')
    common = x.query('(`zscore.dist.cast` >= 2 and `zscore.dist.s129` >= 2)')
    cast_c = x.query('`zscore.dist.cast` >= 2')
    s129_c = x.query('`zscore.dist.s129` >= 2')
    
    if add_prc:
        prc = x.query('bin1_PRC > 0 and bin2_PRC > 0 and (`zscore.dist.cast` >= 2 or `zscore.dist.s129` >= 2) ')
        ax.scatter(prc["zscore.dist.cast"],prc["zscore.dist.s129"],s=10,color="red")
    ax.scatter(x["value.cast"],x["value.s129"],s=5,color="darkgrey",alpha=0.5)
    ax.scatter(cast_c["value.cast"],cast_c["value.s129"],s=8,color="darkgrey",alpha=0.5)
    ax.scatter(s129_c["value.cast"],s129_c["value.s129"],s=8,color="darkgrey",alpha=0.5)
    ax.scatter(cast_specific["value.cast"],cast_specific["value.s129"],s=8,color=cast_color)
    ax.scatter(s129_specific["value.cast"],s129_specific["value.s129"],s=8,color=s129_color)
    ax.scatter(common["value.cast"],common["value.s129"],s=8,color=non_phased_color)

    ax.set_xlabel("CAST npmi")
    ax.set_ylabel("S129 npmi")
    ax.set_title(title)
    ax.set_xlim([0.1,1])
    ax.set_ylim([0.1,1])

#font = {
#        'size'   : 22
#        }
#import matplotlib
#matplotlib.rc('font', **font)


def plot_universe(query,title,type_,save_name,zscore1=2,zscore2 = 2,delta=2,add_prc = False):
    fig,ax = plt.subplots(figsize=(8,8))
    x = contacts_with_features_df.query(query)
    ax.axhline(zscore2,color=s129_color,linestyle='dashed')
    ax.axvline(zscore2,color=cast_color,linestyle='dashed')
    ax.fill_between([zscore2,4], -4, 4, color=cast_color, alpha=0.1)
    ax.fill_between([-4,4], zscore2, 4, color=s129_color, alpha=0.1)
    cast_specific_string = '`zscore.dist.cast` >= {} and `zscore.dist.s129` < {}'.format(zscore1,zscore1)
    s129_specific_string = '`zscore.dist.s129` >= {} and `zscore.dist.cast` < {}'.format(zscore1,zscore1)
    common_strong_string = '`zscore.dist.s129` >= {} and `zscore.dist.cast` >= {}'.format(zscore1,zscore1)
    cast_specific = x.query(cast_specific_string)
    s129_specific = x.query(s129_specific_string)
    common_strong = x.query(common_strong_string)
    cast_c = x.query('`zscore.dist.cast` >= {}'.format(zscore2))
    s129_c = x.query('`zscore.dist.s129` >= {}'.format(zscore2))
    ax.scatter(x["zscore.dist.cast"],x["zscore.dist.s129"],s=5,color="darkgrey",alpha=0.5)
    ax.set_xlabel("CAST zscore")
    ax.set_ylabel("S129 zscore")
    ax.set_title(title)
    ax.set_xlim([-4,4])
    ax.set_ylim([-4,4])
    print ("CAST specific " + cast_specific_string)
    c = print_genes(cast_specific,type_)
    c = len(c)
    print ("CAST strong. zscore >= {}".format(zscore2))
    cc = print_genes(cast_c,type_)
    cc = len(cc)
    print ("S129 specific " + s129_specific_string)
    s = print_genes(s129_specific,type_)
    s = len(s)
    print ("S129 strong. zscore >= {}".format(zscore2))
    ss = print_genes(s129_c,type_)
    ss=len(ss)
    if add_prc:
        cast_prc_c = cast_c.query('bin1_PRC > 0 and bin2_PRC > 0')
        s129_prc_c = s129_c.query('bin1_PRC > 0 and bin2_PRC > 0 ')
        print ("Cast contacts. Cast genes with PRC")
        cpc = print_genes(cast_prc_c,type_)
        cpc = len(cpc)
        print ("S129 contacts. Cast genes with PRC")
        spc = print_genes(s129_prc_c,type_)
        spc = len(spc)
        ax.text(-3.7,2.8,"Genes in S129 strong contacts: {}      w/ PRC: {}".format(ss,spc),color=s129_color)
        ax.text(-0.3,-3.7,"Genes in CAST strong contacts: {}      w/ PRC: {}".format(cc,cpc),color=cast_color)
        cast_prc_c = cast_specific.query('bin1_PRC > 0 and bin2_PRC > 0 ')
        s129_prc_c = s129_specific.query('bin1_PRC > 0 and bin2_PRC > 0 ')
        print ("Cast specific. Cast genes with PRC")
        cpc = print_genes(cast_prc_c,type_)
        cpc = len(cpc)
        print ("S129 specific. Cast genes with PRC")
        spc = print_genes(s129_prc_c,type_)
        spc = len(spc)
        ax.text(-3.7,3,"Genes in S129 differential contacts: {}  w/ PRC {}".format(s,spc),color=s129_color)
        ax.text(-0.3,-3.5,"Genes in CAST differential contacts: {}  w/ PRC {}".format(c,cpc),color=cast_color)
    else:
        #ax.text(-3.7,3,"Genes in S129 differential contacts: {}".format(s),color=s129_color)
        #ax.text(-3.7,2.8,"Genes in S129 strong contacts: {}".format(ss),color=s129_color)
        #ax.text(0,-3.5,"Genes in CAST differential contacts: {}".format(c),color=cast_color)
        #ax.text(0,-3.7,"Genes in CAST strong contacts: {}".format(cc),color=cast_color)
        pass
    if type_ == "cast_contacts": #zscore >= 2
        y = contacts_with_features_df.query(query+ " and `zscore.dist.cast` >= {}".format(zscore2))
        ax.scatter(cast_c["zscore.dist.cast"],cast_c["zscore.dist.s129"],s=8,color=cast_color,alpha=1)
        write_to_pairs(y)
        write_to_file(y)
    elif type_ == "s129_contacts": #zscore >= 2
        y = contacts_with_features_df.query(query+ " and `zscore.dist.s129` >= {}".format(zscore2))
        ax.scatter(s129_c["zscore.dist.cast"],s129_c["zscore.dist.s129"],s=8,color=s129_color,alpha=1)
        write_to_pairs(y)
        write_to_file(y)
    elif type_ == "s129_specific": #zscore >= 1 and delta >= 1
        y = contacts_with_features_df.query(query+ " and "+ s129_specific_string)
        ax.scatter(s129_specific["zscore.dist.cast"],s129_specific["zscore.dist.s129"],s=8,color=s129_color)
        ax.scatter(cast_specific["zscore.dist.cast"],cast_specific["zscore.dist.s129"],s=8,color=s129_color)
        ax.scatter(common_strong["zscore.dist.cast"],common_strong["zscore.dist.s129"],s=8,color=s129_color)
        write_to_pairs(y)
        write_to_file(y)
    elif type_ == "cast_specific": #zscore >= 1 and delta >= 1
        y = contacts_with_features_df.query(query+ " and " + cast_specific_string)
        ax.scatter(cast_specific["zscore.dist.cast"],cast_specific["zscore.dist.s129"],s=8,color=cast_color)
        ax.scatter(s129_specific["zscore.dist.cast"],s129_specific["zscore.dist.s129"],s=8,color=cast_color)
        ax.scatter(common_strong["zscore.dist.cast"],common_strong["zscore.dist.s129"],s=8,color=cast_color)
        write_to_pairs(y)
        write_to_file(y)
    if add_prc:
        #promoter both = red
        prc = contacts_with_features_df.query(query+ " and bin1_PRC > 0 and bin2_PRC > 0 and (`zscore.dist.cast` >= {} or `zscore.dist.s129` >= {})".format(zscore2,zscore2))
        ax.scatter(prc["zscore.dist.cast"],prc["zscore.dist.s129"],s=10,color="red")
        print ("### Number of repressed genes contacting poised enhancer {}".format(len(prc)))
        aux_df = prc.query("`zscore.dist.cast` >= {} and `zscore.dist.s129` < {}".format(zscore1,zscore1))
        print ("CAST specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` < {}".format(zscore1,zscore1))
        print ("S129 specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` >= {}".format(zscore1,zscore1))
        print ("Common strong: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df.to_csv(save_path+"test1.tsv",sep="\t")

        #enhancer only PRC = blue
        prc1 = contacts_with_features_df.query(query+ " and ((bin1_PRC == 0 and bin1_CAST > 0) or (bin1_PRC == 0 and bin1_S129 > 0)) and bin2_PRC > 0 and (`zscore.dist.cast` >= {} or `zscore.dist.s129` >= {})".format(zscore2,zscore2))
        ax.scatter(prc1["zscore.dist.cast"],prc1["zscore.dist.s129"],s=10,alpha=0.8,color="blue")
        prc2 = contacts_with_features_df.query(query+ " and ((bin2_PRC == 0 and bin2_CAST > 0) or (bin2_PRC == 0 and bin2_S129 > 0)) and bin1_PRC > 0 and (`zscore.dist.cast` >= {} or `zscore.dist.s129` >= {})".format(zscore2,zscore2))
        ax.scatter(prc2["zscore.dist.cast"],prc2["zscore.dist.s129"],s=10,alpha=0.8,color="blue")
        prc = pd.concat([prc1,prc2])
        print ("### Number of genes contacting poised enhancer {}".format(len(prc)))
        aux_df = prc.query("`zscore.dist.cast` >= {} and `zscore.dist.s129` < {}".format(zscore1,zscore1))
        print ("CAST specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` < {}".format(zscore1,zscore1))
        print ("S129 specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` >= {}".format(zscore1,zscore1))
        print ("Common strong: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
       
        #promoter only PRC = pink
        prc1 = contacts_with_features_df.query(query+ " and ((bin1_PRC > 0 and bin1_CAST > 0) or (bin1_PRC > 0 and bin1_S129 > 0)) and bin2_PRC == 0 and (`zscore.dist.cast` >= {} or `zscore.dist.s129` >= {})".format(zscore2,zscore2))
        ax.scatter(prc1["zscore.dist.cast"],prc1["zscore.dist.s129"],s=7,alpha=0.8,color="fuchsia")
        prc2 = contacts_with_features_df.query(query+ " and ((bin2_PRC > 0 and bin2_CAST > 0) or (bin2_PRC > 0 and bin2_S129 > 0)) and bin1_PRC == 0 and (`zscore.dist.cast` >= {} or `zscore.dist.s129` >= {})".format(zscore2,zscore2))
        ax.scatter(prc2["zscore.dist.cast"],prc2["zscore.dist.s129"],s=7,alpha=0.8,color="fuchsia")
        prc = pd.concat([prc1,prc2])
        print ("### Number of repressed genes contacting enhancer: {}".format(len(prc)))
        aux_df = prc.query("`zscore.dist.cast` >= {} and `zscore.dist.s129` < {}".format(zscore1,zscore1))
        print ("CAST specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` < {}".format(zscore1,zscore1))
        print ("S129 specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` >= {}".format(zscore1,zscore1))
        print ("Common strong: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
 
        # no PRC
        no_prc = contacts_with_features_df.query(query+ " and bin1_PRC == 0 and bin2_PRC == 0 and (`zscore.dist.cast` >= {} or `zscore.dist.s129` >= {})".format(zscore2,zscore2))
        print ("### Number of contacts with no PRC {}".format(len(no_prc)))
        aux_df = no_prc.query("`zscore.dist.cast` >= {} and `zscore.dist.s129` < {}".format(zscore1,zscore1))
        print ("CAST specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = no_prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` < {}".format(zscore1,zscore1))
        print ("S129 specific: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df = no_prc.query("`zscore.dist.s129` >= {} and `zscore.dist.cast` >= {}".format(zscore1,zscore1))
        print ("Common strong: {}".format(len(aux_df)))
        print_genes(aux_df,type_,True)
        aux_df.to_csv(save_path+"test2.tsv",sep="\t")

        print ("### Total contacts:")
        print ("CAST specific:  {}".format(len(cast_specific)))
        print ("S129 specific:  {}".format(len(s129_specific)))
        print ("Common and strong:  {}".format(len(common_strong)))


    plt.savefig("{}/universe_{}.pdf".format(save_path,save_name),dpi=300)
    plt.savefig("{}/universe_{}.svg".format(save_path,save_name),dpi=300)
    plt.savefig("{}/universe_{}.png".format(save_path,save_name),dpi=300)

    return y,contacts_with_features_df.query(query)
    #texts = [ax.text(a, b, c, ha='center', va='center') for (a,b,c) in zip(cast_specific["zscore.dist.cast"],cast_specific["zscore.dist.s129"],cast_specific.)]
    #adjust_text(texts)
#%%
def add_gene_names(contacts_df,ase_df,type_):
    header = contacts_df.columns.values.tolist()
    cont_pb = pb.BedTool.from_dataframe(contacts_df)
    ase_df = ase_df[["chrom","start_transcript","end_transcript","gene_name_type"]]
    ase_pb = pb.BedTool.from_dataframe(ase_df)
    anchor1_df = cont_pb.intersect(ase_pb,wa=True,wb=True).to_dataframe(names = header+["a","b","c","bin1_{}_{}".format(type_,"gene_names")])
    c_merged = pd.merge(contacts_df,anchor1_df,on=header,how="outer")
    new_header = [item for sublist in [header[3:6],header[0:3]] for item in sublist]
    aux_header = anchor1_df.columns.values.tolist()
    new_header = new_header+aux_header[6:]
    anchor1_df=c_merged[new_header]
    anchor1_pb = pb.BedTool.from_dataframe(anchor1_df)
    anchor1_2_df = anchor1_pb.intersect(ase_pb,wa=True,wb=True).to_dataframe(names = new_header+["d","e","f","bin2_{}_{}".format(type_,"gene_names")])
    anchor1_df = anchor1_df.replace('.', np.NaN)
    anchor1_2_df = anchor1_2_df.replace('.', np.NaN)
    anchor1_df['b'] = anchor1_df['b'].astype('float64')
    anchor1_df['c'] = anchor1_df['c'].astype('float64')
    anchor1_2_df['b'] = anchor1_2_df['b'].astype('float64')
    anchor1_2_df['c'] = anchor1_2_df['c'].astype('float64')
    anchor1_2_df = pd.merge(anchor1_df,anchor1_2_df,on=new_header, how="outer")
    new_header = [item for sublist in [new_header[3:6],new_header[0:3]] for item in sublist]
    aux_header = anchor1_2_df.columns.values.tolist()
    new_header = new_header+aux_header[6:]
    anchor1_2_df=anchor1_2_df[new_header]
    anchor1_2_df.drop(columns=["a","b","c","d","e","f"],inplace=True)
    anchor1_2_df = anchor1_2_df.drop_duplicates()
    print ("{} added.".format("gene names"))
    return anchor1_2_df

def add_feature (feature_name, feature_path_or_df,contacts_df):
    header = contacts_df.columns.values.tolist()
    cont_pb = pb.BedTool.from_dataframe(contacts_df)
    if isinstance(feature_path_or_df, str):
        peaks_pb = pb.BedTool(feature_path_or_df)
    else:
        peaks_pb = pb.BedTool.from_dataframe(feature_path_or_df)

    anchor1_df = cont_pb.intersect(peaks_pb,wa=True,c=True).to_dataframe(names = header+["bin1_{}".format(feature_name)])
    new_header = [item for sublist in [header[3:6],header[0:3]] for item in sublist]
    aux_header = anchor1_df.columns.values.tolist()
    new_header = new_header+aux_header[6:]
    anchor1_df=anchor1_df[new_header]
    anchor1_pb = pb.BedTool.from_dataframe(anchor1_df)
    anchor1_2_df = anchor1_pb.intersect(peaks_pb,wa=True,c=True).to_dataframe(names = new_header+["bin2_{}".format(feature_name)])
    new_header = [item for sublist in [new_header[3:6],new_header[0:3]] for item in sublist]
    aux_header = anchor1_2_df.columns.values.tolist()
    new_header = new_header+aux_header[6:]
    anchor1_2_df=anchor1_2_df[new_header]
    print ("{} added.".format(feature_name))
    return anchor1_2_df

def plot_manhattan(df,text,color):
    fig,axs = plt.subplots(19,1,figsize=(40,20),sharex=True,sharey=True)
    for i,chr_ in enumerate(chrs_):
        ax = axs[i]
        aux_ = df.query('chrom1 == "{}"'.format(chr_))
        values = aux_.start1.values.tolist()+aux_.start2.values.tolist() 
        values_dict = Counter(values)
        w = 50000
        ax.bar(values_dict.keys(),values_dict.values(),width=w,color=color)
        aux1 = cast_df.query('chrom == "{}"'.format(chr_))
        ax.scatter(aux1.start_transcript,len(aux1)*[-1],color=cast_color,s=1)
        aux3 = exp_df.query('chrom == "{}"'.format(chr_))
        ax.scatter(aux3.start_transcript,len(aux3)*[-1.5],color="darkgreen",s=1)
    fig.savefig("{}/manhattan_plot_n_diff_contacts_in_{}.pdf".format(save_path,text))

def plot_manhattan_zoom(df,chrom, start, end, text,color):
    fig,ax = plt.subplots(1,1,figsize=(40,20))
    aux_ = df.query('chrom1 == "{}" and start1 >= {} and start2 >= {} and end1 <= {} and end2 <= {}'.format(chrom,start,start,end,end))
    values = aux_.start1.values.tolist()+aux_.start2.values.tolist() 
    values_dict = Counter(values)
    w = 50000
    ax.bar(values_dict.keys(),values_dict.values(),width=w,color=color)
    aux1 = cast_df.query('chrom == "{}" and start_transcript >= {} and end_transcript <= {}'.format(chrom,start,end))
    ax.scatter(aux1.start_transcript,len(aux1)*[-1],color=cast_color)
    aux2 = s129_df.query('chrom == "{}" and start_transcript >= {} and end_transcript <= {}'.format(chrom,start,end))
    ax.scatter(aux2.start_transcript,len(aux2)*[-1.5],color=s129_color)
    aux3 = exp_df.query('chrom == "{}" and start_transcript >= {} and end_transcript <= {}'.format(chrom,start,end))
    ax.scatter(aux3.start_transcript,len(aux3)*[-2],color="darkgreen")
    fig.savefig("{}/manhattan_plot_n_diff_contacts_in_{}.pdf".format(save_path,text))

def write_to_pairs(df,name="-"):
    df[["chrom1","start1","end1","chrom2","start2","end2"]].to_csv(save_path+"/dynamic_feature.arcs",sep="\t",index=False,header=None)
    if name != "-":
        df[["chrom1","start1","end1","chrom2","start2","end2"]].to_csv(save_path+"/{}.arcs".format(name),sep="\t",index=False,header=None)
def write_to_file(df,name="-"):
    df.to_csv(save_path+"/dynamic_feature.bed",sep="\t",index=False)
    if name != "-":
        df.to_csv(save_path+"/{}.bed".format(name),sep="\t",index=False)

def plot_volcano(master_df,name):
    import datetime
    date_ = datetime.datetime.now()
    date = "{}-{}-{}".format(date_.day,date_.month,date_.year)
    plot_all_text = False
    plot_name = "TPM >= 1"
    figname = "volcano_{}_{}.svg".format(name,date)
    ase_df = master_df.query('TPM_transcript >= 1')
    cast_df = ase_df.query('log2foldchange >= 1 and p_adj_sum <= 0.05')
    s129_df = ase_df.query('log2foldchange <= -1 and p_adj_sum <= 0.05')
    print ("Cast genes:{}".format(len(cast_df)))
    print ("S129 genes:{}".format(len(s129_df)))
    rest_df = pd.concat([ase_df, cast_df, s129_df]).drop_duplicates(keep=False)

    fig,ax = plt.subplots(1,1,figsize=(8,10))
    size = 12
    x = rest_df["log2foldchange"]
    y = np.log10(rest_df["p_adj_sum"])
    ax.scatter(x,-1*y,color="grey",s=size)

    x = cast_df["log2foldchange"]
    y = np.log10(cast_df["p_adj_sum"])
    ax.scatter(x,-1*y,color=cast_color,s=size)
    examples = cast_df.query('(log2foldchange > 6.5  or p_adj_sum < 1E-200) and log2foldchange != inf and p_adj_sum != 0 and  transcript_biotype == "protein_coding"' )
    x = examples.log2foldchange.values.tolist() 
    y = -1*np.log10(examples.p_adj_sum.values.tolist())
    names = examples.gene_name.values.tolist()
    texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)
    print (len(cast_df))
    print ("CAST coding")
    print (len(cast_df.query('transcript_biotype == "protein_coding"')))
    print ("CAST lncRNA")
    print (len(cast_df.query('transcript_biotype == "lncRNA"')))

    x = s129_df["log2foldchange"]
    y = np.log10(s129_df["p_adj_sum"])
    ax.scatter(x,-1*y,color=s129_color,s=size)
    examples = s129_df.query('(log2foldchange < -6  or p_adj_sum < 1E-180) and log2foldchange != -inf and p_adj_sum != 0' )
    x = examples.log2foldchange.values.tolist() 
    y = -1*np.log10(examples.p_adj_sum.values.tolist())
    names = examples.gene_name.values.tolist()
    texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)
    print (len(s129_df))
    print ("S129 coding")
    print (len(s129_df.query('transcript_biotype == "protein_coding"')))
    print ("S129 lncRNA")
    print (len(s129_df.query('transcript_biotype == "lncRNA"')))

    cast_mono_df = cast_df.query('log2foldchange == inf')
    x = [12 for i in range(len(cast_mono_df))]
    y = np.log10(cast_mono_df["p_adj_sum"])
    ax.scatter(x,-1*y,color=cast_color,s=size*2,marker='>')
    examples = cast_mono_df.query('p_adj_sum != 0 and p_adj_sum < 1E-40 and transcript_biotype == "protein_coding"')
    x = [12 for i in range(len(examples))]
    y = -1*np.log10(examples.p_adj_sum.values.tolist())
    names = examples.gene_name.values.tolist()
    texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)

    cast_0_pval = cast_df.query('p_adj_sum == 0')
    x = cast_0_pval["log2foldchange"]
    y = [330 for i in range(len(cast_0_pval))]
    ax.scatter(x,y,color=cast_color,s=size*2,marker='^')
    examples = cast_0_pval.query('transcript_biotype == "protein_coding"')
    x = examples.log2foldchange.values.tolist() 
    y = [330 for i in range(len(examples))]
    names = examples.gene_name.values.tolist()
    texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)

    s129_mono_df = s129_df.query('log2foldchange == -inf')
    x = [-12 for i in range(len(s129_mono_df))]
    y = np.log10(s129_mono_df["p_adj_sum"])
    ax.scatter(x,-1*y,color=s129_color,s=size*2,marker='<')
    examples = s129_mono_df.query('p_adj_sum != 0 and p_adj_sum < 1E-8 and transcript_biotype == "protein_coding"')
    x = [-12 for i in range(len(examples))]
    y = -1*np.log10(examples.p_adj_sum.values.tolist())
    names = examples.gene_name.values.tolist()
    texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)

    s129_0_pval = s129_df.query('p_adj_sum == 0')
    x = s129_0_pval["log2foldchange"]
    y = [330 for i in range(len(s129_0_pval))]
    ax.scatter(x,y,color=s129_color,s=size*2,marker='^')
    examples = s129_0_pval.query('log2foldchange != -inf and transcript_biotype == "protein_coding"')
    x = examples.log2foldchange.values.tolist() 
    y = [330 for i in range(len(examples))]
    names = examples.gene_name.values.tolist()
    texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)
    if plot_all_text:
        examples = master_df.query('(log2foldchange <= -1 or log2foldchange > 1) and p_adj_sum <= 0.05')
        x = examples.log2foldchange.values.tolist() 
        y = -1*np.log10(examples.p_adj_sum.values.tolist())
        names = examples.gene_name.values.tolist()
        texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
        adjust_text(texts)

    s129_corner = s129_df.query('p_adj_sum == 0 and log2foldchange == -inf')
    y = [330 for i in range(len(s129_corner))]
    x = [-12 for i in range(len(s129_corner))]
    ax.scatter(x,y,color=s129_color,s=size*2,marker='<')
    names = s129_corner.gene_name.values.tolist()
    texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
    adjust_text(texts)

    p_adj_cutoff = (-1*np.log10(0.05))
    ax.axhline(p_adj_cutoff,color="black",linestyle="--")
    #ax.axhline(1,color="black",linestyle="--")
    ax.axvline(-1,color="black",linestyle="--")
    ax.axvline(1,color="black",linestyle="--")
    ax.set_ylabel("$-log_{10}$(P-value)")
    ax.set_xlabel("$log_2(FoldChange)$")
    ax.set_title(plot_name)
    ax.set_ylim([-10,350])
    ax.set_xlim([-12,12])
    plt.savefig("{}{}".format(save_path,figname),dpi=300)

#%%

#get WDF of the contacts
seg_cast_path = root+"/Curated_F123_1122_ALL_as_3NP_CAST_at50000.passed_qc_fc5_cw6_s11.table"
seg_s129_path = root+"/Curated_F123_1122_ALL_as_3NP_S129_at50000.passed_qc_fc5_cw6_s11.table"

def get_WDF(x):
    return x.sum()/n_NPs

seg_table_cast = pd.read_csv(seg_cast_path,sep="\t")
seg_table_cast.set_index(['chrom','start','stop'], inplace = True)
seg_table_cast = seg_table_cast.loc[chrs_]
n_NPs = seg_table_cast.shape[1]
seg_table_cast["cast_WDF"] = seg_table_cast.apply(lambda x: get_WDF(x),axis=1)
print (seg_table_cast)

plt.hist(seg_table_cast.cast_WDF,bins=20)

seg_table_s129 = pd.read_csv(seg_s129_path,sep="\t")
seg_table_s129.set_index(['chrom','start','stop'], inplace = True)
seg_table_s129 = seg_table_s129.loc[chrs_]
n_NPs = seg_table_s129.shape[1]
seg_table_s129["s129_WDF"] = seg_table_s129.apply(lambda x: get_WDF(x),axis=1)
seg_table_cast = seg_table_cast["cast_WDF"].reset_index()
seg_table_s129 = seg_table_s129["s129_WDF"].reset_index()
wdf_cast_pb = pb.BedTool.from_dataframe(seg_table_cast).sort()
wdf_s129_pb = pb.BedTool.from_dataframe(seg_table_s129).sort()
plt.hist(seg_table_s129.s129_WDF,bins=20)

#%%

def get_WDF_of_contacts(df,only_promoters = False,cast_promoters=True,add_bins = 0): #add bins means to the left and right of the promoter
    slope = 50000*add_bins
    if not only_promoters:
        wdf_values_cast = []
        #first get for CAST
        aux_pb = pb.BedTool.from_dataframe(df[["chrom1","start1","end1"]])
        cast_df = aux_pb.intersect(wdf_cast_pb,wb=True).to_dataframe(names=["chrom1","start1","end1","cast_WDF"])
        wdf_values_cast.extend(cast_df.cast_WDF.values.tolist())
        aux_pb = pb.BedTool.from_dataframe(df[["chrom2","start2","end2"]])
        cast_df = aux_pb.intersect(wdf_cast_pb,wb=True).to_dataframe(names=["chrom2","start2","end2","cast_WDF"])
        wdf_values_cast.extend(cast_df.cast_WDF.values.tolist())

        wdf_values_s129 = []
        #then S129
        aux_pb = pb.BedTool.from_dataframe(df[["chrom1","start1","end1"]])
        s129_df = aux_pb.intersect(wdf_s129_pb,wb=True).to_dataframe(names=["chrom1","start1","end1","s129_WDF"])
        wdf_values_s129.extend(s129_df.s129_WDF.values.tolist())
        aux_pb = pb.BedTool.from_dataframe(df[["chrom2","start2","end2"]])
        s129_df = aux_pb.intersect(wdf_s129_pb,wb=True).to_dataframe(names=["chrom2","start2","end2","s129_WDF"])
        wdf_values_s129.extend(s129_df.s129_WDF.values.tolist())
        return wdf_values_cast,wdf_values_s129
    if only_promoters:
        if cast_promoters:
            s1 = "bin1_CAST > 0"
            s2 = "bin2_CAST > 0"
        else:
            s1 = "bin1_S129 > 0"
            s2 = "bin2_S129 > 0"
        #first get for CAST
        wdf_values_cast = []
        df1 = df.query(s1)
        df1["start1"] = df1["start1"]-slope
        df1["end1"] = df1["end1"]+slope
        aux_pb = pb.BedTool.from_dataframe(df1[["chrom1","start1","end1"]])
        cast_df = aux_pb.intersect(wdf_cast_pb,wb=True).to_dataframe(names=["chrom1","start1","end1","cast_WDF"])
        wdf_values_cast.extend(cast_df.cast_WDF.values.tolist())
        df2 = df.query(s2)
        df2["start2"] = df2["start2"]-slope
        df2["end2"] = df2["end2"]+slope
        aux_pb = pb.BedTool.from_dataframe(df2[["chrom2","start2","end2"]])
        cast_df = aux_pb.intersect(wdf_cast_pb,wb=True).to_dataframe(names=["chrom2","start2","end2","cast_WDF"])
        wdf_values_cast.extend(cast_df.cast_WDF.values.tolist())

        #then S129
        wdf_values_s129 = []
        df1 = df.query(s1)
        df1["start1"] = df1["start1"]-slope
        df1["end1"] = df1["end1"]+slope
        aux_pb = pb.BedTool.from_dataframe(df1[["chrom1","start1","end1"]])
        s129_df = aux_pb.intersect(wdf_s129_pb,wb=True).to_dataframe(names=["chrom1","start1","end1","s129_WDF"])
        wdf_values_s129.extend(s129_df.s129_WDF.values.tolist())
        df2 = df.query(s2)
        df2["start2"] = df2["start2"]-slope
        df2["end2"] = df2["end2"]+slope
        aux_pb = pb.BedTool.from_dataframe(df2[["chrom2","start2","end2"]])
        s129_df = aux_pb.intersect(wdf_s129_pb,wb=True).to_dataframe(names=["chrom2","start2","end2","s129_WDF"])
        wdf_values_s129.extend(s129_df.s129_WDF.values.tolist())
        return wdf_values_cast,wdf_values_s129


def wdf_boxplots(plot_values_cast,plot_values_s129,title):
    fig,ax1 = plt.subplots(1,1,figsize=(4,8))
    bp1 = ax1.boxplot([plot_values_cast,plot_values_s129],patch_artist=True ,showfliers=False,widths=0.5)
    ax1.set_xticklabels(["cast","s129"])
    colors = [cast_color,s129_color]
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set(alpha=0.5)
    for pos,v in enumerate([plot_values_cast,plot_values_s129]): 
        xv = np.random.normal(pos+1, 0.05, len(v))
        ax1.scatter(xv, v,s=1,alpha=1,color="black")
    plt.setp(bp1["medians"],color="white")
    ax1.set_ylabel("WDF")
    ax1.set_title(title)
    ax1.set_ylim([0.02,0.15])
    plt.savefig("{}/avg_WDF_in_{}.pdf".format(save_path,title),dpi=300)

    fig,ax1 = plt.subplots(1,1,figsize=(2,8))
    plot_values = [x-y for (x,y) in zip(plot_values_cast,plot_values_s129)]
    bp1 = ax1.boxplot(plot_values,patch_artist=True ,showfliers=False,widths=0.5)
    ax1.set_xticklabels(["CAST-S129"])
    ax1.axhline(0)
    colors = ["blue"]
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set(alpha=0.5)
    for pos,v in enumerate([plot_values]): 
        xv = np.random.normal(pos+1, 0.05, len(v))
        ax1.scatter(xv, v,s=1,alpha=1,color="black")
    plt.setp(bp1["medians"],color="white")
    ax1.set_ylabel("Differencia WDF (CAST-S129)")
    ax1.set_title(title)
    ax1.set_ylim([-0.04,0.04])
    plt.savefig("{}/Differential_WDF_in_{}.pdf".format(save_path,title),dpi=300)

#%%

def stacked_bar_plot(ax,type_,x,color):
    ax.bar(type_,x, color=color,label="activation")

def stacked_prc_bar_plot(ax,type_,x,y,z=0):
    bottom = 0 
    ax.bar(type_,x,bottom = bottom, color="darkgreen",label="activation")
    #ax.bar(type_,y,bottom=x,color="crimson",label="PRC in both")
    ax.bar(type_,y,bottom=x,color="crimson",label="PRC in gene")
    summup = x + y
    ax.bar(type_,z,bottom=summup,color="pink",label="PRC in one")

def do_horizontal_bar_prc_plot_normalized(name,cast_all,s129_all,common_all,cast,s129,common):
    print ("-------")
    print (cast_all,s129_all,common_all)
    fig,ax = plt.subplots(figsize=(4,6))
    stacked_prc_bar_plot(ax,"CAST contact",cast[0]/cast_all,cast[1]/cast_all,cast[2])
    stacked_prc_bar_plot(ax,"S129 contact",s129[0]/s129_all,s129[1]/s129_all,s129[2])
    stacked_prc_bar_plot(ax,"Both contact",common[0]/common_all,common[1]/common_all,common[2])
    ax.set_ylabel("Percentage")
    plt.legend(loc="upper right")
    plt.title(name)
    ax.set_ylim([0,1.05])
    fig.savefig("{}{}_normalized.pdf".format(save_path,name),dpi=300)

def do_horizontal_bar_prc_plot(name,cast_all,s129_all,common_all,cast,s129,common):
    print (cast_all,s129_all,common_all)
    fig,ax = plt.subplots(figsize=(4,6))
    stacked_prc_bar_plot(ax,"CAST contact",cast[0],cast[1],cast[2])
    stacked_prc_bar_plot(ax,"S129 contact",s129[0],s129[1],s129[2])
    stacked_prc_bar_plot(ax,"Both contact",common[0],common[1],common[2])
    ax.set_ylabel("Percentage")
    plt.legend(loc="upper right")
    plt.title(name)
    ax.set_ylim([0,140])
    fig.savefig("{}{}.pdf".format(save_path,name),dpi=300)

def do_horizontal_bar_plot(name,cast_all,s129_all,common_all,cast,s129,common):
    print (cast_all,s129_all,common_all)
    fig,ax = plt.subplots(figsize=(4,6))
    stacked_bar_plot(ax,"CAST contact",cast[0],cast_color)
    stacked_bar_plot(ax,"S129 contact",s129[0],s129_color)
    stacked_bar_plot(ax,"Both contact",common[0],both_color)
    ax.set_ylabel("Percentage")
    plt.legend(loc="upper right")
    plt.title(name)
    ax.set_ylim([0,140])
    fig.savefig("{}{}.pdf".format(save_path,name),dpi=300)


def get_gene_promoter_states(list_genes):
    aux_df = ase_df[ase_df["gene_name_type"].isin(list_genes)]
    return dict((x,aux_df.Promoter_state.values.tolist().count(x)) for x in set(aux_df.Promoter_state))
#    return set(aux_df.Promoter_state)
#%%
if recalculate_all:
    #columns: chr1	3000000	3030000	chr1	3330000	3360000	0.5141005133296008
    #cast_contacts_path = root+"/diff_contacts/specific_contacts_cast_NPMI03.pairs"  #OLD
    #cast_c = pd.read_csv(cast_contacts_path,sep="\t",names=["chrom1","start1","end1","chrom2","start2","end2","score"])
    #cast_contacts_path = root+"/diff_contacts/220629/contacts_zdelta34_CAST.zselected_dist.exclusive.1.2.4.10.15.20.50mb.tsv.gz"  #NEW
    #cast_contacts_path = root+"/diff_contacts/220629/contacts_zdelta4_CAST.extremes_dist.exclusive.1.2.4.10.15.20.50mb.tsv.gz"  #NEW
    all_contacts_path = root+"/contacts_full_annotated.zdelta_annotated.dist.exclusive.1.2.4.10.15.20.50mb.tsv.gz" #new 
    cont = pd.read_csv(all_contacts_path,sep="\t")
    cont = cont.rename(columns={"chrom_x":"chrom1","start_x":"start1","end_x":"end1","chrom_y":"chrom2","start_y":"start2","end_y":"end2"})

    peaks = root+"F123.S7p.BCP_peaks.HM_mode.bed"
    contacts_with_features_df = add_feature("S7p",peaks,cont)
    peaks = root+"F123.S5p.BCP_peaks.HM_mode.bed"
    contacts_with_features_df = add_feature("S5p",peaks,contacts_with_features_df)
    peaks = root+"F123.K27me3.only_peaks.bed"
    contacts_with_features_df = add_feature("PRC",peaks,contacts_with_features_df)

    peaks = root+"F123_ATAC_Bing_shared_peaks.bed"
    contacts_with_features_df = add_feature("ATAC_common",peaks,contacts_with_features_df)
    peaks = root+"F123_bing_data030420_genome_CAST.bed"
    contacts_with_features_df = add_feature("ATAC_CAST",peaks,contacts_with_features_df)
    peaks = root+"F123_bing_data030420_genome_S129.bed"
    contacts_with_features_df = add_feature("ATAC_S129",peaks,contacts_with_features_df)

    peaks = root+"H3K4me3.common.peaks.bed"
    contacts_with_features_df = add_feature("H3K4me3_common",peaks,contacts_with_features_df)
    peaks = root+"H3K4me3.CAST.specific.peaks.bed"
    contacts_with_features_df = add_feature("H3K4me3_CAST",peaks,contacts_with_features_df)
    peaks = root+"H3K4me3.S129.specific.peaks.bed"
    contacts_with_features_df = add_feature("H3K4me3_S129",peaks,contacts_with_features_df)

    peaks = root+"H3K27Ac.common.peaks.bed"
    contacts_with_features_df = add_feature("H3K27Ac_common",peaks,contacts_with_features_df)
    peaks = root+"H3K27Ac.CAST.specific.peaks.bed"
    contacts_with_features_df = add_feature("H3K27Ac_CAST",peaks,contacts_with_features_df)
    peaks = root+"H3K27Ac.S129.specific.peaks.bed"
    contacts_with_features_df = add_feature("H3K27Ac_S129",peaks,contacts_with_features_df)

    peaks = root+"CTCF.common.peaks.bed"
    contacts_with_features_df = add_feature("CTCF_common",peaks,contacts_with_features_df)
    peaks = root+"CTCF.CAST.specific.peaks.bed"
    contacts_with_features_df = add_feature("CTCF_CAST",peaks,contacts_with_features_df)
    peaks = root+"CTCF.S129.specific.peaks.bed"
    contacts_with_features_df = add_feature("CTCF_S129",peaks,contacts_with_features_df)

    peaks = root+"Rad21.common.peaks.bed"
    contacts_with_features_df = add_feature("Rad21_common",peaks,contacts_with_features_df)
    peaks = root+"Rad21.CAST.specific.peaks.bed"
    contacts_with_features_df = add_feature("Rad21_CAST",peaks,contacts_with_features_df)
    peaks = root+"Rad21.S129.specific.peaks.bed"
    contacts_with_features_df = add_feature("Rad21_S129",peaks,contacts_with_features_df)
    contacts_with_features_df = add_feature("exp_with_SNP",exp_df,contacts_with_features_df)
    contacts_with_features_df = add_feature("CAST",cast_df,contacts_with_features_df)
    contacts_with_features_df = add_feature("S129",s129_df,contacts_with_features_df)
    contacts_with_features_df = add_feature("CAST_PRC_promoter",cast_prc_df,contacts_with_features_df)
    contacts_with_features_df = add_feature("S129_PRC_promoter",s129_prc_df,contacts_with_features_df)

    aux_cast_df = ase_df.query("(log2foldchange >= 1) and p_adj_sum <= 0.05 and TPM_transcript >= 1")
    aux_s129_df = ase_df.query("(log2foldchange <= -1) and p_adj_sum <= 0.05 and TPM_transcript >= 1")

    #aux_names = add_gene_names(contacts_with_features_df,aux_ase_df)
    aux_names_cast = add_gene_names(contacts_with_features_df,aux_cast_df,"cast")
    header1 = ["chrom1","start1","end1","chrom2","start2","end2"]
    names1 = aux_names_cast.groupby(header1)['bin1_cast_gene_names'].apply(list)
    n1 = names1.reset_index()
    r1 = pd.merge(contacts_with_features_df,n1,on=header1,how="outer")
    header2 = ["chrom2","start2","end2","chrom1","start1","end1"]
    names2 = aux_names_cast.groupby(header2)['bin2_cast_gene_names'].apply(list)
    n2 = names2.reset_index()
    r2 = pd.merge(r1,n2,on=header2,how="outer")
    r2['bin1_cast_gene_names'] = r2['bin1_cast_gene_names'].apply(lambda x: list(set(x)) if isinstance(x, list) else x)
    r2['bin2_cast_gene_names'] = r2['bin2_cast_gene_names'].apply(lambda x: list(set(x)) if isinstance(x, list) else x)
    aux_names_s129 = add_gene_names(contacts_with_features_df,aux_s129_df,"s129")
    header1 = ["chrom1","start1","end1","chrom2","start2","end2"]
    names1 = aux_names_s129.groupby(header1)['bin1_s129_gene_names'].apply(list)
    n1 = names1.reset_index()
    r1 = pd.merge(contacts_with_features_df,n1,on=header1,how="outer")
    header2 = ["chrom2","start2","end2","chrom1","start1","end1"]
    names2 = aux_names_s129.groupby(header2)['bin2_s129_gene_names'].apply(list)
    n2 = names2.reset_index()
    r3 = pd.merge(r1,n2,on=header2,how="outer")
    r3['bin1_s129_gene_names'] = r3['bin1_s129_gene_names'].apply(lambda x: list(set(x)) if isinstance(x, list) else x)
    r3['bin2_s129_gene_names'] = r3['bin2_s129_gene_names'].apply(lambda x: list(set(x)) if isinstance(x, list) else x)

    def process_list(lst):
        if isinstance(lst, list):
            cleaned_lst = [str(x) for x in lst if pd.notnull(x)]
            return ', '.join(cleaned_lst)
        elif np.isnan(lst):
            return np.nan
        else:
            return lst
        
    contacts_with_features_df['bin1_cast_gene_names'] = r2['bin1_cast_gene_names'].apply(lambda x: process_list(x))
    contacts_with_features_df['bin2_cast_gene_names'] = r2['bin2_cast_gene_names'].apply(lambda x: process_list(x))
    contacts_with_features_df['bin1_s129_gene_names'] = r3['bin1_s129_gene_names'].apply(lambda x: process_list(x))
    contacts_with_features_df['bin2_s129_gene_names'] = r3['bin2_s129_gene_names'].apply(lambda x: process_list(x))
    contacts_with_features_df.to_csv(root+"full_contacts_with_features_and_genes.tsv",sep="\t",header=True,index=False)

def get_genes_in_contacts(df,type_,title):
    if type_ == "cast":
        list1_gene_names = df.bin1_cast_gene_names.values.tolist()
        list2_gene_names = df.bin2_cast_gene_names.values.tolist()
    else:
        list1_gene_names = df.bin1_s129_gene_names.values.tolist()
        list2_gene_names = df.bin2_s129_gene_names.values.tolist()
    list_gene_names = list1_gene_names+list2_gene_names
    set_list_gene_names = list(set(list_gene_names))
    set_list_gene_names.remove(".")
    cleaned_list_gene_names = [x.split("_")[0] for x in set_list_gene_names]
    with open('{}{}'.format(save_path,title), 'w') as outfile:
        outfile.write('\n'.join(str(i) for i in cleaned_list_gene_names))
    subset_df = ase_df[ase_df.gene_name.isin(cleaned_list_gene_names)]
    subset_df.to_csv("{}{}.tsv".format(save_path,title),sep="\t")
    return cleaned_list_gene_names

#%%

distance = 2000000

# Now first get CAST EP contacts 

# First CAST promoters

s = "bin1_S5p > 0 and bin2_S5p > 0 and (bin1_ATAC_CAST > 0 or bin1_ATAC_common > 0) and (bin2_ATAC_CAST > 0 or bin2_ATAC_common > 0) and (((bin1_H3K27Ac_CAST > 0 or bin1_H3K27Ac_common > 0) and bin2_CAST > 0) or ((bin2_H3K27Ac_CAST > 0 or bin2_H3K27Ac_common > 0) and bin1_CAST > 0)) and dist <= {}".format(distance)
x, all_ = plot_universe(s,"CAST E-P with common 2Mb","cast_specific","cast_EP",2,2,0,False)
list_names = print_genes(x,"cast_specific")
plot_volcano(ase_df[ase_df["gene_name_type"].isin(list_names)],"CAST_specific_EP")
plot_manhattan(x,"EP_CAST_CAST_contacts",cast_color)

# now get different subsets of contacts in files
x = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` < 2)") # CAST strong and S129 low
print ("Number of contacts like this:")
write_to_pairs(x,"CAST_strong_S129_weak_CAST_EP")
write_to_file(x,"CAST_strong_S129_weak_CAST_EP")
x = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2)") # CAST strong
write_to_pairs(x,"CAST_strong_CAST_EP")
write_to_file(x,"CAST_strong_CAST_EP")
x = contacts_with_features_df.query(s + "and (`zscore.dist.cast` < 2 and `zscore.dist.s129` >= 2)") # CAST weak and S129 strong
write_to_pairs(x,"CAST_weak_S129_strong_CAST_EP")
write_to_file(x,"CAST_weak_S129_strong_CAST_EP")

cast_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` < 2)") # CAST strong and S129 weak
s129_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` < 2 and `zscore.dist.s129` >= 2)") # S129 strong and CAST weak
common_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` >= 2)") # strong in both
#ll = get_genes_in_contacts(cast_c_total,"cast","gene_names_cast_cast_contacts.txt")
#ll = get_genes_in_contacts(s129_c_total,"cast","gene_names_cast_s129_contacts.txt")
#ll = get_genes_in_contacts(common_c_total,"cast","gene_names_cast_common_contacts.txt")
print ("Number of genes involved in cast_contacts: {}".format(len(print_genes(cast_c_total,"cast_specific"))))
print ("Number of genes involved in s129_contacts: {}".format(len(print_genes(s129_c_total,"cast_specific"))))
print ("Number of genes involved in common_contacts: {}".format(len(print_genes(common_c_total,"cast_specific"))))
do_horizontal_bar_plot("CAST_gene_mechanism",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total),0,0],[len(s129_c_total),0,0],[len(common_c_total),0,0]) 

# check if PRC is in the promoters
cast_c_prc_promoter = cast_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter > 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter > 0))") 
s129_c_prc_promoter = s129_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter > 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter > 0))") 
common_c_prc_promoter = common_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter > 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter > 0))")
cast_c_no_prc_promoter = cast_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter == 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter == 0))") 
s129_c_no_prc_promoter = s129_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter == 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter == 0))") 
common_c_no_prc_promoter = common_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter == 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter == 0))")
print ("Numbers:")
print (len(cast_c_total), len(cast_c_prc_promoter),len(cast_c_no_prc_promoter))
print (len(s129_c_total), len(s129_c_prc_promoter),len(s129_c_no_prc_promoter))
print (len(common_c_total), len(common_c_prc_promoter),len(common_c_no_prc_promoter))

do_horizontal_bar_prc_plot("CAST_gene_mechanism",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total)-len(cast_c_prc_promoter),len(cast_c_prc_promoter),0],[len(s129_c_total)-len(s129_c_prc_promoter),len(s129_c_prc_promoter),0],[len(common_c_total)-len(common_c_prc_promoter),len(common_c_prc_promoter),0]) #each bar total and then each list means: no prc, prc both, prc in one of the anchors
do_horizontal_bar_prc_plot_normalized("CAST_gene_mechanism",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total)-len(cast_c_prc_promoter),len(cast_c_prc_promoter),0],[len(s129_c_total)-len(s129_c_prc_promoter),len(s129_c_prc_promoter),0],[len(common_c_total)-len(common_c_prc_promoter),len(common_c_prc_promoter),0]) #each bar total and then each list means: no prc, prc both, prc in one of the anchors

#get WDF of all contacts
c,s = get_WDF_of_contacts(cast_c_total,only_promoters=True,cast_promoters=True,add_bins=2)
wdf_boxplots(c,s,"WDF_cast_contacts_cast_promoter_slope2")
c,s = get_WDF_of_contacts(s129_c_total,only_promoters=True,cast_promoters=True,add_bins=2)
wdf_boxplots(c,s,"WDF_s129_contacts_cast_promoter_slope2")
c,s = get_WDF_of_contacts(common_c_total,only_promoters=True,cast_promoters=True,add_bins=2)
wdf_boxplots(c,s,"WDF_common_contacts_cast_promoter_slope2")

#%%

# First S129 promoters

s = "bin1_S5p > 0 and bin2_S5p > 0 and (bin1_ATAC_S129 > 0 or bin1_ATAC_common > 0) and (bin2_ATAC_S129 > 0 or bin2_ATAC_common > 0) and (((bin1_H3K27Ac_S129 > 0 or bin1_H3K27Ac_common > 0) and bin2_S129 > 0) or ((bin2_H3K27Ac_S129 > 0 or bin2_H3K27Ac_common > 0) and bin1_S129 > 0)) and dist <= {}".format(distance)

x, all_ = plot_universe(s,"S129 E-P with common 2Mb","s129_specific","S129_EP",2,2,0,False)
list_names = print_genes(x,"s129_specific")
plot_volcano(ase_df[ase_df["gene_name_type"].isin(list_names)],"S129_specific_EP")
x = contacts_with_features_df.query(s + "and (`zscore.dist.cast` < 2 and `zscore.dist.s129` >= 2)") 
print ("Number of contacts like this:")
print (len(x))

# now get different subsets of contacts in files
write_to_pairs(x,"S129_strong_CAST_weak_S129_EP")
write_to_file(x,"S129_strong_CAST_weak_S129_EP")
x = contacts_with_features_df.query(s + "and (`zscore.dist.s129` >= 2)") 
write_to_pairs(x,"S129_strong_S129_EP")
write_to_file(x,"S129_strong_S129_EP")
x = contacts_with_features_df.query(s + "and (`zscore.dist.s129` < 2 and `zscore.dist.cast` >= 2)") 
write_to_pairs(x,"S129_weak_CAST_strong_S129_EP")
write_to_file(x,"S129_weak_CAST_strong_S129_EP")

cast_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` < 2)") 
s129_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` < 2 and `zscore.dist.s129` >= 2)") 
common_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` >= 2)") 
#ll = get_genes_in_contacts(cast_c_total,"s129","gene_names_s129_cast_contacts.txt")
#ll = get_genes_in_contacts(s129_c_total,"s129","gene_names_s129_s129_contacts.txt")
#ll = get_genes_in_contacts(common_c_total,"s129","gene_names_s129_common_contacts.txt")
print ("Number of genes involved in cast_contacts: {}".format(len(print_genes(cast_c_total,"s129_specific"))))
print ("Number of genes involved in s129_contacts: {}".format(len(print_genes(s129_c_total,"s129_specific"))))
print ("Number of genes involved in common_contacts: {}".format(len(print_genes(common_c_total,"s129_specific"))))
do_horizontal_bar_plot("S129_gene_mechanism",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total),0,0],[len(s129_c_total),0,0],[len(common_c_total),0,0]) 

# check if PRC is in the promoters
cast_c_prc_promoter = cast_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter > 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter > 0))") 
s129_c_prc_promoter = s129_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter > 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter > 0))") 
common_c_prc_promoter = common_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter > 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter > 0))") 
cast_c_no_prc_promoter = cast_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter == 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter == 0))") 
s129_c_no_prc_promoter = s129_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter == 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter == 0))") 
common_c_no_prc_promoter = common_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter == 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter == 0))") 

print ("Promoter state of genes involved in cast_contacts: {}".format(get_gene_promoter_states(print_genes(cast_c_total,"s129_specific"))))
print ("Promoter state of genes involved in s129_contacts: {}".format(get_gene_promoter_states(print_genes(s129_c_total,"s129_specific"))))
print ("Promoter state of genes involved in common_contacts: {}".format(get_gene_promoter_states(print_genes(common_c_total,"s129_specific"))))


do_horizontal_bar_prc_plot("S129_gene_mechanism",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total)-len(cast_c_prc_promoter),len(cast_c_prc_promoter),0],[len(s129_c_total)-len(s129_c_prc_promoter),len(s129_c_prc_promoter),0],[len(common_c_total)-len(common_c_prc_promoter),len(common_c_prc_promoter),0]) #each bar total and then each list means: no prc, prc both, prc in one of the anchors
do_horizontal_bar_prc_plot_normalized("S129_gene_mechanism",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total)-len(cast_c_prc_promoter),len(cast_c_prc_promoter),0],[len(s129_c_total)-len(s129_c_prc_promoter),len(s129_c_prc_promoter),0],[len(common_c_total)-len(common_c_prc_promoter),len(common_c_prc_promoter),0]) #each bar total and then each list means: no prc, prc both, prc in one of the anchors

#get WDF of all contacts
c,s = get_WDF_of_contacts(cast_c_total,only_promoters=True,cast_promoters=False,add_bins=2)
wdf_boxplots(c,s,"WDF_cast_contacts_s129_promoter_slope2")
c,s = get_WDF_of_contacts(s129_c_total,only_promoters=True,cast_promoters=False,add_bins=2)
wdf_boxplots(c,s,"WDF_s129_contacts_s129_promoter_slope2")
c,s = get_WDF_of_contacts(common_c_total,only_promoters=True,cast_promoters=False,add_bins=2)
wdf_boxplots(c,s,"WDF_common_contacts_s129_promoter_slope2")

##get WDF of all contacts with PRC in promoter
#c,s = get_WDF_of_contacts(cast_c_prc_promoter,only_promoters=True,cast_promoters=False,add_bins=2)
#wdf_boxplots(c,s,"WDF_cast_contacts_s129_promoter_with_PRC_slope2")
#c,s = get_WDF_of_contacts(s129_c_prc_promoter,only_promoters=True,cast_promoters=False,add_bins=2)
#wdf_boxplots(c,s,"WDF_s129_contacts_s129_promoter_with_PRC_slope2")
#c,s = get_WDF_of_contacts(common_c_prc_promoter,only_promoters=True,cast_promoters=False,add_bins=2)
#wdf_boxplots(c,s,"WDF_common_contacts_s129_promoter_with_PRC_slope2")
#
##get WDF of all contacts with no PRC in promoter
#c,s = get_WDF_of_contacts(cast_c_no_prc_promoter,only_promoters=True,cast_promoters=False,add_bins=2)
#wdf_boxplots(c,s,"WDF_cast_contacts_s129_promoter_with_no_PRC_slope2")
#c,s = get_WDF_of_contacts(s129_c_no_prc_promoter,only_promoters=True,cast_promoters=False,add_bins=2)
#wdf_boxplots(c,s,"WDF_s129_contacts_s129_promoter_with_no_PRC_slope2")
#c,s = get_WDF_of_contacts(common_c_no_prc_promoter,only_promoters=True,cast_promoters=False,add_bins=2)
#wdf_boxplots(c,s,"WDF_common_contacts_s129_promoter_with_no_PRC_slope2")


#%%
#### CTCF LOOPS

# The motif I used from Homer is cccccgggggg mstly. that is + strand in our data but, but in the loop is this way "<". So, to get > <  we need bin1 to be negative and bin2 positive. Bin1 is always upstream of bin2
# https://duckduckgo.com/?q=CTCF+motif+orientation&iar=images&iax=images&ia=images&iai=https%3A%2F%2Fopeni.nlm.nih.gov%2Fimgs%2F512%2F254%2F3251562%2FPMC3251562_pone.0028272.g004.png
# https://ars.els-cdn.com/content/image/1-s2.0-S0092867415009150-fx1_lrg.jpg

# and in our data we can see motif callings like: 
#139(CTGGTGCCCTCTTCTGGTGT,+,0.00)
#	375(TGGCCACTAGGTGGCAGCAA,-,0.00

# First with CAST promoters

s = '(bin1_CTCF_CAST > 0 or bin1_CTCF_common) and (bin1_Rad21_CAST > 0 or bin1_Rad21_common) and (bin2_CTCF_CAST > 0 or bin2_CTCF_common) and (bin2_Rad21_CAST > 0 or bin2_Rad21_common) and (bin1_CAST > 0 or bin2_CAST > 0) and dist <= {} and (bin1_negative_CTCF > 0 and bin2_positive_CTCF > 0)'.format(distance)
x,all_ = plot_universe(s,"CAST/common CTCF loop + cohesin with CAST gene in anchor 2Mb","cast_specific","cast_ctcf_loop",2,2,2,False)
list_names = print_genes(x,"cast_specific")
plot_volcano(ase_df[ase_df["gene_name_type"].isin(list_names)],"CAST_specific_loop")
x = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2)") 
print ("Number of contacts like this:")
print (len(x))

# get lists
write_to_pairs(x,"CAST_strong_CAST_CTCF_loop")
write_to_file(x,"CAST_strong_CAST_CTCF_loop")
x = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` < 2)") 
write_to_pairs(x,"CAST_strong_S129_weak_CAST_CTCF_loop")
write_to_file(x,"CAST_strong_S129_weak_CAST_CTCF_loop")
x = contacts_with_features_df.query(s + "and (`zscore.dist.cast` < 2 and `zscore.dist.s129` >= 2)") 
write_to_pairs(x,"CAST_weak_S129_strong_CAST_CTCF_loop")
write_to_file(x,"CAST_weak_S129_strong_CAST_CTCF_loop")

cast_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` < 2)") # Common and strong
s129_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` < 2 and `zscore.dist.s129` >= 2)") # Common and strong
common_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` >= 2)") # Common and strong
#ll = get_genes_in_contacts(cast_c_total,"cast","ctcf_loops_gene_names_cast_cast_contacts.txt")
#ll = get_genes_in_contacts(s129_c_total,"cast","ctcf_loops_gene_names_cast_s129_contacts.txt")
#ll = get_genes_in_contacts(common_c_total,"cast","ctcf_loops_gene_names_cast_common_contacts.txt")
print ("Number of genes involved in cast_contacts: {}".format(len(print_genes(cast_c_total,"s129_specific"))))
print ("Number of genes involved in s129_contacts: {}".format(len(print_genes(s129_c_total,"s129_specific"))))
print ("Number of genes involved in common_contacts: {}".format(len(print_genes(common_c_total,"s129_specific"))))
do_horizontal_bar_plot("cast_CTCF_loop",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total),0,0],[len(s129_c_total),0,0],[len(common_c_total),0,0]) 


# Check also PRC in promoters
cast_c_prc_promoter = cast_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter > 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter > 0))") # Common and strong
s129_c_prc_promoter = s129_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter > 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter > 0))") # Common and strong
common_c_prc_promoter = common_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter > 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter > 0))") # Common and strong

cast_c_no_prc_promoter = cast_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter == 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter == 0))") # Common and strong
s129_c_no_prc_promoter = s129_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter == 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter == 0))") # Common and strong
common_c_no_prc_promoter = common_c_total.query(s + "and ((bin1_CAST > 0 and bin1_CAST_PRC_promoter == 0) or (bin2_CAST > 0 and bin2_CAST_PRC_promoter == 0))") # Common and strong

do_horizontal_bar_prc_plot("cast_CTCF_loop",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total)-len(cast_c_prc_promoter),len(cast_c_prc_promoter),0],[len(s129_c_total)-len(s129_c_prc_promoter),len(s129_c_prc_promoter),0],[len(common_c_total)-len(common_c_prc_promoter),len(common_c_prc_promoter),0]) #each bar total and then each list means: no prc, prc both, prc in one of the anchors

##get WDF of all contacts
#c,s = get_WDF_of_contacts(cast_c_total,only_promoters=True,cast_promoters=True,add_bins=0)
#wdf_boxplots(c,s,"WDF_cast_contacts_cast_promoter_CTCF")
#c,s = get_WDF_of_contacts(s129_c_total,only_promoters=True,cast_promoters=True,add_bins=0)
#wdf_boxplots(c,s,"WDF_s129_contacts_cast_promoter_CTCF")
#c,s = get_WDF_of_contacts(common_c_total,only_promoters=True,cast_promoters=True,add_bins=0)
#wdf_boxplots(c,s,"WDF_common_contacts_cast_promoter_CTCF")


#%%
# now S129 promoters

s = '(bin1_CTCF_S129 > 0 or bin1_CTCF_common) and (bin1_Rad21_S129 > 0 or bin1_Rad21_common) and (bin2_CTCF_S129 > 0 or bin2_CTCF_common) and (bin2_Rad21_S129 > 0 or bin2_Rad21_common) and (bin1_S129 > 0 or bin2_S129 > 0) and dist <= {} and (bin1_negative_CTCF > 0 and bin2_positive_CTCF > 0)'.format(distance)
x,all_ = plot_universe(s,"S129/common CTCF loop + cohesin with S129 gene in anchor 2Mb","s129_contacts","s129_ctcf_loop",2,2,2,False)
list_names = print_genes(x,"s129_specific")
plot_volcano(ase_df[ase_df["gene_name_type"].isin(list_names)],"S129_specific_loop")
x = contacts_with_features_df.query(s + "and (`zscore.dist.s129` >= 2)") 
print ("Number of contacts like this:")
print (len(x))

# get lists
write_to_pairs(x,"S129_strong_S129_CTCF_loop")
write_to_file(x,"S129_strong_S129_CTCF_loop")
x = contacts_with_features_df.query(s + "and (`zscore.dist.s129` >= 2 and `zscore.dist.cast` < 2)") 
print ("Number of contacts like this:")
print (len(x))
write_to_pairs(x,"S129_strong_CAST_weak_S129_CTCF_loop")
write_to_file(x,"S129_strong_CAST_weak_S129_CTCF_loop")
x = contacts_with_features_df.query(s + "and (`zscore.dist.s129` < 2 and `zscore.dist.cast` >= 2)") 
write_to_pairs(x,"S129_weak_CAST_strong_S129_CTCF_loop")
write_to_file(x,"S129_weak_CAST_strong_S129_CTCF_loop")


cast_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` < 2)") # Common and strong
s129_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` < 2 and `zscore.dist.s129` >= 2)") # Common and strong
common_c_total = contacts_with_features_df.query(s + "and (`zscore.dist.cast` >= 2 and `zscore.dist.s129` >= 2)") # Common and strong
print ("Number of genes involved in cast_contacts: {}".format(len(print_genes(cast_c_total,"s129_specific"))))
print ("Number of genes involved in s129_contacts: {}".format(len(print_genes(s129_c_total,"s129_specific"))))
print ("Number of genes involved in common_contacts: {}".format(len(print_genes(common_c_total,"s129_specific"))))
#ll = get_genes_in_contacts(cast_c_total,"s129","ctcf_loops_gene_names_s129_cast_contacts.txt")
#ll = get_genes_in_contacts(s129_c_total,"s129","ctcf_loops_gene_names_s129_s129_contacts.txt")
#ll = get_genes_in_contacts(common_c_total,"s129","ctcf_loops_gene_names_s129_common_contacts.txt")
do_horizontal_bar_plot("S129_CTCF_loop",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total),0,0],[len(s129_c_total),0,0],[len(common_c_total),0,0]) 

# check prc in promoters
cast_c_prc_promoter = cast_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter > 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter > 0))") # Common and strong
s129_c_prc_promoter = s129_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter > 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter > 0))") # Common and strong
common_c_prc_promoter = common_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter > 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter > 0))") # Common and strong

cast_c_no_prc_promoter = cast_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter == 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter == 0))") # Common and strong
s129_c_no_prc_promoter = s129_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter == 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter == 0))") # Common and strong
common_c_no_prc_promoter = common_c_total.query(s + "and ((bin1_S129 > 0 and bin1_S129_PRC_promoter == 0) or (bin2_S129 > 0 and bin2_S129_PRC_promoter == 0))") # Common and strong

do_horizontal_bar_prc_plot("S129_CTCF_loop",len(cast_c_total),len(s129_c_total),len(common_c_total),[len(cast_c_total)-len(cast_c_prc_promoter),len(cast_c_prc_promoter),0],[len(s129_c_total)-len(s129_c_prc_promoter),len(s129_c_prc_promoter),0],[len(common_c_total)-len(common_c_prc_promoter),len(common_c_prc_promoter),0]) #each bar total and then each list means: no prc, prc both, prc in one of the anchors
#%%

