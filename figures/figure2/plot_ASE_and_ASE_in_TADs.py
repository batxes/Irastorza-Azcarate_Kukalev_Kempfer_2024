#%%

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import scipy.stats 
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

ase_path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"

master_df = pd.read_csv(ase_path,sep="\t")
master_df = master_df[master_df.chrom.isin(chrs_)]

aux = master_df.query('TPM_transcript >= 1 and (log2foldchange >= 1 or log2foldchange <= -1) and p_adj_sum <= 0.05')["gene_name"].drop_duplicates()
aux.to_csv(save_path+"ase_genes.txt",header=False, index=False)

aux = master_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')["gene_name"].drop_duplicates()
aux.to_csv(save_path+"cast_genes.txt",header=False, index=False)

aux = master_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')["gene_name"].drop_duplicates()
aux.to_csv(save_path+"s129_genes.txt",header=False, index=False)

aux = master_df.query('TPM_transcript >= 1')["gene_name"].drop_duplicates()
aux.to_csv(save_path+"all_expressed_genes.txt",header=False, index=False)

aux = master_df.query('TPM_transcript >= 1 and number_SNPs > 0')["gene_name"].drop_duplicates()
aux.to_csv(save_path+"all_expressed_with_SNP_genes.txt",header=False, index=False)

imp_file = "{}imprinted_genes_mouse_updated.csv".format(root)
imp_df = pd.read_csv(imp_file,sep="\t")
imp_genes = imp_df.query('Status == "Imprinted"')
imp_genes_cast = imp_genes.query('Expressed_Allele == "Paternal"')["Gene"].values.tolist()
for x in imp_genes.query('Expressed_Allele == "Paternal"')["Aliases"].values:
    y = x.split(",")
    for gene in y:
        imp_genes_cast.append(gene.strip())

imp_genes_s129 = imp_genes.query('Expressed_Allele == "Maternal"')["Gene"].values.tolist()
for x in imp_genes.query('Expressed_Allele == "Maternal"')["Aliases"].values:
    y = x.split(",")
    for gene in y:
        imp_genes_s129.append(gene.strip())
imp_genes_rest = imp_genes.query('Expressed_Allele != "Maternal" and Expressed_Allele != "Paternal"')["Gene"].values.tolist()
for x in imp_genes.query('Expressed_Allele != "Maternal" and Expressed_Allele != "Paternal"')["Aliases"].values:
    y = x.split(",")
    for gene in y:
        imp_genes_rest.append(gene.strip())

print ('Number of imprinted genes: {}'.format(len(imp_genes)))
imp_s129_df = master_df[master_df.gene_name.isin(imp_genes_s129)]
n_imp_s129 = len(imp_s129_df.query('by_DNAmet == 1'))
imp_cast_df = master_df[master_df.gene_name.isin(imp_genes_cast)]
n_imp_cast = len(imp_cast_df.query('by_DNAmet == 1'))
imp_rest_df = master_df[master_df.gene_name.isin(imp_genes_rest)]
n_imp_rest = len(imp_rest_df.query('by_DNAmet == 1'))
print ('From those, differentially methylated promoters: {}'.format(n_imp_cast+n_imp_s129+n_imp_rest))
n_imp_s129 = len(imp_s129_df.query('type == "cast" or type == "s129"'))
n_imp_cast = len(imp_cast_df.query('type == "cast" or type == "s129"'))
n_imp_rest = len(imp_rest_df.query('type == "cast" or type == "s129"'))
print ('From those, ASE genes: {}'.format(n_imp_cast+n_imp_s129+n_imp_rest))
n_imp_s129 = len(imp_s129_df.query('(type == "cast" or type == "s129") and by_DNAmet == 1'))
n_imp_cast = len(imp_cast_df.query('(type == "cast" or type == "s129") and by_DNAmet == 1'))
n_imp_rest = len(imp_rest_df.query('(type == "cast" or type == "s129") and by_DNAmet == 1'))
print ('From those, ASE genes and methylated promoter: {}'.format(n_imp_cast+n_imp_s129+n_imp_rest))
print(imp_s129_df.query('(type == "cast" or type == "s129") and by_DNAmet == 1'))
print(imp_cast_df.query('(type == "cast" or type == "s129") and by_DNAmet == 1'))
print(imp_rest_df.query('(type == "cast" or type == "s129") and by_DNAmet == 1'))
print ('There are {} ASE genes with differentially methylated promoter'.format(len(master_df.query('(type == "cast" or type == "s129") and by_DNAmet == 1'))))

#pie plot with number of diff genes
fig,ax = plt.subplots(1,1,figsize=(5,5))

fracs = [
    len(master_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')),
    len(master_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')),
    len(master_df.query('TPM_transcript >= 1 and ((log2foldchange < 1 and log2foldchange > -1) or p_adj_sum > 0.05)')),
    len(master_df.query('TPM_transcript >= 1 and number_SNPs == 0')),
    len(master_df.query('TPM_transcript < 1'))
]
radius = 10*0.01+0.6
#patches, texts, autotexts = ax.pie(fracs,
ax.pie(fracs,
                                        autopct='%.0f%%',
                                        textprops={'size': 'smaller'},
                                        shadow=False, radius=radius,
                                        colors=[cast_color,s129_color,"darkgreen","lightgreen","grey"])
ax.set_title("")
#plt.setp(autotexts, size='small')
print ("Number of genes for figure 2: {}".format (fracs))
plt.savefig(save_path+"genes_pieplot.pdf",dpi=300)

print (set(master_df.Promoter_state))


def stacked_bar_plot_prc(act,prca,prcr,k27,inac,label,ax): #I swapped positions of exp genes for plotting
    ax.bar([label],act,label="Active",color="darkgreen")
    ax.bar([label],prca,label="PRCa",bottom=act,color="pink")
    summup = act+ prca
    ax.bar([label],prcr,label="PRCr",bottom=summup,color="red")
    summup = summup+prcr
    ax.bar([label],k27,label="K27me_only",bottom=summup,color="darkred")
    summup = summup+k27
    ax.bar([label],1-summup,label="inactive",bottom=summup,color="grey")
    summup = summup+inac
    ax.bar([label],1-summup,label="Rest",bottom=summup,color="black")

# get specific plots

#not ASE
examples = master_df.query('((log2foldchange < 1 and log2foldchange > -1) or p_adj_sum > 0.05) and TPM_transcript >= 1 and (Promoter_state == "Active" or Promoter_state == "Inactive" or Promoter_state == "PRCa" or Promoter_state == "PRCr" or Promoter_state == "K27me3_only")' )

# bar plot
fig, ax1 = plt.subplots(1,1,figsize=(1,4))
all = len(examples)
print(len(examples.query("Promoter_state == 'Active'"))/all,
                    len(examples.query("Promoter_state == 'PRCa'"))/all,
                    len(examples.query("Promoter_state == 'PRCr'"))/all,
                    len(examples.query("Promoter_state == 'K27me3_only'"))/all,
                    len(examples.query("Promoter_state == 'Inactive'"))/all
                    )
stacked_bar_plot_prc(len(examples.query("Promoter_state == 'Active'"))/all,
                    len(examples.query("Promoter_state == 'PRCa'"))/all,
                    len(examples.query("Promoter_state == 'PRCr'"))/all,
                    len(examples.query("Promoter_state == 'K27me3_only'"))/all,
                    len(examples.query("Promoter_state == 'Inactive'"))/all,
                    "ASE genes",ax1)
plt.legend()
plt.savefig(save_path+"promoter_States_in_non_ase_{}.pdf".format(date),dpi=300)

# pie plot, also SI Figure 6b 
fig,ax = plt.subplots(1,1,figsize=(5,5))
fracs = [
        len(examples.query("Promoter_state == 'Active'"))/all,
        len(examples.query("Promoter_state == 'PRCa'"))/all,
        len(examples.query("Promoter_state == 'PRCr'"))/all,
        len(examples.query("Promoter_state == 'K27me3_only'"))/all,
        len(examples.query("Promoter_state == 'Inactive'"))/all,
        len(examples.query("Promoter_state == 'Internal_gene' or Promoter_state == 'Overlapping_TSSs' or Promoter_state == 'S5p_only' or Promoter_state == 'S7p_only' or Promoter_state != Promoter_state"))/all
    ]

radius = 10*0.01+0.6
#patches, texts, autotexts = ax.pie(fracs,
ax.pie(fracs,
                                        autopct='%.0f%%',
                                        textprops={'size': 'smaller'},
                                        shadow=False, radius=radius, explode = [0,0.03,0.05,0.07,0.03,0.03],
                                        colors=["darkgreen","pink","red","darkred","grey","black"])
ax.set_title("")
#plt.setp(autotexts, size='small')
plt.savefig(save_path+"promoter_State_pieplot.pdf",dpi=300)

# get number of promoter states
print (len(examples))
print (len(examples.query("Promoter_state == 'PRCa'")))
print (len(examples.query("Promoter_state == 'Active'")))
print (len(examples.query("Promoter_state == 'Inactive'")))
print (len(examples.query("Promoter_state == 'PRCr'")))
print (len(examples.query("Promoter_state == 'K27me3_only'")))
print (len(examples.query("Promoter_state == 'Internal_gene' or Promoter_state == 'Overlapping_TSSs' or Promoter_state == 'S5p_only' or Promoter_state == 'S7p_only' or Promoter_state != Promoter_state")))

# get promoter states pie plot of the genes in the HIST 1 locus
examples = master_df[master_df["gene_name"].str.contains('^Hist1')]
all = len(examples)
fig,ax = plt.subplots(1,1,figsize=(5,5))
rest = all-sum([
        len(examples.query("Promoter_state == 'Active'")),
        len(examples.query("Promoter_state == 'PRCa'")),
        len(examples.query("Promoter_state == 'PRCr'")),
        len(examples.query("Promoter_state == 'K27me3_only'")),
        len(examples.query("Promoter_state == 'Inactive'"))])                    
fracs = [
        len(examples.query("Promoter_state == 'Active'"))/all,
        len(examples.query("Promoter_state == 'PRCa'"))/all,
        len(examples.query("Promoter_state == 'PRCr'"))/all,
        len(examples.query("Promoter_state == 'K27me3_only'"))/all,
        len(examples.query("Promoter_state == 'Inactive'"))/all,     
        rest/all               
    ]

radius = 10*0.01+0.6
ax.pie(fracs,
                                        autopct='%.0f%%',
                                        textprops={'size': 'smaller'},
                                        shadow=False, radius=radius, explode = [0,0.03,0.05,0.07,0.03,0.03],
                                        colors=["darkgreen","pink","red","darkred","grey","white"])
ax.set_title("")
plt.savefig(save_path+"promoter_State_pieplot_Hist1_genes.pdf",dpi=300)

print ("Hist genes")
print (len(examples))
print (len(examples.query("Promoter_state == 'PRCa'")))
print (len(examples.query("Promoter_state == 'Active'")))
print (len(examples.query("Promoter_state == 'Inactive'")))
print (len(examples.query("Promoter_state == 'PRCr'")))
print (len(examples.query("Promoter_state == 'K27me3_only'")))


#%%

# get boxplot and violinplot showing the TPM of the gene groups

cast_TPM = master_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')["TPM_transcript"]
s129_TPM = master_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')["TPM_transcript"]
exp_TPM = master_df.query('TPM_transcript >= 1 and ((log2foldchange < 1 and log2foldchange >= -1) or p_adj_sum > 0.05) and number_SNPs > 0')["TPM_transcript"]
exp_no_SNP_TPM = master_df.query('TPM_transcript >= 1 and number_SNPs == 0')["TPM_transcript"]
not_exp_TPM = master_df.query('TPM_transcript < 1')["TPM_transcript"]

palette = [cast_color,s129_color, "darkgreen","limegreen",'grey']
fig,ax = plt.subplots(figsize=(4,8))
vals = [cast_TPM,s129_TPM,exp_TPM,exp_no_SNP_TPM,not_exp_TPM]
bp1 = ax.boxplot(vals,showfliers=False,labels=["CAST","S129","EXP","EXP_no_SNP","Not"],patch_artist=True)
for patch, color in zip(bp1['boxes'], palette):
    patch.set_facecolor(color)
plt.xticks(rotation=45)
xs = []
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.4, color=c)
plt.setp(bp1["medians"],color="black")
ax.set_ylabel("TPM")
plt.savefig(save_path+"genes_transcription.pdf",dpi=300)

for n,p in enumerate(palette):
    for m,q in enumerate(palette):
        print ("{} vs {}".format(palette[n],palette[m]))
        print (scipy.stats.ttest_ind(vals[n],vals[m],equal_var=False)[1])

fig,ax = plt.subplots(figsize=(4,8))
ax.violinplot(vals)


#%%

# Volcano plot shown in Figure 2b
# I am also plotting a filtered volcano plot showing imprinted genes

plot_all_text = False
#plot_name = "TPM < 1 & PRCr"
#plot_name = "all not inactive genes"
plot_name = "TPM >= 1"

figname = "volcano_TPM1_{}.svg".format(date)
figname2 = "volcano_TPM1_imprinted_{}.svg".format(date)
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
#examples = s129_df.query('(log2foldchange < -6  or p_adj_sum < 1E-180) and log2foldchange != -inf and p_adj_sum != 0 and  transcript_biotype == "protein_coding"' )
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
plt.savefig("{}{}".format(save_path,figname),dpi=300)

#imprinted genes:
fig,ax = plt.subplots(1,1,figsize=(8,10))
imp_rest_df = master_df[master_df.gene_name.isin(imp_genes_rest)]
x = imp_rest_df["log2foldchange"]
y = np.log10(imp_rest_df["p_adj_sum"])
ax.scatter(x,-1*y,color="black",s=size)

examples = imp_rest_df.query('(log2foldchange >= 1 or log2foldchange <= -1) and p_adj_sum <= 0.05')
x = examples.log2foldchange.values.tolist() 
y = -1*np.log10(examples.p_adj_sum.values.tolist())
names = examples.gene_name.values.tolist()
texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
adjust_text(texts)

imp_cast_df = master_df[master_df.gene_name.isin(imp_genes_cast)]
x = imp_cast_df["log2foldchange"]
y = np.log10(imp_cast_df["p_adj_sum"])
ax.scatter(x,-1*y,color="violet",s=size)

examples = imp_cast_df.query('(log2foldchange >= 1 or log2foldchange <= -1) and p_adj_sum <= 0.05')
x = examples.log2foldchange.values.tolist() 
y = -1*np.log10(examples.p_adj_sum.values.tolist())
names = examples.gene_name.values.tolist()
texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
adjust_text(texts)

imp_s129_df = master_df[master_df.gene_name.isin(imp_genes_s129)]
x = imp_s129_df["log2foldchange"]
y = np.log10(imp_s129_df["p_adj_sum"])
ax.scatter(x,-1*y,color="gold",s=size)

examples = imp_s129_df.query('(log2foldchange >= 1 or log2foldchange <= -1) and p_adj_sum <= 0.05')
x = examples.log2foldchange.values.tolist() 
y = -1*np.log10(examples.p_adj_sum.values.tolist())
names = examples.gene_name.values.tolist()
texts = [ax.text(x[i], y[i], names[i], ha='center', va='center') for i in range(len(x))]
adjust_text(texts)


p_adj_cutoff = (-1*np.log10(0.05))
ax.axhline(p_adj_cutoff,color="black",linestyle="--")
ax.axvline(-1,color="black",linestyle="--")
ax.axvline(1,color="black",linestyle="--")
ax.set_ylabel("$-log_{10}$(P-value)")
ax.set_xlabel("$log_2(FoldChange)$")
#ax.set_ylim(-10,100)
#ax.set_xlim(-3,3)
plt.savefig("{}{}".format(save_path,figname2),dpi=300)

#%%

# These are different plots characterizing the gene types.
# They were not showed in the final article.

import pybedtools as pb
import re
from os import listdir
import operator
from functools import reduce
import math


def sort_df(df, column_idx, key):
    '''Takes a dataframe, a column index and a custom function for sorting, 
    returns a dataframe sorted by that column using that function'''

    col = df.loc[:,column_idx]
    df = df.reindex([i[1] for i in sorted(zip(col,range(len(col))), key=key)])
    return df

def draw_boxplot(vals,names,ytext,palette,path,ylim):
    #font = {'family' : 'normal',
    #        'weight' : 'bold',
    #        'size'   : 22}
    #matplotlib.rc('font', **font)
    fig,ax = plt.subplots(figsize=(4,8))
    bp1 = plt.boxplot(vals,labels=names,showfliers=False,notch=True)
    xs = []
    for i, col in enumerate(range(len(vals))):
        xs.append(np.random.normal(i + 1, 0.04, len(vals[i])))
    for x, val, c in zip(xs, vals, palette):
        plt.scatter(x, val, alpha=0.4, color=c)
    ax.set_ylabel(ytext)
    #font = {'family': 'sansserif',
    #        'color':  'black',
    #        'weight': 'normal',
    #        'size': 16,
    #        }
    #plt.text(1.1, 70,"n: {}".format(len(vals[0])) , fontdict=font)
    #plt.text(2.1, 70,"n: {}".format(len(vals[1])) , fontdict=font)
    plt.ylim(0,ylim)
    plt.setp(bp1["medians"],color="black")
    plt.savefig(path)
    #plt.savefig("/home/ibai/Number_of_exp_genes_per_{}_binsize_{}_in_ASE_and_non_ASE_bins.png".format(bin_size,name),dpi=300)

cols = master_df.columns.values.tolist()
#set gene_name column in the 5th place
#cols = [cols[1:5],[cols[0]],cols[5:]]
cols = [[cols[0]],cols[3:5],cols[1:3],cols[5:]]
cols = reduce(operator.concat, cols)
master_df_aux = master_df[cols]
#ase_df_aux.rename(columns={"chr": "chrom"},inplace=True)
master_df_aux.rename(columns={"start_transcript": "start","end_transcript":"end"},inplace=True)
# function for sorting
# take only protein coding if there are both for a gene
#cmp = lambda x:2 if 'protein_coding' in x else 1 if 'lncRNA' in x else 0     
#master_df_aux = sort_df(master_df_aux,'transcript_biotype',cmp).drop_duplicates('gene_name', keep="last")
master_df_aux = master_df_aux.sort_values(by=['chrom', 'start'])
master_df_aux.dropna(subset=["start"],inplace=True)
master_df_aux["start"] = master_df_aux["start"].astype(int)
master_df_aux["end"] = master_df_aux["end"].astype(int)
master_df_aux["size"] = master_df_aux["end"]- master_df_aux["start"]
exp_df = master_df_aux.query('TPM_transcript >= 1')
cast_genes_df = exp_df.query('log2foldchange >= 1 & p_adj_sum <= 0.05')
print ("cast genes: {}".format(len(cast_genes_df)))
s129_genes_df = exp_df.query('log2foldchange <= -1 & p_adj_sum <= 0.05')
print ("s129 genes: {}".format(len(s129_genes_df)))
mono_cast_genes_df = exp_df.query('log2foldchange == "inf" & p_adj_sum <= 0.05')
mono_s129_genes_df = exp_df.query('log2foldchange == "-inf" & p_adj_sum <= 0.05')
exp_SNP_df = exp_df.query('number_SNPs >= 1 and ((log2foldchange < 1 and log2foldchange > -1) or p_adj_sum > 0.05)')
print ("exp SNP genes: {}".format(len(exp_SNP_df)))
#exp_no_SNP_df = exp_df.query('number_SNPs == "NaN"')
exp_no_SNP_df = exp_df.query('number_SNPs == 0')
print ("exp no SNP genes: {}".format(len(exp_no_SNP_df)))
not_exp_df = master_df_aux.query('TPM_transcript < 1')
print ("not exp genes: {}".format(len(not_exp_df)))

cast_pb = pb.BedTool.from_dataframe(cast_genes_df[["chrom","start","end"]]).sort()
s129_pb = pb.BedTool.from_dataframe(s129_genes_df[["chrom","start","end"]]).sort()
exp_pb = pb.BedTool.from_dataframe(exp_SNP_df[["chrom","start","end"]]).sort()
exp_no_SNP_pb = pb.BedTool.from_dataframe(exp_no_SNP_df[["chrom","start","end"]]).sort()
not_exp_pb = pb.BedTool.from_dataframe(not_exp_df[["chrom","start","end"]]).sort()

closest_cast_to_exp_df = cast_pb.closest(exp_pb,d=True,t="first",io=True).to_dataframe()
closest_s129_to_exp_df = s129_pb.closest(exp_pb,d=True,t="first",io=True).to_dataframe()
closest_exp_no_SNP_to_exp_df = exp_no_SNP_pb.closest(exp_pb,d=True,t="first",io=True).to_dataframe()
print (closest_exp_no_SNP_to_exp_df)
closest_not_exp_to_exp_df = not_exp_pb.closest(exp_pb,d=True,t="first",io=True).to_dataframe()

closest_cast_to_cast_df = cast_pb.closest(cast_pb,d=True,t="first",io=True).to_dataframe()
closest_s129_to_s129_df = s129_pb.closest(s129_pb,d=True,t="first",io=True).to_dataframe()
closest_exp_no_SNP_to_exp_no_SNP_df = exp_no_SNP_pb.closest(exp_no_SNP_pb,d=True,t="first",io=True).to_dataframe()
closest_exp_SNP_to_exp_SNP_df = exp_pb.closest(exp_pb,d=True,t="first",io=True).to_dataframe()
closest_not_exp_to_not_exp_df = not_exp_pb.closest(not_exp_pb,d=True,t="first",io=True).to_dataframe()

chrom_size_path = root+"mm10.chrom.sizes"
bin_size = 1000000
a = pb.example_bedtool('a.bed')
aux_pb = pb.BedTool.window_maker(a,g=chrom_size_path,w=bin_size)
density_cast_df = aux_pb.intersect(cast_pb,c=True).to_dataframe()
density_cast_df = density_cast_df[density_cast_df["chrom"].isin(chrs_)]
density_s129_df = aux_pb.intersect(s129_pb,c=True).to_dataframe()
density_s129_df = density_s129_df[density_s129_df["chrom"].isin(chrs_)]
density_exp_SNP_df = aux_pb.intersect(exp_pb,c=True).to_dataframe()
density_exp_SNP_df = density_exp_SNP_df[density_exp_SNP_df["chrom"].isin(chrs_)]
density_exp_no_SNP_df = aux_pb.intersect(exp_no_SNP_pb,c=True).to_dataframe()
density_exp_no_SNP_df = density_exp_no_SNP_df[density_exp_no_SNP_df["chrom"].isin(chrs_)]
density_not_exp_df = aux_pb.intersect(not_exp_pb,c=True).to_dataframe()
density_not_exp_df = density_not_exp_df[density_not_exp_df["chrom"].isin(chrs_)]

#another type of boxplots

font = {'family': 'sansserif',
        'color':  'black',
        'weight': 'normal',
        'size': 8,
        }
labels = ["CAST","S129","Expressed SNP","Expressed no SNP","Not Expressed"]
fig,ax = plt.subplots(figsize=(2,4))
#bp1 = ax.boxplot([cast_genes_df["size"].values.tolist(),s129_genes_df["size"].values.tolist(),
bp1 = ax.boxplot([cast_genes_df["size"].values.tolist(),s129_genes_df["size"].values.tolist(),
exp_SNP_df["size"].values.tolist(),exp_no_SNP_df["size"],not_exp_df["size"].values.tolist()],
showfliers=False,patch_artist=True,labels=labels)
colors = [cast_color,s129_color, "lightgreen","green",'grey']
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp1["medians"],color="black")
plt.xticks(rotation=45)
ax.set_xlabel("Genes")
ax.set_ylabel("Gene length")
#ax.text(2.5, 18,"n: {}".format(len(s129_b_exp)) , fontdict=font)
fig.savefig(save_path+"gene_length_{}.png".format(date),dpi=300,bbox_inches='tight')

fig,ax = plt.subplots(figsize=(2,4))
bp1 = ax.boxplot([density_cast_df["name"].values.tolist(),density_s129_df["name"].values.tolist(),
density_exp_SNP_df["name"].values.tolist(),density_exp_no_SNP_df["name"].values.tolist(),density_not_exp_df["name"].values.tolist()],
showfliers=False,patch_artist=True,labels=labels)
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp1["medians"],color="black")
plt.xticks(rotation=45)
ax.set_xlabel("Genes")
ax.set_ylabel("Gene density per 1Mb")
#ax.text(2.5, 18,"n: {}".format(len(s129_b_exp)) , fontdict=font)
fig.savefig(save_path+"gene_density_{}.png".format(date),dpi=300,bbox_inches='tight')

fig,ax = plt.subplots(figsize=(2,4))
bp1 = ax.boxplot([cast_genes_df["TPM_transcript"].values.tolist(),s129_genes_df["TPM_transcript"].values.tolist(),
exp_SNP_df["TPM_transcript"].values.tolist(),exp_no_SNP_df["TPM_transcript"].values.tolist(),not_exp_df["TPM_transcript"].values.tolist()],
showfliers=False,patch_artist=True,labels=labels)
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp1["medians"],color="black")
plt.xticks(rotation=45)
ax.set_xlabel("Genes")
ax.set_ylabel("TPM")
#ax.text(2.5, 18,"n: {}".format(len(s129_b_exp)) , fontdict=font)
fig.savefig(save_path+"TPM_{}.png".format(date),dpi=300,bbox_inches='tight')

fig,ax = plt.subplots(figsize=(2,4))
bp1 = ax.boxplot([closest_cast_to_exp_df["thickStart"].values.tolist(),closest_s129_to_exp_df["thickStart"].values.tolist(),
closest_exp_SNP_to_exp_SNP_df["thickStart"].values.tolist(),
closest_exp_no_SNP_to_exp_df["thickStart"].values.tolist(),
closest_not_exp_to_exp_df["thickStart"].values.tolist()],
showfliers=False,patch_artist=True,labels=labels)
ax.set_yscale('log')
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp1["medians"],color="black")
plt.xticks(rotation=45)
ax.set_xlabel("Genes")
ax.set_ylabel("Distance to expressed gene in bps (log)")
#ax.text(2.5, 18,"n: {}".format(len(s129_b_exp)) , fontdict=font)
fig.savefig(save_path+"distance_{}.png".format(date),dpi=300,bbox_inches='tight')

fig,ax = plt.subplots(figsize=(2,4))
bp1 = ax.boxplot([closest_cast_to_cast_df["thickStart"].values.tolist(),closest_s129_to_s129_df["thickStart"].values.tolist(),
closest_exp_SNP_to_exp_SNP_df["thickStart"].values.tolist(),closest_exp_no_SNP_to_exp_no_SNP_df["thickStart"].values.tolist(),closest_not_exp_to_not_exp_df["thickStart"].values.tolist()],
showfliers=False,patch_artist=True,labels=labels)
ax.set_yscale('log')
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp1["medians"],color="black")
plt.xticks(rotation=45)
ax.set_xlabel("Genes")
ax.set_ylabel("Distance to same type gene in bps (log)")
#ax.text(2.5, 18,"n: {}".format(len(s129_b_exp)) , fontdict=font)
fig.savefig(save_path+"distance_same_type_{}.png".format(date),dpi=300,bbox_inches='tight')


#%%
####the same with violinplots

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

def plot_violin(data_list,ylabel,savepath,log=False):
    fig,ax = plt.subplots(figsize=(2,4))
    bp1 = ax.violinplot(data_list,
    showmedians=True,showextrema=False)
    for patch, color in zip(bp1['bodies'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(1)
    #for partname in ('crplot
    # bars','cmins','cmaxes','cmedians'):
    for partname in ('','cmedians'):
        try:
            vp = bp1[partname]
            vp.set_edgecolor("black")
            #vp.set_linewidth(1)
        except:
            continue
    set_axis_style(ax, labels)
    plt.xticks(rotation=45)
    ax.set_xlabel("Genes")
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale('log')
    fig.savefig(savepath,dpi=300,bbox_inches='tight')
    
plot_violin([cast_genes_df["size"].values.tolist(),s129_genes_df["size"].values.tolist(),
    exp_SNP_df["size"].values.tolist(),exp_no_SNP_df["size"],not_exp_df["size"].values.tolist()],
    "Gene length",save_path+"Gene_length_violin_{}.pdf".format(date))

plot_violin([density_cast_df["name"].values.tolist(),density_s129_df["name"].values.tolist(),
    density_exp_SNP_df["name"].values.tolist(),density_exp_no_SNP_df["name"].values.tolist(),density_not_exp_df["name"].values.tolist()],
    "Gene_density per 1Mb",save_path+"gene_density_violin_{}.pdf".format(date))
plot_violin([cast_genes_df["TPM_transcript"].values.tolist(),s129_genes_df["TPM_transcript"].values.tolist(),
    exp_SNP_df["TPM_transcript"].values.tolist(),exp_no_SNP_df["TPM_transcript"].values.tolist(),not_exp_df["TPM_transcript"].values.tolist()],
    "TPM",save_path+"TPM_violin.pdf")
plot_violin([closest_cast_to_exp_df["thickStart"].values.tolist(),closest_s129_to_exp_df["thickStart"].values.tolist(),
    closest_exp_SNP_to_exp_SNP_df["thickStart"].values.tolist(),closest_exp_no_SNP_to_exp_df["thickStart"].values.tolist(),closest_not_exp_to_exp_df["thickStart"].values.tolist()],
    "Distance to expressed gene in bps (log)",save_path+"distance_violin_{}.pdf".format(date),True)
plot_violin([closest_cast_to_cast_df["thickStart"].values.tolist(),closest_s129_to_s129_df["thickStart"].values.tolist(),
    closest_exp_SNP_to_exp_SNP_df["thickStart"].values.tolist(),closest_exp_no_SNP_to_exp_df["thickStart"].values.tolist(),closest_not_exp_to_not_exp_df["thickStart"].values.tolist()],
    "Distance to same type gene in bps (log)",save_path+"distance_same_type_{}.pdf".format(date),True)

#%%

# plots shown in Figure2a and Figure2d

# homogenize the gene expression data to analyze specific gene groups
master_df = pd.read_csv(ase_path,sep="\t")
master_df = master_df[master_df.chrom.isin(chrs_)]
master_df = master_df.drop_duplicates()
cols = master_df.columns.values.tolist()
cols = [[cols[0]],cols[3:5],cols[1:3],cols[5:]]
cols = reduce(operator.concat, cols)
master_df = master_df[cols]
master_df.rename(columns={"start_transcript": "start","end_transcript":"end"},inplace=True)
cmp = lambda x:2 if 'protein_coding' in x else 1 if 'lncRNA' in x else 0  
master_df = master_df.sort_values(by=['chrom', 'start'])
master_df.dropna(subset=["start"],inplace=True)
master_df["start"] = master_df["start"].astype(int)
master_df["end"] = master_df["end"].astype(int)
master_df["gene_size"] = master_df["end"] - master_df["start"]

rest_df = master_df.query('TPM_transcript >= 1 and number_SNPs > 0 and ((log2foldchange < 1 and log2foldchange > -1) or p_adj_sum > 0.05)')
rest_no_SNP_df = master_df.query('TPM_transcript >= 1 and number_SNPs == 0')
cast_df = master_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')
s129_df = master_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')
lnc_df = master_df.query('transcript_biotype == "lncRNA"')

#get metabolic genes
#metaboli genes from GO http://www.informatics.jax.org/go/term/GO:0008152
meta_list = root+"GO_term_summary_20221117_070849.csv"
meta_df = pd.read_csv(meta_list,sep="\t")
meta_list = meta_df.Symbol.values.tolist()
meta_list = list(set(meta_list))

#HK genes housekeeping genes plot
hk_list = root+"Housekeeping_TranscriptsMouse.csv"
hk_df = pd.read_csv(hk_list,sep=";")
hk_list = hk_df.Genes.values.tolist()

#ŧhe list was changed for an updated one (taking into account genes in the ASE gene list):
ribo_list = root+"Mouse_ribosomal_proteins.txt"
ribo_list = root+"Mouse_ribosomal_proteins_updated.txt"
ribo_df = pd.read_csv(ribo_list,sep="\t",names=["subunit","gene_name"])
ribo_list = ribo_df.gene_name.values.tolist()

# cell cycle gene list
cc_list = root+"Mouse_cell_cycle.txt"
cc_df = pd.read_csv(cc_list,sep="\t",names=["gene_name"])
cc_list = cc_df.gene_name.values.tolist()

#histone genes
hist_list = root+"mouse_histone_genes.txt"
hist_df = pd.read_csv(hist_list,sep="\t",names=["gene_name"])
hist_list = hist_df.gene_name.values.tolist()

# map kinases
mapk_list = root+"mouse_mapk_signalling_genes.txt"
mapk_df = pd.read_csv(mapk_list,sep="\t",names=["gene_name"])
mapk_list = mapk_df.gene_name.values.tolist()

lnc_rest_df = lnc_df.query('TPM_transcript >= 1 and number_SNPs > 0 and ((log2foldchange < 1 and log2foldchange > -1) or p_adj_sum > 0.05)')
lnc_rest_no_SNP_df = lnc_df.query('TPM_transcript >= 1 and number_SNPs == 0')
lnc_cast_df = lnc_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')
lnc_s129_df = lnc_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')
n_lnc = len(lnc_df)
n_lnc_rest = len(lnc_rest_df)
n_lnc_no_SNP_rest = len(lnc_rest_no_SNP_df)
n_lnc_cast = len(lnc_cast_df)
n_lnc_s129 = len(lnc_s129_df)
lnc_list = lnc_df.gene_name.values.tolist()

def intersection(lst1, lst2):
    lst2 = [lst.lower() for lst in lst2]
    lst3 = [value for value in lst1 if value.lower() in lst2]
    return lst3

rest_list = rest_df.gene_name.values.tolist() 
rest_no_SNP_list = rest_no_SNP_df.gene_name.values.tolist() 

cast_list = cast_df.gene_name.values.tolist() 
s129_list = s129_df.gene_name.values.tolist() 
n_rest = len(rest_list)
n_rest_no_SNP = len(rest_no_SNP_list)
n_cast = len(cast_list)
n_s129 = len(s129_list)
n_exp = n_rest+n_cast+n_s129+n_rest_no_SNP
n_genes = len(master_df)
n_hist = len(hist_list)
n_cc = len(cc_list)
n_ribo = len(ribo_list)
n_hk = len(hk_list)
n_mapk = len(mapk_list)
n_meta = len(meta_list)

rest_no_SNP_list = [x for x in rest_no_SNP_list if x not in rest_list]

n_hk_rest = len(intersection(hk_list, rest_list))
n_hk_rest_no_SNP = len(intersection(hk_list, rest_no_SNP_list))
n_hk_cast = len(intersection(hk_list, cast_list))
n_hk_s129 = len(intersection(hk_list,s129_list))
print (n_hk_rest,n_hk_rest_no_SNP,n_hk_cast,n_hk_s129,n_hk)
print (len(intersection(rest_list,rest_no_SNP_list)))

n_ribo_rest = (len(intersection(ribo_list, rest_list)))
n_ribo_rest_no_SNP = len(intersection(ribo_list, rest_no_SNP_list))
n_ribo_cast = (len(intersection(ribo_list, cast_list)))
n_ribo_s129 = (len(intersection(ribo_list, s129_list)))

n_cc_rest = (len(intersection(cc_list, rest_list)))
n_cc_rest_no_SNP = len(intersection(cc_list, rest_no_SNP_list))
n_cc_cast = (len(intersection(cc_list, cast_list)))
n_cc_s129 = (len(intersection(cc_list, s129_list)))

n_hist_rest = (len(intersection(hist_list, rest_list)))
n_hist_rest_no_SNP = len(intersection(hist_list, rest_no_SNP_list))
n_hist_cast = (len(intersection(hist_list, cast_list)))
n_hist_s129 = (len(intersection(hist_list, s129_list)))

n_mapk_rest = (len(intersection(mapk_list, rest_list)))
n_mapk_rest_no_SNP = len(intersection(mapk_list, rest_no_SNP_list))
n_mapk_cast = (len(intersection(mapk_list, cast_list)))
n_mapk_s129 = (len(intersection(mapk_list, s129_list)))

n_meta_rest = (len(intersection(meta_list, rest_list)))
n_meta_rest_no_SNP = len(intersection(meta_list, rest_no_SNP_list))
n_meta_cast = (len(intersection(meta_list, cast_list)))
n_meta_s129 = (len(intersection(meta_list, s129_list)))


def stacked_bar_plot(cast_per,s129_per,exp_no_SNP_per,exp_per,label,ax): #I swapped positions of exp genes for plotting
    print (label)
    print (cast_per,s129_per,exp_no_SNP_per,exp_per)
    ax.bar([label],cast_per,label="CAST",color=cast_color)
    ax.bar([label],s129_per,label="S129",bottom=cast_per,color=s129_color)
    summup =cast_per+ s129_per
    ax.bar([label],exp_no_SNP_per,label="exp_wo_SNP",bottom=summup,color="lightgreen")
    summup = summup+exp_no_SNP_per
    ax.bar([label],exp_per,label="exp_with_SNP",bottom=summup,color="darkgreen")
    summup = summup+exp_per
    ax.bar([label],1-summup,label="not_exp",bottom=summup,color="grey")


fig,axs = plt.subplots(1,8,figsize=(5,6),sharey=True)
stacked_bar_plot(n_cast/n_genes,
                    n_s129/n_genes,
                    n_rest_no_SNP/n_genes,
                    n_rest/n_genes,
                    "All genes",axs[0])
stacked_bar_plot(n_hist_cast/n_hist,
                    n_hist_s129/n_hist,
                    n_hist_rest_no_SNP/n_hist,
                    n_hist_rest/n_hist,
                    "Histone",axs[4])
stacked_bar_plot(n_cc_cast/n_cc,
                    n_cc_s129/n_cc,
                    n_cc_rest_no_SNP/n_cc,
                    n_cc_rest/n_cc,
                    "Cell cycle",axs[2])
stacked_bar_plot(n_ribo_cast/n_ribo,
                    n_ribo_s129/n_ribo,
                    n_ribo_rest_no_SNP/n_ribo,
                    n_ribo_rest/n_ribo,
                    "Ribosomal",axs[3])
stacked_bar_plot(n_meta_cast/n_meta,
                    n_meta_s129/n_meta,
                    n_meta_rest_no_SNP/n_meta,
                    n_meta_rest/n_meta,
                    "Metabolic",axs[7])
stacked_bar_plot(n_hk_cast/n_hk,
                    n_hk_s129/n_hk,
                    n_hk_rest_no_SNP/n_hk,
                    n_hk_rest/n_hk,
                    "HK",axs[1])
stacked_bar_plot(n_mapk_cast/n_mapk,
                    n_mapk_s129/n_mapk,
                    n_mapk_rest_no_SNP/n_mapk,
                    n_mapk_rest/n_mapk,
                    "MAPK",axs[5])
stacked_bar_plot(n_lnc_cast/n_lnc,
                    n_lnc_s129/n_lnc,
                    n_lnc_no_SNP_rest/n_lnc,
                    n_lnc_rest/n_lnc,
                    "lncRNA",axs[6])
axs[0].set_ylabel("Percentage")
axs[0].set_title(n_genes)
axs[1].set_title(n_hk)
axs[2].set_title(n_cc)
axs[3].set_title(n_ribo)
axs[4].set_title(n_hist)
axs[5].set_title(n_mapk)
axs[6].set_title(n_lnc)
axs[7].set_title(n_meta)
axs[0].set_yticklabels([0,20,40,60,80,100])
plt.legend(loc="upper right")
plt.savefig(save_path+"types_of_genes_percentage_{}.pdf".format(date),dpi=300,bbox_inches='tight')

meta_df = meta_df.rename(columns={"Symbol":"gene_name"})
meta_df = meta_df.drop_duplicates(subset=["gene_name"])
master_and_meta_df = pd.merge(master_df,meta_df,on=["gene_name"]) 
master_and_meta_df.to_csv("{}ase_and_meta_df.tsv".format(save_path),sep="\t")
print (len(master_df))
print (len(master_df.query("Promoter_state == 'PRCa' or Promoter_state == 'PRCr' or Promoter_state == 'K27me3_only'")))

print(len (master_and_meta_df))
print (len(master_and_meta_df.query("Promoter_state == 'PRCa' or Promoter_state == 'PRCr' or Promoter_state == 'K27me3_only'")))

allelic_df = master_df.query('TPM_transcript >= 1 and (log2foldchange >= 1 or log2foldchange <= -1) and p_adj_sum <= 0.05')
print (len(allelic_df))
print (len(allelic_df.query("Promoter_state == 'PRCa' or Promoter_state == 'PRCr' or Promoter_state == 'K27me3_only'")))

allelic_meta_df = master_and_meta_df.query('TPM_transcript >= 1 and (log2foldchange >= 1 or log2foldchange <= -1) and p_adj_sum <= 0.05')
print (len(allelic_meta_df))
print (len(allelic_meta_df.query("Promoter_state == 'PRCa' or Promoter_state == 'PRCr' or Promoter_state == 'K27me3_only'")))


# %%

# similar plot but showing number of HK, ribo histone genes and other genes in the differen groups

def stacked_bar_plot_genes(cast_per,s129_per,exp_per,label,ax): #names are also changed 
    ax.bar([label],cast_per,label="HK",color="steelblue")
    ax.bar([label],s129_per,label="Ribo",bottom=cast_per,color="salmon")
    summup =cast_per+ s129_per
    ax.bar([label],exp_per,label="Hist",bottom=summup,color="gold")
    summup = summup+exp_per
    ax.bar([label],1-summup,label="Other",bottom=summup,color="darkgrey")

fig,axs = plt.subplots(1,4,figsize=(3,6),sharey=True)

stacked_bar_plot_genes(n_hk_rest/n_rest,
                    n_ribo_rest/n_rest,
                    n_hist_rest/n_rest,
                    "EXP",axs[0])
stacked_bar_plot_genes(n_hk_rest_no_SNP/n_rest_no_SNP,
                    n_ribo_rest_no_SNP/n_rest_no_SNP,
                    n_hist_rest_no_SNP/n_rest_no_SNP,
                    "EXP WO SNP",axs[1])
stacked_bar_plot_genes(n_hk_cast/n_cast,
                    n_ribo_cast/n_cast,
                    n_hist_cast/n_cast,
                    "CAST",axs[2])
stacked_bar_plot_genes(n_hk_s129/n_s129,
                    n_ribo_s129/n_s129,
                    n_hist_s129/n_s129,
                    "S129",axs[3])

axs[0].set_ylabel("Percentage")
axs[0].set_title(n_rest)
axs[1].set_title(n_rest_no_SNP)
axs[2].set_title(n_cast)
axs[3].set_title(n_s129)
plt.legend()
plt.savefig(save_path+"types_of_genes_in_gene_groups_{}.pdf".format(date),dpi=300,bbox_inches='tight')

# %%
# Violin plot showing Number of SNPs in different groups of genes
# pretty interesting that ribosomal genes have so many SNPs.

vector = []
info = "gene_size"
info = "number_SNPs"
for list_ in [master_df.gene_name.values.tolist(),hk_list,cc_list,ribo_list,hist_list,mapk_list]:
    list_ = [value.lower() for value in list_]
    aux_df = master_df.loc[master_df.gene_name.str.lower().isin(list_)]
    #aux_df = aux_df.loc[aux_df.transcript_biotype == "protein_coding"]
    print (np.median(aux_df[info].values.tolist()))
    aux_df["n_SNPs_normalized"] = aux_df.number_SNPs/aux_df.gene_size
    #aux_df["n_SNPs_normalized"] = aux_df.number_exon_SNPs/aux_df.gene_size
    #aux_df["n_SNPs_normalized"] = aux_df.number_intron_SNPs/aux_df.gene_size
    v =aux_df["n_SNPs_normalized"].values.tolist()
    vector.append(v)

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

fig,ax = plt.subplots(1,1,figsize=(12,12))
parts = ax.violinplot(vector,showmeans=False, showmedians=True,showextrema=True)
for pc in parts['bodies']:
    pc.set_alpha(1)
    pc.set_edgecolor([0,0,0,1])
    pc.set_facecolor('#FFFFFF')
    pc.set_facecolor([1,1,1,0])
for pos,v in enumerate(vector): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax.scatter(xv, v,s=1,alpha=1,color="black")
labels = ["all genes","HK","Cell cycle","Ribo","Hist","MAPK"]
set_axis_style(ax, labels)
ax.set_ylabel("number of SNPs in genes normalized by size")
fig.savefig(save_path+"Number_of_SNPs_in_gene_groups.pdf",dpi=300,bbox_inches='tight')


# %%
#check ratios with different gene sizes

fig,axs = plt.subplots(1,15,figsize=(10,6),sharey=True)
stacked_bar_plot(n_cast/n_genes,
                    n_s129/n_genes,
                    n_rest_no_SNP/n_genes,
                    n_rest/n_genes,
                    "All genes",axs[0])
axs[0].set_title(n_genes)
sizes = [50000,25000,10000,5000,2000,1500,1000,975,950,925,900,750,500,250]
sizes_str = ["50000","25000","10000","5000","2000","1500","1000","975","950","925","900","750","500","250"]
for n,s in enumerate(sizes):
    m_df = master_df.loc[master_df.gene_size < s]
    print (len(m_df))
    rest_df = m_df.query('TPM_transcript >= 1 and number_SNPs > 0 and (log2foldchange < 1 and log2foldchange > -1)')
    rest_no_SNP_df = m_df.query('TPM_transcript >= 1 and number_SNPs == 0')
    cast_df = m_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')
    s129_df = m_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')
    n_cast_ = len(cast_df) 
    n_s129_ = len(s129_df) 
    n_rest_ = len(rest_df)
    n_rest_no_SNP_ = len(rest_no_SNP_df)
    n_genes_ = len(m_df)

    stacked_bar_plot(n_cast_/n_genes_,
                        n_s129_/n_genes_,
                        n_rest_no_SNP_/n_genes_,
                        n_rest_/n_genes_,
                        sizes_str[n],axs[n+1])
    axs[n+1].set_title(n_genes_)
fig.savefig(save_path+"proportion_of_genes_depending_on_their_size.pdf",dpi=300,bbox_inches='tight')
# %%

# Now I started adding TAD information to do the analysis and comapre with ASE genes.

#allelic function canonical
def allelic_func(row):
    if row['TPM_transcript'] >= 1:
        if row['number_SNPs'] > 0:
            if row["p_adj_sum"] <= 0.05:
                if row["log2foldchange"] <= -1:
                    val = "S129"
                elif row["log2foldchange"] >= 1:
                    val = "CAST"
                else:
                    val = 'exp'
            else:
                val = 'exp'
        else:
            val = 'exp_no_SNP'
    else:
        val = 'no'
    return val

# calling allelic more relaxed
def allelic_relaxed_func12312_wronggg(row):
    if row['TPM_transcript'] >= 1:
        if row["p_adj_sum"] <= 0.05:
            if row["log2foldchange"] < -0:
                val = "S129"
            elif row["log2foldchange"] > 0:
                val = "CAST"
            else:
                val = 'exp'
        else:
            val = 'exp'
    else:
        val = 'no'
    return val

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

def find_cast_TAD_func_gene(row):
    start = row['start']
    tads = tads_cast_df.query('chrom == "{}"'.format(row['chrom']))
    previous_start = 3000000
    previous_end = 10000000
    for start_tad,end_tad in zip(tads.start.values.tolist(),tads.end.values.tolist()):
        if start < start_tad:
            #if start <= previous_end:
            #    return int(previous_start)
            #else: return -1
            return int(previous_start)
        else:
            previous_start = start_tad
            previous_end = end_tad
            continue
    return previous_start
def find_s129_TAD_func_gene(row):
    start = row['start']
    tads = tads_s129_df.query('chrom == "{}"'.format(row['chrom']))
    previous_start = 3000000
    previous_end = 10000000
    for start_tad,end_tad in zip(tads.start.values.tolist(),tads.end.values.tolist()):
        if start < start_tad:
            #if start <= previous_end:
            #    return int(previous_start)
            #else: return -1
            return int(previous_start)
        else:
            previous_start = start_tad
            previous_end = end_tad
            continue
    return previous_start
def find_TAD_func_gene(row):
    start = row['start']
    tads = tads_df.query('chrom == "{}"'.format(row['chrom']))
    previous_start = 3000000
    previous_end = 10000000
    for start_tad,end_tad in zip(tads.start.values.tolist(),tads.end.values.tolist()):
        if start < start_tad:
            #if start <= previous_end:
            #    return int(previous_start)
            #else: return -1
            return int(previous_start)
        else:
            previous_start = start_tad
            previous_end = end_tad
            continue
    return previous_start

# HERE I START TO USE TSS_sTART INSTEAD OF TRANSCRIPT START

# add TAD information to the master table.

master_df['ASE'] = master_df.apply(allelic_func, axis=1)
print ("ASE applied")
master_df["TAD_start_cast"] = master_df.apply(find_cast_TAD_func_gene,axis=1)
print ("TAD1 applied")
master_df["TAD_start_s129"] = master_df.apply(find_s129_TAD_func_gene,axis=1)
print ("TAD2 applied")
master_df["TAD_start"] = master_df.apply(find_TAD_func_gene,axis=1)
print ("TAD3 applied")

#%%
# TAD size distribution depending on number of genes inside the TADs.

aux_df = master_df_aux.query('TPM_transcript >= 1')[["chrom","start","end"]]
aux_pb = pb.BedTool.from_dataframe(aux_df)
tads_pb = pb.BedTool.from_dataframe(tads_df)
inter = tads_pb.intersect(aux_pb,c=True).to_dataframe(names=["chrom","start","end","gene_count"])
inter = inter[inter.chrom.isin(chrs_)]
inter["TAD_size"] = inter.end - inter.start
print (inter)
tad_size_list = []
medians = [0]
for i in range(90):
    aux_df = inter.query('gene_count == {}'.format(i))
    lst = aux_df.TAD_size.values.tolist()
    medians.append(np.median(lst))
    tad_size_list.append(lst)

print (len(medians))
print (len(tad_size_list))
fig,ax = plt.subplots(figsize=(20,6))
ax.boxplot(tad_size_list,showfliers=False)
ax.plot(medians,color="orange")
ax.set_xlabel("Number of exp genes in TADs")
ax.set_ylabel("TAD sizes in bps")
fig.savefig(save_path+"TAD_size_distribution_depending_on_number_of_genes_inside.pdf",dpi=300,bbox_inches='tight')



#%%
# bubble plot for TADs
# SHown in SI figure 5d 

# Bubble sizes are modified so they show an exponential number of the current value, check lines:
#s = aux_mix_df.n_cast**3
#s = aux_mix_df.n_s129**3
#s = aux_mix_df.n_ase**2
# So, the legend of the plots should be square or cube rooted 

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

chrs_ = [ 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19','chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']
ase_path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"
ase_df = pd.read_csv(ase_path, sep="\t")
ase_df = ase_df[ase_df.chrom.isin(chrs_)]
#ase_df = ase_df.drop_duplicates()
ase_df_ = ase_df.rename(columns={"TSS_start":"start","TSS_end":"end"})
cast_color = "#611163"
s129_color = "#EA8D1F"
both_color = "#7A7A7A"
non_phased_color = "#09AA9E"

cast_genes = ase_df_.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')[["chrom","start_transcript","end_transcript"]]
s129_genes = ase_df_.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')[["chrom","start_transcript","end_transcript"]]
exp_genes = ase_df_.query('TPM_transcript >= 1 and ((log2foldchange < 1 and log2foldchange > -1) or p_adj_sum > 0.05) and number_SNPs >= 1')[["chrom","start_transcript","end_transcript"]]

random_ = False
for_main_figure = True
if random_:
    aux_ = ase_df_.query('TPM_transcript >= 1 and number_SNPs >= 1')
    aux_ase = aux_.sample(len(cast_genes)+len(s129_genes))
    aux_cast = aux_ase.sample(len(cast_genes))
    aux_s129 = aux_ase.drop(aux_cast.index, axis=0)
    aux_cast_pb = pb.BedTool.from_dataframe(aux_cast)
    aux_s129_pb = pb.BedTool.from_dataframe(aux_s129)
    aux_exp_pb = pb.BedTool.from_dataframe(exp_genes)
else:
    aux_cast_pb = pb.BedTool.from_dataframe(cast_genes)
    aux_s129_pb = pb.BedTool.from_dataframe(s129_genes)
    aux_exp_pb = pb.BedTool.from_dataframe(exp_genes)


if for_main_figure:
    #fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(6,12),sharex=True)
    fig,axs = plt.subplots(3,2,figsize=(8,12),sharex=True)
else:
    #fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(12,2,figsize=(6,12),sharex=True,sharey=True)
    fig,axs = plt.subplots(3,2,figsize=(8,12),sharex=True,sharey=True)
#first cast
aux_tads = []
label_=["cast","s129"]
scs = []
handles = []
labels = []
hist_values = []
for i,tads_ in enumerate([tads_cast_df,tads_s129_df]):
    ax1 = axs[0][i]
    ax2 = axs[1][i]
    ax3 = axs[2][i]
    tads_pb = pb.BedTool.from_dataframe(tads_)
    tads_pb = tads_pb.intersect(aux_cast_pb,c=True)
    tads_pb = tads_pb.intersect(aux_s129_pb,c=True)
    tads_pb = tads_pb.intersect(aux_exp_pb,c=True)
    tads_df = tads_pb.to_dataframe(names=["chrom","start","end","n_cast","n_s129","n_exp"])
    tads_df["TADsize"] = tads_df.end - tads_df.start
    aux_tads.append(tads_df)
    aux_cast_df = tads_df.query('n_cast != 0 and n_s129 == 0')
    aux_s129_df = tads_df.query('n_cast == 0 and n_s129 != 0')
    aux_mix_df = tads_df.query('n_cast != 0 and n_s129 != 0')
    aux_mix_df["n_ase"] = aux_mix_df.n_cast + aux_mix_df.n_s129
    sc = ax1.scatter(
        x = aux_cast_df.TADsize,
        y = aux_cast_df.n_exp,
        s = aux_cast_df.n_cast**3,
        c = cast_color,
        alpha=0.5, 
        edgecolors="white", 
        linewidth=2
    )
    handle, label = sc.legend_elements(prop="sizes", alpha=0.5)
    legend2 = ax1.legend(handle, label, loc="upper right", title="Sizes")
    scs.append(sc)
    handles.append(handle)
    labels.append(label)
    print ("CAST genes per {} TAD: {}".format(label_[i],sum(aux_cast_df.n_cast)/len(aux_cast_df))) 
    print ("Number of TADs: {}".format(len(aux_cast_df)))
    ax1.text(4000000,20,"CAST genes per {} TAD: {:0.2f}".format(label_[i],sum(aux_cast_df.n_cast)/len(aux_cast_df)))
    sc = ax2.scatter(
        x = aux_s129_df.TADsize,
        y = aux_s129_df.n_exp,
        s = aux_s129_df.n_s129**3,
        c = s129_color,
        alpha=0.5, 
        edgecolors="white", 
        linewidth=2
    )
    handle, label = sc.legend_elements(prop="sizes", alpha=0.5)
    #legend2 = ax2.legend(handle, label, loc="upper right", title="Sizes")
    scs.append(sc)
    handles.append(handle)
    labels.append(label)
    print ("S129 genes per {} TAD: {}".format(label_[i],sum(aux_s129_df.n_s129)/len(aux_s129_df))) 
    print ("Number of TADs: {}".format(len(aux_s129_df)))
    ax2.text(4000000,20,"S129 genes per {} TAD: {:0.2f}".format(label_[i],sum(aux_s129_df.n_s129)/len(aux_s129_df)))
    sc = ax3.scatter(
        x = aux_mix_df.TADsize,
        y = aux_mix_df.n_exp,
        s = aux_mix_df.n_ase**2,
        c = both_color,
        alpha=0.5, 
        edgecolors="white", 
        linewidth=2
    )
    handle, label = sc.legend_elements(prop="sizes", alpha=0.5)
    legend2 = ax3.legend(handle, label, loc="upper right", title="Sizes")
    scs.append(sc)
    handles.append(handle)
    labels.append(label)
    print ("ASE genes per {} TAD: {}".format(label_[i],sum(aux_mix_df.n_ase)/len(aux_mix_df))) 
    print ("Number of TADs: {}".format(len(aux_mix_df)))
    ax3.text(4000000,20,"ASE genes per {} TAD: {:0.2f}".format(label_[i],sum(aux_mix_df.n_ase)/len(aux_mix_df)))
    print (max(aux_cast_df.n_cast))
    if i == 0:
        hist_values.append(aux_cast_df.n_cast.values.tolist())
    elif i == 1:
        hist_values.append(aux_s129_df.n_s129.values.tolist())
#plt.legend(handles,labels)
#plt.xscale('log')
plt.xlabel("TAD size")
plt.ylabel("Number of Biallelic genes")
plt.title("TAD overview. Random: {}".format(random_))
#plt.ylim(0,50000)
#plt.xlim(30, 75)
plt.savefig(save_path+"bubble_plots_{}_random_{}_main_figure_{}.pdf".format(label_[i],random_,for_main_figure))
plt.show()



#%%

# permutation test to check if number of genes in TADs is significantly different from the real numbers.
import random
permuted_hist_values_cast = []
permuted_hist_values_s129 = []
cast_permuted_tads = []
s129_permuted_tads = []
cast_permuted_ratio = []
s129_permuted_ratio = []
ase_permuted_ratio_cast = []
ase_permuted_ratio_s129 = []
n_cast_tads_real = 0
n_s129_tads_real = 0
cast_true_ratio_cast = 0
s129_true_ratio_s129 = 0
ase_true_ratio_cast = 0
ase_true_ratio_s129 = 0
get_significance = True
if get_significance:
    random.seed(3)
    permutes = 10000
    label=["cast","s129"]
    for j,tads_ in enumerate([tads_cast_df,tads_s129_df]):
        cast_counter, s129_counter, mix_counter = 0,0,0
        tads_pb = pb.BedTool.from_dataframe(tads_)
        aux_cast_pb = pb.BedTool.from_dataframe(cast_genes)
        aux_s129_pb = pb.BedTool.from_dataframe(s129_genes)
        tads_aux_pb = tads_pb.intersect(aux_cast_pb,c=True)
        tads_aux_pb = tads_aux_pb.intersect(aux_s129_pb,c=True)
        tads_df = tads_aux_pb.to_dataframe(names=["chrom","start","end","n_cast","n_s129"])
        aux_cast_df = tads_df.query('n_cast != 0 and n_s129 == 0')
        aux_s129_df = tads_df.query('n_cast == 0 and n_s129 != 0')
        aux_mix_df = tads_df.query('n_cast != 0 and n_s129 != 0')
        aux_mix_df["n_ase"] = aux_mix_df.n_cast + aux_mix_df.n_s129
        cast_true_ratio = sum(aux_cast_df.n_cast)/len(aux_cast_df)
        s129_true_ratio = sum(aux_s129_df.n_s129)/len(aux_s129_df)
        mix_true_ratio = sum(aux_mix_df.n_ase)/len(aux_mix_df)
        aux_ = ase_df_.query('TPM_transcript >= 1 and number_SNPs >= 1')
        if j == 0:
            cast_true_ratio_cast = cast_true_ratio
            ase_true_ratio_cast = mix_true_ratio
            n_cast_tads_real = len(aux_cast_df)
        elif j == 1:
            s129_true_ratio_s129 = s129_true_ratio
            ase_true_ratio_s129 = mix_true_ratio
            n_s129_tads_real = len(aux_s129_df)

        for i in range(permutes):
            if i % 100 == 0:
                print (i)
            aux_ase = aux_.sample(len(cast_genes)+len(s129_genes))
            aux_cast = aux_ase.sample(len(cast_genes))
            aux_s129 = aux_ase.drop(aux_cast.index, axis=0)
            aux_cast_pb = pb.BedTool.from_dataframe(aux_cast)
            aux_s129_pb = pb.BedTool.from_dataframe(aux_s129)
            tads_aux_pb = tads_pb.intersect(aux_cast_pb,c=True)
            tads_aux_pb = tads_aux_pb.intersect(aux_s129_pb,c=True)
            tads_df = tads_aux_pb.to_dataframe(names=["chrom","start","end","n_cast","n_s129"])
            aux_cast_df = tads_df.query('n_cast != 0 and n_s129 == 0')
            aux_s129_df = tads_df.query('n_cast == 0 and n_s129 != 0')
            aux_mix_df = tads_df.query('n_cast != 0 and n_s129 != 0')
            aux_mix_df["n_ase"] = aux_mix_df["n_cast"] + aux_mix_df["n_s129"]
            cast_ratio = sum(aux_cast_df.n_cast)/len(aux_cast_df)
            s129_ratio = sum(aux_s129_df.n_s129)/len(aux_s129_df)
            mix_ratio = sum(aux_mix_df.n_ase)/len(aux_mix_df)
            if cast_ratio > cast_true_ratio:
                cast_counter += 1
            if s129_ratio > s129_true_ratio:
                s129_counter += 1
            if mix_ratio > mix_true_ratio:
                mix_counter += 1
            if j == 0:
                permuted_hist_values_cast.append(aux_cast_df.n_cast.values.tolist())
                cast_permuted_tads.append(len(aux_cast_df))
                cast_permuted_ratio.append(cast_ratio)
                ase_permuted_ratio_cast.append(mix_ratio)
            elif j == 1:
                permuted_hist_values_s129.append(aux_s129_df.n_s129.values.tolist())
                s129_permuted_tads.append(len(aux_s129_df))
                s129_permuted_ratio.append(s129_ratio)
                ase_permuted_ratio_s129.append(mix_ratio)
        print ("p value cast: {}/{} = {}".format(cast_counter,permutes,cast_counter/permutes))
        print ("p value s129: {}/{} = {}".format(s129_counter,permutes,s129_counter/permutes))
        print ("p value mix: {}/{} = {}".format(mix_counter,permutes,mix_counter/permutes))


fig, ax = plt.subplots(figsize=(4,4))
vector = [[n_cast_tads_real],cast_permuted_tads,[n_s129_tads_real],s129_permuted_tads]
parts = ax.violinplot(vector,showmeans=False, showmedians=True,showextrema=True)
for pc in parts['bodies']:
    pc.set_alpha(1)
    pc.set_edgecolor([0,0,0,1])
    pc.set_facecolor('#FFFFFF')
    pc.set_facecolor([1,1,1,0])
for partname in ('','cmedians'):
    try:
        vp = parts[partname]
        vp.set_edgecolor("black")
        #vp.set_linewidth(1)
    except:
        continue
colors = [cast_color,cast_color,s129_color,s129_color]
for pos,v in enumerate(vector): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax.scatter(xv, v,s=1,alpha=1,color=colors[pos])
labels = ["CAST TADs","Permuted CAST TADs","S129 TADs", "Permuted S129 TADs"]
set_axis_style(ax, labels)
ax.set_ylabel("number of TADs")
plt.savefig(save_path+"tads_clustered.pdf",dpi=300)

fig, ax = plt.subplots(figsize=(4,4))
vector = [[cast_true_ratio_cast],cast_permuted_ratio,[s129_true_ratio],s129_permuted_ratio,
          [ase_true_ratio_cast],ase_permuted_ratio_cast,[ase_true_ratio_s129],ase_permuted_ratio_s129]
parts = ax.violinplot(vector,showmeans=False, showmedians=True,showextrema=True)
for pc in parts['bodies']:
    pc.set_alpha(1)
    pc.set_edgecolor([0,0,0,1])
    pc.set_facecolor('#FFFFFF')
    pc.set_facecolor([1,1,1,0])
for partname in ('','cmedians'):
    try:
        vp = parts[partname]
        vp.set_edgecolor("black")
        #vp.set_linewidth(1)
    except:
        continue
colors = [cast_color,cast_color,s129_color,s129_color,both_color,both_color,both_color,both_color]
for pos,v in enumerate(vector): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax.scatter(xv, v,s=1,alpha=1,color=colors[pos])
labels = ["CAST","Permuted CAST","S129", "Permuted S129","ASE","Permuted ASE CAST","ASE","Permuted ASE S129"]
set_axis_style(ax, labels)
ax.set_ylabel("genes per TAD")
plt.savefig(save_path+"tads_ratios_clustered.pdf",dpi=300)
fig, ax = plt.subplots(figsize=(4,4))
vector = [[cast_true_ratio_cast],cast_permuted_ratio,[s129_true_ratio],s129_permuted_ratio,
          [ase_true_ratio_cast],ase_permuted_ratio_cast,[ase_true_ratio_s129],ase_permuted_ratio_s129]
parts = ax.violinplot(vector,showmeans=False, showmedians=True,showextrema=True)
for pc in parts['bodies']:
    pc.set_alpha(1)
    pc.set_edgecolor([0,0,0,1])
    pc.set_facecolor('#FFFFFF')
    pc.set_facecolor([1,1,1,0])
for partname in ('','cmedians'):
    try:
        vp = parts[partname]
        vp.set_edgecolor("black")
        #vp.set_linewidth(1)
    except:
        continue

labels = ["CAST","Permuted CAST","S129", "Permuted S129","ASE","Permuted ASE CAST","ASE","Permuted ASE S129"]
set_axis_style(ax, labels)
ax.set_ylabel("genes per TAD")
plt.savefig(save_path+"tads_ratios_clustered_2.pdf",dpi=300)



df_cast = pd.DataFrame(permuted_hist_values_cast)
df_s129 = pd.DataFrame(permuted_hist_values_s129)
  
        
#%%
# histogram to show significance 

fig,([ax1,ax2]) = plt.subplots(1,2,figsize=(8,4),sharey=True,sharex=True)
#ax1.hist(hist_values[0],bins=10)
#ax2.hist(hist_values[1],bins=10)
#ax3.hist(permuted_hist_values[0],bins=10)
#ax4.hist(permuted_hist_values[1],bins=10)

width=0.4
size=0.4
alpha=1

b = np.bincount(hist_values[0])[1:]
ind = np.arange(len(b))
ax1.bar(ind-width/2., b,color=cast_color,alpha=alpha,width=size)
#ax1.plot(ind-width/2., b,color=cast_color)
ax1.axvline(1.77,color=cast_color)

b = np.bincount(hist_values[1])[1:]
ind = np.arange(len(b))
ax2.bar(ind-width/2., b,color=s129_color,alpha=alpha,width=size)
ax2.axvline(1.57,color=s129_color)

b = np.bincount(permuted_hist_values_cast[0])[1:]
ind = np.arange(len(b))
ax1.bar(ind+width/2., b,color="grey",alpha=alpha,width=size)
#ax1.plot(ind-width/2., b,color="grey")
ax1.axvline(1.48,color="grey")

b = np.bincount(permuted_hist_values_s129[0])[1:]
ind = np.arange(len(b))
ax2.bar(ind+width/2., b,color="grey",alpha=alpha,width=size)
ax2.axvline(1.34,color="grey")

ax1.set_xticks(np.arange(9))
ax1.set_xticklabels(np.arange(10)[1:])
fig.savefig(save_path+"histogram_significance_ASEs_in_TADs.pdf",dpi=300)

import scipy
print (scipy.stats.ks_2samp(hist_values[0],permuted_hist_values_cast[0],alternative='less'))
print (scipy.stats.ks_2samp(hist_values[1],permuted_hist_values_s129[0],alternative='less'))


fig,([ax1,ax2]) = plt.subplots(1,2,figsize=(8,4),sharey=True,sharex=True)
b = np.bincount(hist_values[0])
ind = np.arange(len(b))
ax1.plot(ind,np.cumsum(b/sum(b)),color=cast_color)
b = np.bincount(permuted_hist_values_cast[0])
ind = np.arange(len(b))
ax1.plot(ind,np.cumsum(b/sum(b)),color="grey")
b = np.bincount(hist_values[1])
ind = np.arange(len(b))
ax2.plot(ind,np.cumsum(b/sum(b)),color=s129_color)
b = np.bincount(permuted_hist_values_s129[0])
ind = np.arange(len(b))
ax2.plot(ind,np.cumsum(b/sum(b)),color="grey")



#%%
# Figure 3c and SI Figure 5c upset plots.
# depending on the variable, we can plot CAST TADs or S129 TADs.

from upsetplot import UpSet
upset_df = pd.DataFrame(columns=["no","exp","CAST","S129"])
counter = 0
#tad_type = "TAD_start_s129"
tad_type = "TAD_start_cast"
#tad_type = "TAD_start"
for chr_ in chrs_:
    aux_df = master_df.query("chrom == '{}' and ASE != 'no'".format(chr_))[["ASE",tad_type]]
    tads = set(aux_df[tad_type].values)
    for tad in tads:
        genes = set(aux_df.query('{} == {}'.format(tad_type,tad)).ASE.values)
        for gene in genes:
            if gene == "exp_no_SNP":
                gene = "exp"
            upset_df.at[counter,gene] = True
        counter += 1 
upset_df = upset_df.fillna(0)
binary_categories_df = upset_df[["exp","CAST","S129"]]
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '>')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))
facecolor = "black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True)
upset.style_subsets(present="exp>",
                    facecolor="darkgreen",
                    label="exp")
upset.style_subsets(present="CAST>",
                    facecolor=cast_color,
                    label="CAST")
upset.style_subsets(present="S129>",
                    facecolor=s129_color,
                    label="S129")
upset.style_subsets(present=["CAST>", "S129>"],
                    facecolor="teal",
                    label="ASE")

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(tad_type)
plt.plot()
fig.savefig("{}upset_plot_tads_genes_{}_no_not_expressed_{}.pdf".format(save_path,tad_type,date),dpi=300)


#first convert df into upset plot like df
upset_df = pd.DataFrame(columns=["no","exp","CAST","S129"])
counter = 0
for chr_ in chrs_:
    aux_df = master_df.query("chrom == '{}'".format(chr_))[["ASE",tad_type]]
    tads = set(aux_df[tad_type].values)
    for tad in tads:
        genes = set(aux_df.query('{} == {}'.format(tad_type,tad)).ASE.values)
        for gene in genes:
            if gene == "exp_no_SNP":
                gene = "exp"
            upset_df.at[counter,gene] = True
        counter += 1 
upset_df = upset_df.fillna(0)

from upsetplot import UpSet
binary_categories_df = upset_df[["no","exp","CAST","S129"]]
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '>')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))
facecolor = "black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor,
                        show_counts=True)
upset.style_subsets(present="no>",
                    facecolor="grey",
                    label="not expressed")
upset.style_subsets(present="exp>",
                    facecolor="darkgreen",
                    label="exp")
upset.style_subsets(present="CAST>",
                    facecolor=cast_color,
                    label="CAST")
upset.style_subsets(present="S129>",
                    facecolor=s129_color,
                    label="S129")
upset.style_subsets(present=["CAST>", "S129>"],
                    facecolor="teal",
                    label="ASE")

#first convert df into upset plot like df
upset_df = pd.DataFrame(columns=["no","exp","CAST","S129"])
counter = 0
for chr_ in chrs_:
    aux_df = master_df.query("chrom == '{}'".format(chr_))[["ASE",tad_type]]
    tads = set(aux_df[tad_type].values)
    for tad in tads:
        genes = set(aux_df.query('{} == {}'.format(tad_type,tad)).ASE.values)
        for gene in genes:
            if gene == "exp_no_SNP":
                gene = "exp"
            upset_df.at[counter,gene] = True
        counter += 1 
upset_df = upset_df.fillna(0)

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(tad_type)
fig.savefig("{}upset_plot_tads_genes_{}_{}.pdf".format(save_path,tad_type,date),dpi=300)

#now SUP with no SNP exp
upset_df = pd.DataFrame(columns=["no","exp","exp_no_SNP","CAST","S129"])
counter = 0
for chr_ in chrs_:
    aux_df = master_df.query("chrom == '{}'".format(chr_))[["ASE",tad_type]]
    tads = set(aux_df[tad_type].values)
    for tad in tads:
        genes = set(aux_df.query('{} == {}'.format(tad_type,tad)).ASE.values)
        for gene in genes:
            upset_df.at[counter,gene] = True
        counter += 1 
upset_df = upset_df.fillna(0)

binary_categories_df = upset_df[["no","exp","exp_no_SNP","CAST","S129"]]
binary_categories_df = binary_categories_df.rename(columns=lambda x: x + '>')
upset_final_df = pd.concat([upset_df, binary_categories_df],axis=1)
upset_final_df = upset_final_df.set_index(list(binary_categories_df.columns))
#facecolor=["darkgreen","grey","darkgreen","teal",cast_color,cast_color,s129_color,s129_color,"teal","darkgreen","darkgreen",cast_color,s129_color,cast_color,cast_color,"darkgreen","teal",s129_color,s129_color,"darkgreen","teal",s129_color,"teal",cast_color]
facecolor = "black"
upset = UpSet(upset_final_df, subset_size='count', 
                        intersection_plot_elements=3,
                        sort_by="cardinality",
                        facecolor=facecolor)
upset.style_subsets(present="no>",
                    facecolor="grey",
                    label="not expressed")
upset.style_subsets(present="exp_no_SNP>",
                    facecolor="lightgreen",
                    label="exp_no_SNP")
upset.style_subsets(present="exp>",
                    facecolor="darkgreen",
                    label="exp")
upset.style_subsets(present="CAST>",
                    facecolor=cast_color,
                    label="CAST")
upset.style_subsets(present="S129>",
                    facecolor=s129_color,
                    label="S129")
upset.style_subsets(present=["CAST>", "S129>"],
                    facecolor="teal",
                    label="ASE")

fig = plt.figure(figsize=(20,40))
upset.plot(fig=fig)
plt.title(tad_type)
fig.savefig("{}upset_plot_tads_genes_supp_{}_{}.pdf".format(save_path,tad_type,date),dpi=300)


#%%

import scipy.stats

# Boxplots steaming from Figure3c and SI Figure 5c.
# These figures are Figure3d, 3e, SI 5f and 5g

if tad_type == "TAD_start_cast":
    tads_df = tads_cast_df[tads_cast_df["chrom"].isin(chrs_)]
    tad_type = "cast_TADs"
elif tad_type == "TAD_start_s129":
    tads_df = tads_s129_df[tads_s129_df["chrom"].isin(chrs_)]
    tad_type = "s129_TADs"
else:
    tads_df = tads_df[tads_df["chrom"].isin(chrs_)]
    tad_type = "unphased_TADs"

tads_pb = pb.BedTool.from_dataframe(tads_df)


inter = tads_pb.intersect(cast_pb,c=True)
inter = inter.intersect(s129_pb,c=True)
inter = inter.intersect(exp_pb,c=True)
inter = inter.intersect(exp_no_SNP_pb,c=True)
inter = inter.intersect(not_exp_pb,c=True)
inter_df = inter.to_dataframe(names = ["n_cast","n_s129","n_exp_SNP","n_exp_no_SNP","n_not_exp"])
inter_df["n_total_exp_genes"] = inter_df.n_cast +inter_df.n_s129 +inter_df.n_exp_SNP +inter_df.n_exp_no_SNP 
inter_df = inter_df.reset_index()
inter_df["TADsize"] = inter_df.level_2 - inter_df.level_1 

# depending on the type of TADs. We will choose 2 or more genes
not_exp_TAD_df = inter_df.query("n_not_exp != 0 and n_cast == 0 and n_s129 == 0 and n_exp_SNP == 0 and n_exp_no_SNP == 0").reset_index()
exp_SNP_TAD_df = inter_df.query("n_cast == 0 and n_s129 == 0 and n_exp_SNP != 0 and n_exp_no_SNP == 0").reset_index()
exp_wo_SNP_TAD_df = inter_df.query("n_cast == 0 and n_s129 == 0 and n_exp_SNP == 0 and n_exp_no_SNP != 0").reset_index()
ase_TAD_df = inter_df.query("n_cast != 0 and n_s129 != 0 ").reset_index()
cast_TAD_df = inter_df.query("n_cast != 0 and n_s129 == 0").reset_index()
s129_TAD_df = inter_df.query("n_cast == 0 and n_s129 != 0").reset_index()
#not_exp_TAD_df = inter_df.query("n_not_exp != 0 and n_cast == 0 and n_s129 == 0 and n_exp_SNP == 0 and n_exp_no_SNP == 0").reset_index()
#exp_SNP_TAD_df = inter_df.query("n_cast == 0 and n_s129 == 0 and n_exp_SNP >= 2 and n_exp_no_SNP == 0").reset_index()
#exp_wo_SNP_TAD_df = inter_df.query("n_cast == 0 and n_s129 == 0 and n_exp_SNP == 0 and n_exp_no_SNP >= 2").reset_index()
#ase_TAD_df = inter_df.query("n_cast >= 2 and n_s129 >= 2 ").reset_index()
#cast_TAD_df = inter_df.query("n_cast >= 2 and n_s129 == 0").reset_index()
#s129_TAD_df = inter_df.query("n_cast == 0 and n_s129 >= 2").reset_index()
not_exp_TAD_df = not_exp_TAD_df.iloc[: , 1:]
exp_SNP_TAD_df = exp_SNP_TAD_df.iloc[: , 1:]
exp_wo_SNP_TAD_df = exp_wo_SNP_TAD_df.iloc[: , 1:]
ase_TAD_df = ase_TAD_df.iloc[: , 1:]
cast_TAD_df = cast_TAD_df.iloc[: , 1:]
s129_TAD_df = s129_TAD_df.iloc[: , 1:]
not_exp_TAD_df.rename(columns={"level_0":"chrom","level_1":"start","level_2":"end"},inplace=True)
exp_SNP_TAD_df.rename(columns={"level_0":"chrom","level_1":"start","level_2":"end"},inplace=True)
exp_wo_SNP_TAD_df.rename(columns={"level_0":"chrom","level_1":"start","level_2":"end"},inplace=True)
ase_TAD_df.rename(columns={"level_0":"chrom","level_1":"start","level_2":"end"},inplace=True)
cast_TAD_df.rename(columns={"level_0":"chrom","level_1":"start","level_2":"end"},inplace=True)
s129_TAD_df.rename(columns={"level_0":"chrom","level_1":"start","level_2":"end"},inplace=True)
print (not_exp_TAD_df)
print (exp_SNP_TAD_df)
print (exp_wo_SNP_TAD_df)
print (ase_TAD_df)
print (cast_TAD_df)
print (s129_TAD_df)

labels = ["exp_w_SNP","not_exp","cast","s129","ase","exp_no_SNP"]
fig,ax = plt.subplots()
plot_values = [exp_SNP_TAD_df.TADsize.values.tolist(),not_exp_TAD_df.TADsize.values.tolist(),cast_TAD_df.TADsize.values.tolist(),s129_TAD_df.TADsize.values.tolist(),ase_TAD_df.TADsize.values.tolist(),exp_wo_SNP_TAD_df.TADsize.values.tolist()]

bp1 = ax.boxplot(plot_values,patch_artist=True ,showfliers=False,widths=0.5)
ax.set_xticklabels(labels)
colors = ["darkgreen","grey",cast_color,s129_color,non_phased_color,"limegreen"]
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(alpha=0.5)

for pos,v in enumerate(plot_values): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax.scatter(xv, v,s=1,alpha=1,color="black")
plt.setp(bp1["medians"],color="white")
ax.set_xlabel("TAD groups")
ax.set_ylabel("TAD size")
plt.savefig("{}TAD_size_of_TAD_groups_{}.pdf".format(save_path,tad_type),dpi=300)

fig,ax = plt.subplots()
plot_values = [exp_SNP_TAD_df.n_total_exp_genes.values.tolist(),not_exp_TAD_df.n_total_exp_genes.values.tolist(),cast_TAD_df.n_total_exp_genes.values.tolist(),s129_TAD_df.n_total_exp_genes.values.tolist(),ase_TAD_df.n_total_exp_genes.values.tolist(),exp_wo_SNP_TAD_df.n_total_exp_genes.values.tolist()]

bp1 = ax.boxplot(plot_values,patch_artist=True ,showfliers=False,widths=0.5)
ax.set_xticklabels(labels)
colors = ["darkgreen","grey",cast_color,s129_color,non_phased_color,"limegreen"]
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(alpha=0.5)

for pos,v in enumerate(plot_values): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax.scatter(xv, v,s=1,alpha=1,color="black")
plt.setp(bp1["medians"],color="white")
ax.set_xlabel("TAD groups")
ax.set_ylabel("Number of exp. genes")
plt.savefig("{}exp_genes_in_TAD_groups_{}.pdf".format(save_path,tad_type),dpi=300)

#normalized by tad size
fig,ax = plt.subplots()
plot_values = [[a/b for a,b in zip(exp_SNP_TAD_df.n_total_exp_genes.values.tolist(),exp_SNP_TAD_df.TADsize.values.tolist())]
            ,[a/b for a,b in zip(not_exp_TAD_df.n_total_exp_genes.values.tolist(),not_exp_TAD_df.TADsize.values.tolist())]
            ,[a/b for a,b in zip(cast_TAD_df.n_total_exp_genes.values.tolist(),cast_TAD_df.TADsize.values.tolist())]
            ,[a/b for a,b in zip(s129_TAD_df.n_total_exp_genes.values.tolist(),s129_TAD_df.TADsize.values.tolist())]
            ,[a/b for a,b in zip(ase_TAD_df.n_total_exp_genes.values.tolist(),ase_TAD_df.TADsize.values.tolist())]
            ,[a/b for a,b in zip(exp_wo_SNP_TAD_df.n_total_exp_genes.values.tolist(),exp_wo_SNP_TAD_df.TADsize.values.tolist())]
            ]

bp1 = ax.boxplot(plot_values,patch_artist=True ,showfliers=False,widths=0.5)
ax.set_xticklabels(labels)
colors = ["darkgreen","grey",cast_color,s129_color,non_phased_color,"limegreen"]
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(alpha=0.5)

for pos,v in enumerate(plot_values): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax.scatter(xv, v,s=1,alpha=1,color="black")
plt.setp(bp1["medians"],color="white")
ax.set_xlabel("TAD groups")
ax.set_ylabel("Number of exp. genes normalized by TAD size")
plt.savefig("{}exp_genes_in_TAD_groups_normalized_{}.pdf".format(save_path,tad_type),dpi=300)

no_exp_TAD_pb = pb.BedTool.from_dataframe(not_exp_TAD_df).sort()
exp_SNP_TAD_pb = pb.BedTool.from_dataframe(exp_SNP_TAD_df).sort()
exp_wo_SNP_TAD_pb = pb.BedTool.from_dataframe(exp_wo_SNP_TAD_df).sort()
ase_TAD_pb = pb.BedTool.from_dataframe(ase_TAD_df).sort()
cast_TAD_pb = pb.BedTool.from_dataframe(cast_TAD_df).sort()
s129_TAD_pb = pb.BedTool.from_dataframe(s129_TAD_df).sort()


#get number of H3K27me3 peaks in these different groups
peaks_path = root+"F123.K27me3.BCP_peaks.HM_mode.bed"
h3k27me3_df = pd.read_csv(peaks_path,sep="\t",names=["chrom","start","end","size","score"])
h3k27me3_df = h3k27me3_df[h3k27me3_df["chrom"].isin(chrs_)]
h3k27me3_pb = pb.BedTool.from_dataframe(h3k27me3_df).sort()

pbs = [exp_SNP_TAD_pb,no_exp_TAD_pb,cast_TAD_pb,s129_TAD_pb,ase_TAD_pb,exp_wo_SNP_TAD_pb]
plot_values = []
for x in pbs:
    aux = x.intersect(h3k27me3_pb,c=True,wa=True).to_dataframe(names=["chrom"	,"start",	"end"	,"n_cast",	"n_s129",	"n_exp_SNP"	,"n_exp_no_SNP"	,"n_not_exp"	,"n_total_exp_genes"	,"TADsize","n_h3k27me3_peaks"])
    plot_values.append([a/b for a,b in zip(aux.n_h3k27me3_peaks.values.tolist(),aux.TADsize.values.tolist())]) #these are prc2 peaks

fig,ax = plt.subplots()
bp1 = ax.boxplot(plot_values,patch_artist=True ,showfliers=False,widths=0.5)
ax.set_xticklabels(labels)
colors = ["darkgreen","grey",cast_color,s129_color,non_phased_color,"limegreen"]
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(alpha=0.5)

for pos,v in enumerate(plot_values): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax.scatter(xv, v,s=1,alpha=1,color="black")
plt.setp(bp1["medians"],color="white")
ax.set_xlabel("TAD groups")
ax.set_ylabel("number of H3K27me3 peaks normalized by TADsize")
#ax.set_ylabel("number of H3K27me3 peaks normalized by number of genes")
plt.savefig("{}PRC2_peaks_in_TAD_groups_{}.pdf".format(save_path,tad_type),dpi=300)

print (scipy.stats.ttest_ind(plot_values[0],plot_values[1]))
print (scipy.stats.ttest_ind(plot_values[0],plot_values[2]))
print (scipy.stats.ttest_ind(plot_values[0],plot_values[3]))
print (scipy.stats.ttest_ind(plot_values[0],plot_values[4]))
print (scipy.stats.ttest_ind(plot_values[0],plot_values[5]))
print (scipy.stats.ttest_ind(plot_values[2],plot_values[3]))
print (scipy.stats.ttest_ind(plot_values[2],plot_values[4]))
print (scipy.stats.ttest_ind(plot_values[3],plot_values[4]))

#%%

# now get the WDF values and the differential WDF for those groups
# Figure 3e and SI Figure 5g

seg_cast_path = "{}Curated_F123_1122_ALL_as_3NP_CAST_at50000.passed_qc_fc5_cw6_s11.table".format(root)
seg_s129_path = "{}Curated_F123_1122_ALL_as_3NP_S129_at50000.passed_qc_fc5_cw6_s11.table".format(root)
chromosomes = chrs_

def get_WDF(x):
    return x.sum()/n_NPs

seg_table_cast = pd.read_csv(seg_cast_path,sep="\t")
seg_table_cast.set_index(['chrom','start','stop'], inplace = True)
seg_table_cast = seg_table_cast.loc[chromosomes]
n_NPs = seg_table_cast.shape[1]
seg_table_cast["cast_WDF"] = seg_table_cast.apply(lambda x: get_WDF(x),axis=1)
print (seg_table_cast)

plt.hist(seg_table_cast.cast_WDF)

seg_table_s129 = pd.read_csv(seg_s129_path,sep="\t")
seg_table_s129.set_index(['chrom','start','stop'], inplace = True)
seg_table_s129 = seg_table_s129.loc[chromosomes]
n_NPs = seg_table_s129.shape[1]
seg_table_s129["s129_WDF"] = seg_table_s129.apply(lambda x: get_WDF(x),axis=1)
seg_table_cast = seg_table_cast["cast_WDF"]
seg_table_s129 = seg_table_s129["s129_WDF"]

seg_table_cast = seg_table_cast.reset_index()
seg_table_s129 = seg_table_s129.reset_index()
#seg_table = pd.merge(seg_table_cast,seg_table_s129,right_index=True,left_index=True)
seg_table = pd.merge(seg_table_cast,seg_table_s129)
seg_table = seg_table.reset_index()
print (seg_table)
seg_table["diff_WDF"] = seg_table.cast_WDF-seg_table.s129_WDF
seg_table = seg_table[["chrom","start","stop","diff_WDF"]]
wdf_cast_pb = pb.BedTool.from_dataframe(seg_table_cast)
wdf_s129_pb = pb.BedTool.from_dataframe(seg_table_s129)
wdf_diff_pb = pb.BedTool.from_dataframe(seg_table)

labels = ["exp_w_SNP","not_exp","cast","s129","ase","exp_no_SNP"]

exp_SNP_TAD_df = exp_SNP_TAD_df.iloc[:,[0,1,2]]
exp_SNP_TAD_df.columns = ["chrom","start","end"]
exp_SNP_TAD_pb = pb.BedTool.from_dataframe(exp_SNP_TAD_df)
not_exp_TAD_df = not_exp_TAD_df.iloc[:,[0,1,2]]
not_exp_TAD_df.columns = ["chrom","start","end"]
no_exp_TAD_pb = pb.BedTool.from_dataframe(not_exp_TAD_df)
cast_TAD_df = cast_TAD_df.iloc[:,[0,1,2]]
cast_TAD_df.columns = ["chrom","start","end"]
cast_TAD_pb = pb.BedTool.from_dataframe(cast_TAD_df)
s129_TAD_df = s129_TAD_df.iloc[:,[0,1,2]]
s129_TAD_df.columns = ["chrom","start","end"]
s129_TAD_pb = pb.BedTool.from_dataframe(s129_TAD_df)
ase_TAD_df = ase_TAD_df.iloc[:,[0,1,2]]
ase_TAD_df.columns = ["chrom","start","end"]
ase_TAD_pb = pb.BedTool.from_dataframe(ase_TAD_df)
exp_wo_SNP_TAD_df = exp_wo_SNP_TAD_df.iloc[:,[0,1,2]]
exp_wo_SNP_TAD_df.columns = ["chrom","start","end"]
exp_wo_SNP_TAD_pb = pb.BedTool.from_dataframe(exp_wo_SNP_TAD_df)

pbs = [exp_SNP_TAD_pb,no_exp_TAD_pb,cast_TAD_pb,s129_TAD_pb,ase_TAD_pb,exp_wo_SNP_TAD_pb]
plot_values_cast = []
for n,x in enumerate(pbs):
    aux = x.intersect(wdf_cast_pb,wb=True,wa=True).to_dataframe()
    mean_wdf_df = (aux.groupby(['chrom','start','end'])["thickStart"].mean()).values.tolist()
    plot_values_cast.append(mean_wdf_df)
plot_values_s129 = []
for n,x in enumerate(pbs):
    aux = x.intersect(wdf_s129_pb,wb=True,wa=True).to_dataframe()
    mean_wdf_df = (aux.groupby(['chrom','start','end'])["thickStart"].mean()).values.tolist()
    plot_values_s129.append(mean_wdf_df)

fig,(ax1,ax2) = plt.subplots(1,2,sharey=True)
bp1 = ax1.boxplot(plot_values_cast,patch_artist=True ,showfliers=False,widths=0.5)
ax1.set_xticklabels(labels)
colors = ["darkgreen","grey",cast_color,s129_color,non_phased_color,"limegreen"]
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(alpha=0.5)
for pos,v in enumerate(plot_values_cast): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax1.scatter(xv, v,s=1,alpha=1,color="black")
plt.setp(bp1["medians"],color="white")
ax1.set_xlabel("TAD groups")
ax1.set_ylabel("average WDF of TADs")
bp2 = ax2.boxplot(plot_values_s129,patch_artist=True ,showfliers=False,widths=0.5)
ax2.set_xticklabels(labels)
colors = ["darkgreen","grey",cast_color,s129_color,non_phased_color,"limegreen"]
for patch, color in zip(bp2['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(alpha=0.5)
for pos,v in enumerate(plot_values_s129): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax2.scatter(xv, v,s=1,alpha=1,color="black")
plt.setp(bp2["medians"],color="white")
ax2.set_xlabel("TAD groups")
ax2.set_ylabel("average WDF of TADs")
ax1.set_title("CAST WDF")
ax2.set_title("S129 WDF")
plt.savefig("{}/avg_WDF_in_TAD_groups_{}.pdf".format(save_path,tad_type),dpi=300)
plot_values_diff = []
for n,x in enumerate(pbs):
    aux = x.intersect(wdf_diff_pb,wb=True,wa=True).to_dataframe()
    mean_wdf_df = (aux.groupby(['chrom','start','end'])["thickStart"].mean()).values.tolist()
    plot_values_diff.append(mean_wdf_df)
fig,ax1 = plt.subplots(1,1,sharey=True)
bp1 = ax1.boxplot(plot_values_diff,patch_artist=True ,showfliers=False,widths=0.5)
ax1.set_xticklabels(labels)
colors = ["darkgreen","grey",cast_color,s129_color,non_phased_color,"limegreen"]
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)
    patch.set(alpha=0.5)
for pos,v in enumerate(plot_values_diff): 
    xv = np.random.normal(pos+1, 0.05, len(v))
    ax1.scatter(xv, v,s=1,alpha=1,color="black")
plt.setp(bp1["medians"],color="white")
ax1.set_xlabel("TAD groups")
ax1.set_ylabel("average diff WDF of TADs")

print (scipy.stats.ttest_ind(plot_values_diff[0],plot_values_diff[1]))
print (scipy.stats.ttest_ind(plot_values_diff[0],plot_values_diff[2]))
print (scipy.stats.ttest_ind(plot_values_diff[0],plot_values_diff[3]))
print (scipy.stats.ttest_ind(plot_values_diff[0],plot_values_diff[4]))
print (scipy.stats.ttest_ind(plot_values_diff[0],plot_values_diff[5]))
print (scipy.stats.ttest_ind(plot_values_diff[2],plot_values_diff[3]))
print (scipy.stats.ttest_ind(plot_values_diff[2],plot_values_diff[4]))
print (scipy.stats.ttest_ind(plot_values_diff[3],plot_values_diff[4]))

plt.savefig("{}diff_avg_WDF_in_TAD_groups_{}.pdf".format(save_path,tad_type),dpi=300)
