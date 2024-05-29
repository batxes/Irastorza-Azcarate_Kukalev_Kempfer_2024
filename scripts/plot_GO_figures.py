#%%

# Figures showing scatterplots for Gene Ontology
# Figures 2c and 4e


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
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

def plot_fig (x,y,s,c,genes,title,save_title):
    fig,ax = plt.subplots(figsize=(4,4))
    sc=ax.scatter(x,range(len(y)),s=s,c=c,cmap="RdYlBu_r",linewidths=1,edgecolors="black")
    ax.set_yticks(range(len(y)))
    y = [m+"\n"+n for (m,n) in zip(y,genes)]
    ax.set_yticklabels(y)
    ax.set_xlabel("P value")
    handles, labels = sc.legend_elements(prop="sizes", alpha=0.6,num=6)
    legend2 = ax.legend(handles, labels, loc="upper right",title="Number of genes",prop={'size': 10})
    v1 = np.linspace(c.min(), c.max(), 4, endpoint=True)
    print(v1)
    cbar = fig.colorbar(sc,label="% of genes from pathway",shrink=0.5,orientation="horizontal",ticks=v1)
    ax.set_title(title)
    fig.savefig("{}GO_{}_without_fdr.pdf".format(save_path,save_title),bbox_inches="tight")

def plot_fig_fdr (x,y,s,c,genes,title,save_title,fdr):
    fig,ax = plt.subplots(figsize=(4,4))
    ax.axhline(0.05,ls="--",color="black")
    sc=ax.scatter(x,fdr,s=s,c=c,cmap="RdYlBu_r",linewidths=1,edgecolors="black")
    y = [m+"\n"+n for (m,n) in zip(y,genes)]
    plt.gca().invert_yaxis()
    ax.set_xlabel("P value")
    handles, labels = sc.legend_elements(prop="sizes", alpha=0.6,num=6)
    legend2 = ax.legend(handles, labels, loc="upper right",title="Number of genes",prop={'size': 10})
    v1 = np.linspace(c.min(), c.max(), 4, endpoint=True)
    print(v1)
    cbar = fig.colorbar(sc,label="% of genes from pathway",shrink=0.5,orientation="horizontal",ticks=v1)
    ax.set_title(title)
    fig.savefig("{}GO_{}_with_fdr.pdf".format(save_path,save_title),bbox_inches="tight")

# %%

def plot_(GO,save_title):
    df = pd.read_csv(GO,sep="\t")
    #df = df.iloc[[0,4,6,7,9,10,12,13]]
    df = df.iloc[[0,1,2,3,4,5,6,7,8]]
    #df = df[df["Ontology Type"] == "biological_process"]
    #df = df.sort_values(by=["PermuteP"],ascending=False)
    df = df.sort_values(by=["pValue"],ascending=False)
    #p_value = df["PermuteP"]
    p_value = df["pValue"]
    #name = df["Ontology Name"]
    name = df["description"]
    #n_genes = df["Number Changed"]
    n_genes = df["overlap"]*10
    fdr = df["FDR"]
    total_in_path = df["size"]
    gene_ratio = df.overlap/total_in_path*100
    #gene_ratio = df["Percent Changed"]
    #genes = df["gene symbols"]
    genes = df["userId"]
    aux_genes = []
    for g in genes:
        #l = g.split("|")[:5]
        l = g.split(";")[:5]
        aux_genes.append(", ".join(l))
    plot_fig(p_value,name,n_genes,gene_ratio,aux_genes,"GO analysis. {}. UNI = Expressed".format(save_title),save_title)
    plot_fig_fdr(p_value,name,n_genes,gene_ratio,aux_genes,"GO analysis. {}. UNI = Expressed".format(save_title),save_title,fdr)

plot_(root+"enrichment_results_EZH_minus.txt","EZH_minus")
plot_(root+"enrichment_results_EZH_plus.txt","EZH_plus")
plot_(root+"enrichment_results_ERT_minus.txt","ERT_minus")
plot_(root+"enrichment_results_ERT_plus.txt","ERT_plus")
plot_(root+"CAST_GO_webgestalt.txt","CAST_GO")
plot_(root+"S129_GO_webgestalt.txt","S129_GO")
plot_(root+"ASE_GO.csv","ASE_GO")

# %%
