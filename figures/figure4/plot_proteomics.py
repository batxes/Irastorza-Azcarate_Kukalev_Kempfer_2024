#%%

# plots for the analysis of mass spectrometry data

import re
from adjustText import adjust_text
import pandas as pd
import pybedtools as pb
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
import numpy as np
import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *
import plotly.graph_objects as go
from scipy.stats import zscore 
prot = pd.read_csv("{}/proteinGroups_LondonSamples.txt".format(root),sep="\t")
prot.rename(columns={"Gene names":"gene_names","Ratio H/L ert-":"ert_minus","Ratio H/L ert+":"ert_plus","Ratio H/L ezh-":"ezh_minus","Ratio H/L ezh+":"ezh_plus","Ratio H/L OS25-":"os25_minus","Ratio H/L OS25+":"os25_plus"},inplace=True)
prot = prot[~prot["Protein IDs"].str.contains("CON__")]
prot = prot[~prot["Protein IDs"].str.contains("REV__")]

uniprot = pd.read_csv("{}uniprot_gene_names.tsv".format(root),sep="\t")[["Entry","Gene Names"]] #Downloaded from Uniprot
aux_dict = dict(zip(uniprot["Entry"], uniprot["Gene Names"]))
uniprot_dict = {}
for key,value in aux_dict.items():
    if value is not np.nan:
        values = value.split(" ")
        if len(values) == 1:
            uniprot_dict[key] = value
        else:
            uniprot_dict[key] = ";".join(values)
    else:
        uniprot_dict[key] = "not_gene_name_found"

ase_path = root+"300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed"
master_df = pd.read_csv(ase_path,sep="\t")
master_df = master_df[master_df.chrom.isin(chrs_)]

def set_type(row):
    if row["log2foldchange"] >= 1 and row["p_adj_sum"] <= 0.05 and row["TPM_transcript"] >= 1:
        return "cast"
    elif row["log2foldchange"] <= -1 and row["p_adj_sum"] <= 0.05 and row["TPM_transcript"] >= 1:
        return "s129"
    elif ((row["log2foldchange"] > -1 and row["log2foldchange"] < 1) or row["p_adj_sum"] > 0.05) and row["TPM_transcript"] >= 1:
        return "exp_with_SNP"
    elif row["number_SNPs"] and row["TPM_transcript"] >= 1:
        return "exp_no_SNP"
    elif row["TPM_transcript"] < 1:
        return "not_exp"
    else:
        return "no_type"

master_df["type"] = master_df.apply(set_type,axis=1)
    
cast_names = master_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')[["gene_name"]].values.tolist()
cast_names = [x[0] for x in cast_names]
cast_names = list(set(cast_names))
s129_names = master_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')[["gene_name"]].values.tolist()
s129_names = [x[0] for x in s129_names]
s129_names = list(set(s129_names))
exp_names = master_df.query('TPM_transcript >= 1 and ((log2foldchange > -1 and log2foldchange < 1) or p_adj_sum > 0.05)')[["gene_name"]].values.tolist()
exp_names = [x[0] for x in exp_names]
exp_names = list(set(exp_names))
not_exp_names = master_df.query('TPM_transcript < 1')[["gene_name"]].values.tolist()
not_exp_names = [x[0] for x in not_exp_names]
not_exp_names = list(set(not_exp_names))
all_names = master_df[["gene_name"]].values.tolist()
all_names = [x[0] for x in all_names]
all_names = list(set(all_names))

print ("{} = {}".format(len(cast_names)+len(s129_names)+len(exp_names)+len(not_exp_names),len(all_names)))

names = master_df.query('Promoter_state == "PRCa"')[["gene_name"]].values.tolist()
names = [x[0] for x in names]
prca_names = list(set(names))
names = master_df.query('Promoter_state == "PRCr"')[["gene_name"]].values.tolist()
names = [x[0] for x in names]
prcr_names = list(set(names))
names = master_df.query('Promoter_state == "Active"')[["gene_name"]].values.tolist()
names = [x[0] for x in names]
active_names = list(set(names))
names = master_df.query('Promoter_state == "Inactive"')[["gene_name"]].values.tolist()
names = [x[0] for x in names]
inactive_names = list(set(names))
names = master_df.query('Promoter_state == "K27me3_only"')[["gene_name"]].values.tolist()
names = [x[0] for x in names]
k27me3only_names = list(set(names))
print ("{} = {}".format(len(prca_names)+len(prcr_names)+len(active_names)+len(inactive_names)+len(k27me3only_names),len(all_names)))

log2fc_dict = dict(zip(master_df["gene_name"], master_df["log2foldchange"]))

def fill_names(row):
    name = row["gene_names"]
    if name is np.nan:
        uniprot_ids = row["Protein IDs"]
        #print (uniprot_ids)
        #print (type(uniprot_ids))
        uniprot_ids = uniprot_ids.split(";")
        if len(uniprot_ids) == 1:
            id = uniprot_ids[0].split("__")
            try:
                if len(id) == 1:
                    return uniprot_dict[id[0]]
                else:
                    if id[0] != "CON":
                        return uniprot_dict[id[1]]
            except:
                return "not_gene_name_found"
        else:
            names = []
            for uniprot_id in uniprot_ids:
                id = uniprot_id.split("__")
                try:
                    if len(id) == 1:
                        names.append(uniprot_dict[id[0]])
                    else:
                        if id[0] != "CON":
                            names.append(uniprot_dict[id[1]])
                except:
                    names.append("not_gene_name_found")
            return ';'.join(names)
    else:
        return name

prot['gene_names'] = prot.apply(fill_names, axis=1)

#allelic function canonical
def allelic_func(row):
    if row["gene_names"] is np.nan:
        return "NAN",0
    else:
        names = row['gene_names'].split(";")
        to_return = []
        l2fc = []
        for name in names:
            if name in (all_names):
                if name in (cast_names):
                    l2fc.append(log2fc_dict[name])
                    to_return.append("CAST")
                elif name in (s129_names):
                    l2fc.append(log2fc_dict[name])
                    to_return.append("S129")
                elif name in (exp_names):
                    try:
                        l2fc.append(log2fc_dict[name])
                    except:
                        l2fc.append(0)
                    to_return.append("EXP")
                elif name in (not_exp_names):
                    try:
                        l2fc.append(log2fc_dict[name])
                    except:
                        l2fc.append(0)
                    to_return.append("NOT_EXP")
                else:
                    try:
                        l2fc.append(log2fc_dict[name])
                    except:
                        l2fc.append(0)
                    to_return.append("OTHER")
            else:
                try:
                    l2fc.append(log2fc_dict[name])
                except:
                    l2fc.append(0)
                to_return.append("NOT_FOUND")
        if len(to_return) == 1:
            try:
                return pd.Series([to_return[0],log2fc_dict[name]])
            except:
                return pd.Series([to_return[0],0])
        else:
            if "CAST" in to_return:
                return pd.Series(["CAST",l2fc[to_return.index("CAST")]])
            elif "S129" in to_return:
                return pd.Series(["S129",l2fc[to_return.index("S129")]])
            elif "EXP" in to_return:
                return pd.Series(["EXP",l2fc[to_return.index("EXP")]])
            elif "NOT_EXP" in to_return:
                return pd.Series(["NOT_EXP",l2fc[to_return.index("NOT_EXP")]])
            elif "OTHER" in to_return:
                return pd.Series(["OTHER",l2fc[to_return.index("OTHER")]])
            else:
                return pd.Series(["NOT_FOUND",l2fc[to_return.index("NOT_FOUND")]])
prot[['ASE',"l2fc_gene"]] = prot.apply(allelic_func, axis=1)

#allelic function canonical
def promoter_state_func(row):
    if row["gene_names"] is np.nan:
        return "NAN"
    else:
        names = row['gene_names'].split(";")
        to_return = []
        for name in names:
            if name in (all_names):
                if name in (prca_names):
                    to_return.append("PRCa")
                elif name in (prcr_names):
                    to_return.append("PRCr")
                elif name in (k27me3only_names):
                    to_return.append("K27me3_only")
                elif name in (active_names):
                    to_return.append("Active")
                elif name in (inactive_names):
                    to_return.append("Inactive")
                else:
                    to_return.append("Other")
            else:
                to_return.append("NOT_FOUND")
        if len(to_return) == 1:
            return to_return[0]
        else:
            if "PRCa" in to_return:
                return "PRCa"
            elif "PRCr" in to_return:
                return "PRCr"
            elif "K27me3_only" in to_return:
                return "K27me3_only"
            elif "Active" in to_return:
                return "Active"
            elif "Inactive" in to_return:
                return "Inactive"
            elif "Other" in to_return:
                return "Other"
            else:
                return "NOT_FOUND"
prot['Promoter_state'] = prot.apply(promoter_state_func, axis=1)

prot["L2FC_ERT"] = np.log2(prot.ert_plus/prot.ert_minus)
prot["L2FC_EZH"] = np.log2(prot.ezh_plus/prot.ezh_minus)
prot["L2FC_OS25"] = np.log2(prot.os25_plus/prot.os25_minus)

prot["L2FC_ERT_NORMALIZED"] = zscore(prot.L2FC_ERT,nan_policy="omit")
prot["L2FC_EZH_NORMALIZED"] = zscore(prot.L2FC_EZH,nan_policy="omit")
prot["L2FC_OS25_NORMALIZED"] = zscore(prot.L2FC_OS25,nan_policy="omit")

prot.to_csv("{}prot_df.tsv".format(root),sep="\t")

#do alluvial diagrams plots to see which proteins do not appear in the mass spec anymore

def pro_state_exp(row):
    minus = ""
    plus = ""
    if row[1] == row[1]:
        plus = row["Promoter_state"]
    else:
        plus = "NO_PROTEIN_FOUND"
    if row[0] == row[0]:
        minus = row["Promoter_state"]
    else:
        minus = "NO_PROTEIN_FOUND"
    return pd.Series([minus,plus])
         

ert_exp = prot.query('ert_minus == ert_minus or ert_plus == ert_plus')[["ert_minus","ert_plus","Promoter_state"]]
ert_exp[['Promoter_state_minus','Promoter_state_plus']] = ert_exp.apply(pro_state_exp, axis=1)
ezh_exp = prot.query('ezh_minus == ezh_minus or ezh_plus == ezh_plus')[["ezh_minus","ezh_plus","Promoter_state"]]
ezh_exp[['Promoter_state_minus','Promoter_state_plus']] = ezh_exp.apply(pro_state_exp, axis=1)
os25_exp = prot.query('os25_minus == os25_minus or os25_plus == os25_plus')[["os25_minus","os25_plus","Promoter_state"]]
os25_exp[['Promoter_state_minus','Promoter_state_plus']] = os25_exp.apply(pro_state_exp, axis=1)

def create_alluvial(df,label):
    class_minus = go.parcats.Dimension(
    values=df.Promoter_state_minus,
    #categoryorder='category ascending',
    label="minus experiment",
    #categoryarray=["NO_PROTEIN_FOUND", "NOT_FOUND", "Active","Other","PRCa","Inactive","K27me3_only","PRCr"]
    #categoryarray=[0,1,2,3,4,5,6,7]
    )
    class_plus = go.parcats.Dimension(
    values=df.Promoter_state_plus, label="plus experiment"
    )
    color = df.Promoter_state_minus
    colorscale = [[0, 'black'], [0.14, 'white'],[0.28,'green'],[0.42,'darkgrey'],[0.56,'pink'],[0.7,'grey'],[0.84,'red'],[1,'crimson']]
    #colorscale = [[0, 'darkblue'], [1,'crimson']]
    fig = go.Figure(data = [go.Parcats(dimensions=[class_minus,  class_plus],
        line={'color': color, 'colorscale': colorscale},),
        ])
    fig.update_layout(
        autosize=False,
        width=400,
        height=800,
        title=label,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
            )
    )

    fig.write_image("{}Promoter_state_changes_in_{}.pdf".format(save_path,label))
    fig.show()
mapping = {'NO_PROTEIN_FOUND': 0, 'NOT_FOUND': 1, 'Active':2, 'Other':3, 'PRCa':4, 'Inactive':5, 'K27me3_only':6, 'PRCr':7}
ert_exp = ert_exp.replace({'Promoter_state_minus': mapping})
ezh_exp = ezh_exp.replace({'Promoter_state_minus': mapping})
os25_exp = os25_exp.replace({'Promoter_state_minus': mapping})
create_alluvial(ert_exp[["Promoter_state_minus","Promoter_state_plus"]],"ert")
create_alluvial(ezh_exp[["Promoter_state_minus","Promoter_state_plus"]],"ezh")
create_alluvial(os25_exp[["Promoter_state_minus","Promoter_state_plus"]],"os25")


import matplotlib.pyplot as plt

fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(12,12),sharex=True)
ax1.hist(prot.ert_minus,bins=100)
ax2.hist(prot.ert_plus,bins=300)
ax3.hist(prot.ezh_minus,bins=200)
ax4.hist(prot.ezh_plus,bins=100)
ax5.hist(prot.os25_minus,bins=100)
ax6.hist(prot.os25_plus,bins=100)
ax1.set_title("ERT -")
ax2.set_title("ERT +")
ax3.set_title("EZH -")
ax4.set_title("EZH +")
ax5.set_title("OS25 -")
ax6.set_title("OS25 +")
ax5.set_xlabel("RATIO H/L ")
ax6.set_xlabel("RATIO H/L ")
ax1.set_xlim([0,3])


fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(12,12),sharex=True)
ax1.hist(prot.L2FC_ERT,bins=200)
ert_per_95 = np.percentile(prot.L2FC_ERT.dropna(),95)
ert_per_5 = np.percentile(prot.L2FC_ERT.dropna(),5)
ax1.axvline(ert_per_5,color="black")
ax1.axvline(ert_per_95,color="black")
ax1.axvline(0,color="black")
print (ert_per_5,ert_per_95)
ax2.hist(prot.L2FC_EZH,bins=200)
ezh_per_95 = np.percentile(prot.L2FC_EZH.dropna(),95)
ezh_per_5 = np.percentile(prot.L2FC_EZH.dropna(),5)
ax2.axvline(ezh_per_5,color="black")
ax2.axvline(ezh_per_95,color="black")
ax2.axvline(0,color="black")
ax3.hist(prot.L2FC_OS25,bins=200)
os25_per_95 = np.percentile(prot.L2FC_OS25.dropna(),95)
os25_per_5 = np.percentile(prot.L2FC_OS25.dropna(),5)
ax3.axvline(os25_per_5,color="black")
ax3.axvline(os25_per_95,color="black")
ax3.axvline(0,color="black")
ax1.set_title("ERT")
ax2.set_title("EZH")
ax3.set_title("OS25")
ax3.set_xlabel("L2FC")
ax1.set_xlim([-3,3])
ax3.set_xlim([-3,3])
ax2.set_xlim([-3,3])

fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(12,12),sharex=True)
ax1.hist(prot.L2FC_ERT_NORMALIZED,bins=400)
ert_per_95 = np.percentile(prot.L2FC_ERT_NORMALIZED.dropna(),95)
ert_per_5 = np.percentile(prot.L2FC_ERT_NORMALIZED.dropna(),5)
ax1.axvline(ert_per_5,color="black")
ax1.axvline(ert_per_95,color="black")
ax1.axvline(0,color="black")
print (ert_per_5,ert_per_95)
ax2.hist(prot.L2FC_EZH_NORMALIZED,bins=400)
ezh_per_95 = np.percentile(prot.L2FC_EZH_NORMALIZED.dropna(),95)
ezh_per_5 = np.percentile(prot.L2FC_EZH_NORMALIZED.dropna(),5)
ax2.axvline(ezh_per_5,color="black")
ax2.axvline(ezh_per_95,color="black")
ax2.axvline(0,color="black")
ax3.hist(prot.L2FC_OS25_NORMALIZED,bins=400)
os25_per_95 = np.percentile(prot.L2FC_OS25_NORMALIZED.dropna(),95)
os25_per_5 = np.percentile(prot.L2FC_OS25_NORMALIZED.dropna(),5)
ax3.axvline(os25_per_5,color="black")
ax3.axvline(os25_per_95,color="black")
ax3.axvline(0,color="black")
ax1.set_title("ERT")
ax2.set_title("EZH")
ax3.set_title("OS25")
ax3.set_xlabel("L2FC_NORMALIZED")
ax1.set_xlim([-5,5])
ax3.set_xlim([-5,5])
ax2.set_xlim([-5,5])

fig,ax = plt.subplots(figsize=(24,12))
ax.scatter(prot.L2FC_EZH,np.log10(prot["Intensity ezh+"]),color="grey")
ax.axvline(ezh_per_5,color="black")
ax.axvline(ezh_per_95,color="black")
ax.axvline(0,color="black")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Hist')]
ax.scatter(aux.L2FC_EZH,np.log10(aux["Intensity ezh+"]),color="red")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Rps.*|Rpl.* ')]
ax.scatter(aux.L2FC_EZH,np.log10(aux["Intensity ezh+"]),color="orange")
ax.set_xlabel("L2FC EZH")
ax.set_ylabel("Log10 Intensity ezh+")

fig,ax = plt.subplots(figsize=(24,12))
ax.scatter(prot.L2FC_ERT,np.log10(prot["Intensity ert+"]),color="grey")
ax.axvline(ert_per_5,color="black")
ax.axvline(ert_per_95,color="black")
ax.axvline(0,color="black")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Hist')]
ax.scatter(aux.L2FC_ERT,np.log10(aux["Intensity ert+"]),color="red")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Rps.*|Rpl.* ')]
ax.scatter(aux.L2FC_ERT,np.log10(aux["Intensity ert+"]),color="orange")
ax.set_xlabel("L2FC ERT")
ax.set_ylabel("Log10 Intensity ert+")

fig,ax = plt.subplots(figsize=(24,12))
ax.scatter(prot.L2FC_OS25,np.log10(prot["Intensity OS25+"]),color="grey")
ax.axvline(os25_per_5,color="black")
ax.axvline(os25_per_95,color="black")
ax.axvline(0,color="black")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Hist')]
ax.scatter(aux.L2FC_OS25,np.log10(aux["Intensity OS25+"]),color="red")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Rps.*|Rpl.* ')]
ax.scatter(aux.L2FC_OS25,np.log10(aux["Intensity OS25+"]),color="orange")
ax.set_xlabel("L2FC OS25")
ax.set_ylabel("Log10 Intensity os25+")
#ax.set_ylim([0,0.01])

fig,ax = plt.subplots(figsize=(24,12))
ax.scatter(prot.L2FC_EZH_NORMALIZED,np.log10(prot["Intensity ezh+"]),color="grey")
ax.axvline(ezh_per_5,color="black")
ax.axvline(ezh_per_95,color="black")
ax.axvline(0,color="black")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Hist')]
ax.scatter(aux.L2FC_EZH_NORMALIZED,np.log10(aux["Intensity ezh+"]),color="red")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Rps.*|Rpl.* ')]
ax.scatter(aux.L2FC_EZH_NORMALIZED,np.log10(aux["Intensity ezh+"]),color="orange")
ax.set_xlabel("L2FC EZH_NORMALIZED")
ax.set_ylabel("Log10 Intensity ezh+")

fig,ax = plt.subplots(figsize=(24,12))
ax.scatter(prot.L2FC_ERT_NORMALIZED,np.log10(prot["Intensity ert+"]),color="grey")
ax.axvline(ert_per_5,color="black")
ax.axvline(ert_per_95,color="black")
ax.axvline(0,color="black")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Hist')]
ax.scatter(aux.L2FC_ERT_NORMALIZED,np.log10(aux["Intensity ert+"]),color="red")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Rps.*|Rpl.* ')]
ax.scatter(aux.L2FC_ERT_NORMALIZED,np.log10(aux["Intensity ert+"]),color="orange")
ax.set_xlabel("L2FC ERT_NORMALIZED")
ax.set_ylabel("Log10 Intensity ert+")

fig,ax = plt.subplots(figsize=(24,12))
ax.scatter(prot.L2FC_OS25_NORMALIZED,np.log10(prot["Intensity OS25+"]),color="grey")
ax.axvline(os25_per_5,color="black")
ax.axvline(os25_per_95,color="black")
ax.axvline(0,color="black")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Hist')]
ax.scatter(aux.L2FC_OS25_NORMALIZED,np.log10(aux["Intensity OS25+"]),color="red")
aux = prot.dropna(subset=["gene_names"])
aux = aux[aux["gene_names"].str.contains('^Rps.*|Rpl.* ')]
ax.scatter(aux.L2FC_OS25_NORMALIZED,np.log10(aux["Intensity OS25+"]),color="orange")
ax.set_xlabel("L2FC OS25_NORMALIZED")
ax.set_ylabel("Log10 Intensity os25+")
#ax.set_ylim([0,0.01])

def plot_tree(x,y,cutoff1,cutoff2):
    fig,ax = plt.subplots(figsize=(24,12))
    ax.scatter(prot[x],np.log10(prot[y]),color="lightgrey")
    #ax.axvline(cutoff1,color="black")
    #ax.axvline(cutoff2,color="black")
    ax.axvline(0,color="black")
    aux = prot.query("ASE == 'CAST'")
    ax.scatter(aux[x],np.log10(aux[y]),color=cast_color)
    aux = prot.query("ASE == 'S129'")
    ax.scatter(aux[x],np.log10(aux[y]),color=s129_color)
    ax.set_xlabel(x)
    ax.set_ylabel("Log 10 {}".format(y))

def plot_tree_promoter_states(x,y,cutoff1,cutoff2):
    fig,ax = plt.subplots(figsize=(24,12))
    ax.scatter(prot[x],np.log10(prot[y]),color="lightgrey")
    #ax.axvline(cutoff1,color="black")
    #ax.axvline(cutoff2,color="black")
    ax.axvline(0,color="black")
    #aux = prot.query("Promoter_state == 'Active'")
    #ax.scatter(aux[x],np.log10(aux[y]),color="palegreen")
    aux = prot.query("Promoter_state == 'PRCa'")
    ax.scatter(aux[x],np.log10(aux[y]),color="crimson")
    aux = prot.query("Promoter_state == 'PRCr'")
    ax.scatter(aux[x],np.log10(aux[y]),color="blue")
    aux = prot.query("Promoter_state == 'K27me3_only'")
    ax.scatter(aux[x],np.log10(aux[y]),color="lime")
    ax.set_xlabel(x)
    ax.set_ylabel("Log 10 {}".format(y))

cast_color = "#611163"
s129_color = "#EA8D1F"
both_color = "#7A7A7A"
non_phased_color = "#09AA9E"
plot_tree("L2FC_ERT","Intensity ert+",ert_per_5,ert_per_95)
plot_tree("L2FC_EZH","Intensity ezh+",ezh_per_5,ezh_per_95)
plot_tree("L2FC_OS25","Intensity OS25+",os25_per_5,os25_per_95)
plot_tree("L2FC_ERT_NORMALIZED","Intensity ert+",ert_per_5,ert_per_95)
plot_tree("L2FC_EZH_NORMALIZED","Intensity ezh+",ezh_per_5,ezh_per_95)
plot_tree("L2FC_OS25_NORMALIZED","Intensity OS25+",os25_per_5,os25_per_95)
plot_tree_promoter_states("L2FC_ERT_NORMALIZED","Intensity ert+",ert_per_5,ert_per_95)
plot_tree_promoter_states("L2FC_EZH_NORMALIZED","Intensity ezh+",ezh_per_5,ezh_per_95)
plot_tree_promoter_states("L2FC_OS25_NORMALIZED","Intensity OS25+",os25_per_5,os25_per_95)



#%%

# tree and boxplots shown in figures 4d, 4f and SI figures 6e, 6f

def plot_pretty_tree(x,y,z,type_):
    vmin = 5
    vmax = 10
    fig,ax = plt.subplots(figsize=(18,10))
    pprot = prot[["gene_names",x,y,z]]
    pprot = pprot.dropna()
    abundance_index = np.log10(pprot[y]/pprot[z])
    intensities = np.log10(pprot[y])    
    colormap = plt.cm.get_cmap('viridis')
    sc = ax.scatter(pprot[x],abundance_index,c = intensities,edgecolors="black",alpha=0.7,s=100,cmap=colormap,vmin=vmin, vmax=vmax)
    ax.axvline(0,color="black")
    ax.set_xlabel(x)
    ax.set_ylabel("log10({}/{})".format(y,z))
    cax = fig.add_axes([0.91, 0.25, 0.025, 0.5])
    cbar = fig.colorbar(sc, cax=cax, orientation='vertical')
    cbar.ax.set_title('log10({})'.format(y))
    annotate = True
    pprot["abundance_index"] = abundance_index
    pprot["intensity"] = intensities
    annotate_df = pprot.copy()
    
    annotate_df1 = annotate_df[annotate_df["gene_names"].str.contains('^Hist') ]
    annotate_df2 = annotate_df[annotate_df["gene_names"].str.contains('^Rps.*|Rpl.* ')]
    annotate_df = pd.concat([annotate_df1,annotate_df2])
    x_label = annotate_df[x].values
    y_label = np.log10(annotate_df[y]/annotate_df[z])
    z_label = annotate_df.gene_names.values
    if annotate:
        for (a,b,c) in zip(x_label,y_label,z_label):
            print (c)
            ax.annotate(c, xy=(a-0.25,b+0.1))
    ax.set_xlim([-5,5])
    fig.savefig("{}annotated_tree_{}_plot.pdf".format(save_path,type_),dpi=300)
    fig.savefig("{}annotated_tree_{}_plot.png".format(save_path,type_),dpi=300)

    # additional plot for hist and ribo proteins
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(6,6))
    aux = pprot.dropna(subset=["gene_names"])
    aux = aux[aux["gene_names"].str.contains('^Hist')]
    aux2 = pprot.dropna(subset=["gene_names"])
    aux2 = aux2[aux2["gene_names"].str.contains('^Rps.*|Rpl.* ')]
    aa = [aux["abundance_index"],aux2["abundance_index"]]
    ll = [aux[x],aux2[x]]
    bp1 = ax1.boxplot(aa,showfliers=False,patch_artist=True,widths=0.7)

    for patch in bp1['boxes']:
        patch.set_facecolor("none")
        patch.set(alpha=0.9)
    plt.setp(bp1["medians"],color="black")
    bp2 = ax2.boxplot(ll,showfliers=False,patch_artist=True,widths=0.7)
    for patch in bp2['boxes']:
        patch.set_facecolor("none")
        patch.set(alpha=0.9)
    plt.setp(bp2["medians"],color="black")
    for pos,v in enumerate(aa): 
        xv = np.random.normal(pos+1, 0.05, len(v))
        ax1.scatter(xv, v,s=6,alpha=1,color="blue")
    for pos,v in enumerate(ll): 
        xv = np.random.normal(pos+1, 0.05, len(v))
        ax2.scatter(xv, v,s=6,alpha=1,color="crimson")
    ax1.set_ylabel("log10({}/{})".format(y,z))
    ax1.set_xticklabels(["Hist","Ribo"])
    ax2.set_ylabel("{}".format(x))
    ax2.set_xticklabels(["Hist","Ribo"])
    ax2.set_ylim([-4.5,4.5])
    fig.savefig("{}box_hist_ribo_{}_plot.pdf".format(save_path,type_),dpi=300)
    fig.savefig("{}box_hist_ribo_{}_plot.png".format(save_path,type_),dpi=300)

def plot_pretty_tree_promotor_state(x,y,z):
    vmin = 5
    vmax = 10
    fig,ax = plt.subplots(figsize=(18,10))
    pprot = prot[["gene_names",x,y,z,"Promoter_state"]]
    pprot = pprot.dropna()
    aux = prot.query("Promoter_state == 'PRCa'")
    abundance_index = np.log10(aux[y]/aux[z])
    ax.scatter(aux[x],np.log10(aux[y]/aux[z]),color="crimson",vmin=vmin, vmax=vmax,alpha=0.7, s=100)
    aux = prot.query("Promoter_state == 'PRCr'")
    ax.scatter(aux[x],np.log10(aux[y]/aux[z]),color="blue",vmin=vmin, vmax=vmax,alpha=0.7, s=100)
    aux = prot.query("Promoter_state == 'K27me3_only'")
    ax.scatter(aux[x],np.log10(aux[y]/aux[z]),color="darkgreen",vmin=vmin, vmax=vmax,alpha=0.7, s=100)

    ax.axvline(0,color="black")
    ax.set_xlabel(x)
    ax.set_ylabel("Log 10 {}".format(y))
    cax = fig.add_axes([0.91, 0.25, 0.025, 0.5])

    annotate = True
    pprot["abundance_index"] = abundance_index
    annotate_df = pprot.query('abundance_index > 6.8 or (abundance_index >= 5 and L2FC_EZH_NORMALIZED < -1)')
    
    x = annotate_df[x].values
    y = np.log10(annotate_df[y]/annotate_df[z])
    z = annotate_df.gene_names.values
    if annotate:
        for (a,b,c) in zip(x,y,z):
            print (c)
            ax.annotate(c, xy=(a-0.25,b+0.1))
    ax.set_xlim([-5,5])
    fig.savefig(save_path+"annotated_tree_plot_promoter_state.pdf",dpi=300)


def plot_pretty_tree_ase(x,y,z):
    vmin = 5
    vmax = 10
    fig,ax = plt.subplots(figsize=(18,10))
    pprot = prot[["gene_names",x,y,z,"ASE"]]
    pprot = pprot.dropna()
    aux = prot.query("ASE == 'CAST' or ASE == 'S129'")
    abundance_index = np.log10(aux[y]/aux[z])
    aux = prot.query("ASE == 'CAST'")
    ax.scatter(aux[x],np.log10(aux[y]/aux[z]),color=cast_color,vmin=vmin, vmax=vmax,alpha=0.7, s=100)
    aux = prot.query("ASE == 'S129'")
    ax.scatter(aux[x],np.log10(aux[y]/aux[z]),color=s129_color,vmin=vmin, vmax=vmax,alpha=0.7, s=100)
    ax.axvline(0,color="black")
    ax.set_xlabel(x)
    ax.set_ylabel("Log 10 {}".format(y))
    cax = fig.add_axes([0.91, 0.25, 0.025, 0.5])
    annotate = True
    pprot["abundance_index"] = abundance_index
    annotate_df = pprot.query('abundance_index > 6.8 or (abundance_index >= 5 and L2FC_EZH_NORMALIZED < -1)')
    x = annotate_df[x].values
    y = np.log10(annotate_df[y]/annotate_df[z])
    z = annotate_df.gene_names.values
    if annotate:
        for (a,b,c) in zip(x,y,z):
            print (c)
            ax.annotate(c, xy=(a-0.25,b+0.1))
    ax.set_xlim([-5,5])
    fig.savefig(save_path+"annotated_tree_plot_ASE.pdf",dpi=300)

def plot_pretty_tree_ase_PRC(x,y,z):
    vmin = 5
    vmax = 10
    fig,ax = plt.subplots(figsize=(18,10))
    pprot = prot[["gene_names",x,y,z,"ASE","Promoter_state"]]
    pprot = pprot.dropna()
    aux = prot.query("(Promoter_state == 'PRCa' or Promoter_state == 'PRCr' or Promoter_state == 'K27me3_only' ) and (ASE == 'CAST' or ASE == 'S129')")
    abundance_index = np.log10(aux[y]/aux[z])
    aux = prot.query("(Promoter_state == 'PRCa' or Promoter_state == 'PRCr' or Promoter_state == 'K27me3_only' ) and (ASE == 'CAST')")
    ax.scatter(aux[x],np.log10(aux[y]/aux[z]),color=cast_color,vmin=vmin, vmax=vmax,alpha=0.7, s=100)
    aux = prot.query("(Promoter_state == 'PRCa' or Promoter_state == 'PRCr' or Promoter_state == 'K27me3_only' ) and (ASE == 'S129')")
    ax.scatter(aux[x],np.log10(aux[y]/aux[z]),color=s129_color,vmin=vmin, vmax=vmax,alpha=0.7, s=100)
    ax.axvline(0,color="black")
    ax.set_xlabel(x)
    ax.set_ylabel("Log 10 {}".format(y))
    cax = fig.add_axes([0.91, 0.25, 0.025, 0.5])
    annotate = True
    pprot["abundance_index"] = abundance_index
    annotate_df = pprot.query('abundance_index > 6.8 or (abundance_index >= 5 and L2FC_EZH_NORMALIZED < -1)')
    x = annotate_df[x].values
    y = np.log10(annotate_df[y]/annotate_df[z])
    z = annotate_df.gene_names.values
    if annotate:
        for (a,b,c) in zip(x,y,z):
            print (c)
            ax.annotate(c, xy=(a-0.25,b+0.1))
    ax.set_xlim([-5,5])
    fig.savefig(save_path+"annotated_tree_plot_ASE_and_PRC.pdf",dpi=300)

plot_pretty_tree("L2FC_EZH_NORMALIZED","Intensity ezh+","Peptides ezh-","EZH")
plot_pretty_tree_promotor_state("L2FC_EZH_NORMALIZED","Intensity ezh-","Peptides ezh+")
plot_pretty_tree_ase("L2FC_EZH_NORMALIZED","Intensity ezh+","Peptides ezh-")
plot_pretty_tree_ase_PRC("L2FC_EZH_NORMALIZED","Intensity ezh+","Peptides ezh-")
plot_pretty_tree("L2FC_ERT_NORMALIZED","Intensity ert+","Peptides ert-","ERT")
plot_pretty_tree("L2FC_OS25_NORMALIZED","Intensity OS25+","Peptides OS25-","OS25")
#%%
# Hist and/or Ribo table

table = prot[prot["gene_names"].str.contains('^Hist')][["gene_names","ASE","Promoter_state",'L2FC_ERT_NORMALIZED', 'L2FC_EZH_NORMALIZED',
       'L2FC_OS25_NORMALIZED']]
table.to_csv(save_path+"mass_spec_table_updated.csv",sep="\t")


#%%


cast_ =  (len(prot.query('ASE == "CAST"')))
s129_ = (len(prot.query('ASE == "S129"')))
exp_ = (len(prot.query('ASE == "EXP"')))
# (len(prot.query('ASE == "NOT_ASE"')))
# (len(prot.query('ASE == "NOT_FOUND"')))
all_ = (len(prot))
prca_ = (len(prot.query('Promoter_state == "PRCa"')))
prcr_ =  (len(prot.query('Promoter_state == "PRCr"')))
k27meonly_ = (len(prot.query('Promoter_state == "K27me3_only"')))

fig,ax = plt.subplots(figsize=(8,15))
ax.bar([0,1,2,3,4,5,6],[len(cast_names),len(s129_names),len(exp_names),len(prca_names),len(prcr_names),len(k27me3only_names),len(all_names)],color="grey")
ax.bar([0,1,2,3,4,5,6],[cast_,s129_,exp_,prca_,prcr_,k27meonly_, all_],color=[cast_color,s129_color,"darkgreen","crimson","blue","lime","black"])
ax.set_ylabel("Number of proteins")
ax.set_xticklabels(["","CAST","S129","Expressed","PRCa","PRCr","K27me3_only","All"])

ax.text(-0.2,30500,"{:.2f}".format(100*cast_/len(cast_names)))
ax.text(0.8,30500,"{:.2f}".format(100*s129_/len(s129_names)))
ax.text(1.8,30500,"{:.2f}".format(100*exp_/len(exp_names)))
ax.text(2.8,30500,"{:.2f}".format(100*prca_/len(prca_names)))
ax.text(3.8,30500,"{:.2f}".format(100*prcr_/len(prcr_names)))
ax.text(4.8,30500,"{:.2f}".format(100*k27meonly_/len(k27me3only_names)))
ax.text(5.8,30500,"{:.2f}".format(100*all_/len(all_names)))

#scatter plot for log2foldchange of ASE and l2fc of peptides
def plot_vs_l2fc(exp):
    fig,ax = plt.subplots(figsize=(12,12))
    aux = prot.query('l2fc_gene >= 1')
    ax.scatter(aux[exp],aux.l2fc_gene,s=10,color=cast_color)
    aux = prot.query('l2fc_gene <= -1')
    ax.scatter(aux[exp],aux.l2fc_gene,s=10,color=s129_color)
    aux = prot.query('l2fc_gene > -1 and l2fc_gene < 1')
    ax.scatter(aux[exp],aux.l2fc_gene,s=10,color="grey")
    ax.axvline(0,color="black")
    ax.axhline(1,color="black")
    ax.axhline(-1,color="black")
    ax.set_xlabel(exp)
    ax.set_ylabel("L2FC RNA-seq")
    ax.grid()

    fig,ax = plt.subplots(figsize=(12,12))
    aux = prot.query('l2fc_gene >= 1')
    ax.scatter(aux[exp],aux.l2fc_gene,s=10,color=cast_color)
    aux = prot.query('l2fc_gene <= -1')
    ax.scatter(aux[exp],aux.l2fc_gene,s=10,color=s129_color)
    aux = prot.query('l2fc_gene > -1 and l2fc_gene < 1')
    ax.scatter(aux[exp],aux.l2fc_gene,s=10,color="grey")
    ax.axvline(0,color="black")
    ax.axhline(1,color="black")
    ax.axhline(-1,color="black")
    ax.set_xlim([-5,5])
    ax.grid()
    ax.set_xlabel(exp)
    ax.set_ylabel("L2FC RNA-seq")

plot_vs_l2fc("L2FC_ERT_NORMALIZED")
plot_vs_l2fc("L2FC_EZH_NORMALIZED")
plot_vs_l2fc("L2FC_OS25_NORMALIZED")

aux = prot[["gene_names","L2FC_EZH_NORMALIZED","L2FC_ERT_NORMALIZED","L2FC_OS25_NORMALIZED","Intensity ezh+","Sequence length"]]
dict_ezh = {}
dict_ert = {}
dict_os25 = {}
dict_ezh_intensity = {}
dict_seq_length = {}
for i, row in aux.iterrows():
    gene_names = row["gene_names"]
    l2fc_ezh = row["L2FC_EZH_NORMALIZED"]
    l2fc_ert = row["L2FC_ERT_NORMALIZED"]
    l2fc_os25 = row["L2FC_OS25_NORMALIZED"]
    l2fc_ezh_int = row["Intensity ezh+"]
    l2fc_seq_len = row["Sequence length"]
    gene_names = gene_names.split(";")
    if len(gene_names) == 1:
        dict_ezh[gene_names[0]] = l2fc_ezh
        dict_ert[gene_names[0]] = l2fc_ert
        dict_os25[gene_names[0]] = l2fc_os25
        dict_ezh_intensity[gene_names[0]] = l2fc_ezh_int
        dict_seq_length[gene_names[0]] = l2fc_seq_len
    else:
        for x in gene_names:
            dict_ezh[x] = l2fc_ezh
            dict_ert[x] = l2fc_ert
            dict_os25[x] = l2fc_os25
            dict_ezh_intensity[x] = l2fc_ezh_int
            dict_seq_length[x] = l2fc_seq_len

def add_ratio_ezh(row):
    try:
        value = dict_ezh[row["gene_name"]]
    except:
        value = np.nan
    return value
def add_ratio_ert(row):
    try:
        value = dict_ert[row["gene_name"]]
    except:
        value = np.nan
    return value
def add_ratio_os25(row):
    try:
        value = dict_os25[row["gene_name"]]
    except:
        value = np.nan
    return value
def add_intensity(row):
    try:
        value = dict_ezh_intensity[row["gene_name"]]
    except:
        value = np.nan
    return value
def add_seq_len(row):
    try:
        value = dict_seq_length[row["gene_name"]]
    except:
        value = np.nan
    return value

master_df["L2FC_EZH_NORMALIZED"] = master_df.apply(add_ratio_ezh,axis=1)
master_df["L2FC_ERT_NORMALIZED"] = master_df.apply(add_ratio_ert,axis=1)
master_df["L2FC_OS25_NORMALIZED"] = master_df.apply(add_ratio_os25,axis=1)
master_df["Intensity ezh+"] = master_df.apply(add_intensity,axis=1)
master_df["Sequence length"] = master_df.apply(add_seq_len,axis=1)

h3k27me3_path = root+"F123.K27me3.BCP_peaks.HM_mode.bed"
h3k27me3 = pd.read_csv(h3k27me3_path,sep="\t",names=["chrom","start","end","size","score"])
s5p_path = root+"F123.S5p.BCP_peaks.HM_mode.bed"
s5p = pd.read_csv(s5p_path,sep="\t",names=["chrom","start","end","size","score"])
s5p = s5p[["chrom","start","end"]]
cast_genes = master_df.query('TPM_transcript >= 1 and log2foldchange >= 1 and p_adj_sum <= 0.05')[["chrom","start_transcript","end_transcript"]]
cast_genes = cast_genes.rename(columns={"start_transcript":"start","end_transcript":"end"})
s129_genes = master_df.query('TPM_transcript >= 1 and log2foldchange <= -1 and p_adj_sum <= 0.05')[["chrom","start_transcript","end_transcript"]]
s129_genes = s129_genes.rename(columns={"start_transcript":"start","end_transcript":"end"})
exp_genes = master_df.query('TPM_transcript >= 1 and ((log2foldchange > -1 and log2foldchange < 1) or p_adj_sum > 0.05)')[["chrom","start_transcript","end_transcript"]]
exp_genes = exp_genes.rename(columns={"start_transcript":"start","end_transcript":"end"})
prca_genes = master_df.query('Promoter_state == "PRCa"')[["chrom","start_transcript","end_transcript"]]
prca_genes = prca_genes.rename(columns={"start_transcript":"start","end_transcript":"end"})

lads_path = root+"mESC.mm10.liftover.LADs.bed"
lads_df = pd.read_csv(lads_path,sep="\t")[["chrom","start","end"]]
#%%
table = master_df[master_df["gene_name"].str.contains('^Hist')][['chrom','start_transcript',
       'end_transcript',"gene_name","type","Promoter_state",'L2FC_ERT_NORMALIZED', 'L2FC_EZH_NORMALIZED',
       'L2FC_OS25_NORMALIZED']]

table.to_csv(save_path+"mass_spec_table_updated_2.csv",sep="\t")

#%%


def plot_pretty_tree_close_to_promoter(x,y,z):
    master_pb = pb.BedTool.from_dataframe(master_df).sort()
    fig,ax = plt.subplots(figsize=(18,10))
    h3k27me3_df = h3k27me3[h3k27me3["chrom"].isin(chrs_)][["chrom","start","end"]]
    h3k27me3_pb = pb.BedTool.from_dataframe(h3k27me3_df).sort()
    closest = master_pb.closest(h3k27me3_pb,d=True)
    pprot = closest.to_dataframe(names=master_df.columns.values.tolist()+["prc_chrom","prc_start","prc_end","distance"])
    pprot = pprot.replace(".",np.nan)
    pprot[y] = pprot[y].astype(float)
    pprot[z] = pprot[z].astype(float)
    pprot["L2FC_EZH_NORMALIZED"] = pprot["L2FC_EZH_NORMALIZED"].astype(float)

    #pprot = pprot.dropna()
    print (pprot[y])
    print (pprot[z])
    abundance_index = np.log10((pprot[y])/(pprot[z]))

    colormap = plt.cm.get_cmap('viridis_r')
    vmin = 0
    vmax = 100000
    sc = ax.scatter(pprot[x],abundance_index,c = pprot["distance"],edgecolors="black",alpha=0.7,s=50,cmap=colormap,vmin=vmin, vmax=vmax)
    #ax.scatter(pprot[x],np.log10(pprot[y]),c = sizes,s=sizes,edgecolors="black")
    ax.axvline(0,color="black")
    ax.set_xlabel(x)
    ax.set_ylabel("Log 10 {}".format(y))
    #sm = plt.cm.ScalarMappable(cmap=colormap)
    #sm.set_clim(vmin=vmin, vmax=vmax)
    cax = fig.add_axes([0.91, 0.25, 0.025, 0.5])
    #plt.colorbar(sm)
    #fig.colorbar()
    fig.colorbar(sc, cax=cax, orientation='vertical')

    annotate = False
    pprot["abundance_index"] = abundance_index
    pprot["abundance_index"] = pprot["abundance_index"].astype(float)
    annotate_df = pprot.query('abundance_index > 6.8 or (abundance_index >= 5 and L2FC_EZH_NORMALIZED < -1)')
    #annotate_df = pprot[pprot["gene_names"].isin(plot_test_names)]
    #master_df = master_df[master_df.chrom.isin(chrs_)]
    
    x = annotate_df[x].values
    y = np.log10(annotate_df[y]/annotate_df[z])
    #y = annotate_df[y].values
    z = annotate_df.gene_name.values
    if annotate:
        for (a,b,c) in zip(x,y,z):
            print (c)
            ax.annotate(c, xy=(a-0.25,b+0.1))
    ax.set_xlim([-5,5])
    fig.savefig(save_path+"annotated_tree_plot_close_to_promoters.pdf",dpi=300)

plot_pretty_tree_close_to_promoter("L2FC_EZH_NORMALIZED","Intensity ezh+","Sequence length")
