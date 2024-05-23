#%%

#script that reads insulation score in boundaries and compares between datasets

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
from scipy.stats.stats import pearsonr,spearmanr
import seaborn as sns
import datetime
date_ = datetime.datetime.now()
date = "{}-{}-{}".format(date_.day,date_.month,date_.year)

paths = [
        root+"F123_all_as_3NPs_CAST_50kb_ins400K_400000.insulation.bed",
        root+"F123_all_as_3NPs_S129_50kb_ins400K_400000.insulation.bed",
        root+"F123_all_as_3NPs_50kb_ins400K_400000.insulation.bed"
        ]
def get_chr(x):
    chr_ = x.split(":")[0]
    return chr_

dfs = []
for dataset in paths:
    df = pd.read_csv(dataset,sep="\t",names=["region","bin1","bin2","midpoint","n1","n2","n_middle","score","other","other2"])
    df = df.fillna(0)
    df["chrom"] = df["region"].apply(lambda x : get_chr(x))
    dfs.append(df)

def compare_two_datasets (df1,df2,chr_):
    aux_df1 = df1[df1.chrom == chr_].copy()
    aux_df2 = df2[df2.chrom == chr_].copy()
    values1 = aux_df1["score"].values.tolist()
    values2 = aux_df2["score"].values.tolist()
    return pearsonr(values1,values2)[0]


# %%

# scatter plots between alleles and bulk

fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(18,12)) 
aux = dfs[2]
nonphased = aux[aux.chrom.isin(chrs_)].copy()
ax1.scatter(dfs[0]["score"].values,dfs[1]["score"].values)
ax2.scatter(dfs[0]["score"].values,nonphased["score"].values)
ax3.scatter(dfs[1]["score"].values,nonphased["score"].values)


# %%
import seaborn as sns
cast_df = dfs[0][["region","score"]]
cast_df = cast_df.rename(columns={"score":"score_cast"})
s129_df = dfs[1][["region","score"]]
s129_df = s129_df.rename(columns={"score":"score_s129"})
np_df = nonphased[["region","score"]]
np_df = np_df.rename(columns={"score":"score_np"})
is_df = pd.merge(cast_df,s129_df,how="outer")
is_df = pd.merge(is_df,np_df,how="outer")
print (is_df)

is_df = is_df.query('score_cast != 0 and score_s129 != 0 and score_np != 0')
print (is_df)


#%%
fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(18,6),sharex=True,sharey=True) 
ax1.plot([-1,1],[-1,1],ls="--",color="black")
ax2.plot([-1,1],[-1,1],ls="--",color="black")
ax3.plot([-1,1],[-1,1],ls="--",color="black")
ax1.set_title("r = {}\n r_s = {}".format(pearsonr(is_df.score_cast.values,is_df.score_s129.values)[0],spearmanr(is_df.score_cast.values,is_df.score_s129.values)[0]))
ax2.set_title("r = {}\n r_s = {}".format(pearsonr(is_df.score_cast.values,is_df.score_np.values)[0],spearmanr(is_df.score_cast.values,is_df.score_np.values)[0]))
ax3.set_title("r = {}\n r_s = {}".format(pearsonr(is_df.score_s129.values,is_df.score_np.values)[0],spearmanr(is_df.score_s129.values,is_df.score_np.values)[0]))
sns.kdeplot(
    data=is_df, x="score_cast", y="score_s129",
    fill=True, thresh=0, levels=100, cmap="cubehelix_r",ax=ax1
)
sns.kdeplot(
    data=is_df, x="score_cast", y="score_np",
    fill=True, thresh=0, levels=100, cmap="cubehelix_r",ax=ax2
)
sns.kdeplot(
    data=is_df, x="score_s129", y="score_np",
    fill=True, thresh=0, levels=100, cmap="cubehelix_r",ax=ax3
)
ax1.set_xlim(-0.5,0.5)
ax1.set_ylim(-0.5,0.5)
plt.savefig("{}IS_scores_correlation_{}.pdf".format(save_path,date),dpi=300)


# %%
import scipy.stats as stats


diff = [x-y for (x,y) in zip(is_df.score_cast.values,is_df.score_s129.values)]
z_scores = stats.zscore(diff)


fig,ax=plt.subplots(figsize=(8,8))
ax.hist(diff,bins=100,color="grey")
#0.05 pvalue = 1.64 z-score
#z = (X-mean)/stdev
#X = z*stedv+mean
line = 1.64*np.std(diff)+np.mean(diff)
ax.axvline(line)
ax.axvline(line*-1)

outliers = []
for d in diff:
    if d <= line*-1 or d >= line:
        outliers.append(d)
print ("Number of outliers = {} out of {}. {}%".format(len(outliers),len(diff),len(outliers)/len(diff)))
print ("In regions = {}bp".format(len(outliers)*50000))
plt.savefig("{}IS_scores_distribution_{}.pdf".format(save_path,date),dpi=300)



# %%
