#%%

#script that reads insulation score in boundaries and compares between datasets
import pandas as pd
from scipy.stats.stats import pearsonr,spearmanr
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import datetime
date_ = datetime.datetime.now()
date = "{}-{}-{}".format(date_.day,date_.month,date_.year)

root = "../data/insulation_s../data/insulation_scores"
chrs_ = ['chr1','chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18','chr19', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']

# value 8 is the insulation score

paths = ["{}F123_R2_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_all_as_3NP_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_1NP_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_3NP_50kb_400000.insulation.bed".format(root),
        "{}F123_panhisto_CAST_50kb_400000.insulation.bed".format(root),
        "{}F123_panhisto_S129_50kb_400000.insulation.bed".format(root),
        "{}F123_R2_3NP_CAST_50kb_400000.insulation.bed".format(root),
        "{}F123_R2_3NP_S129_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_1NP_CAST_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_1NP_S129_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_3NP_CAST_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_3NP_S129_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_all_as_3NP_CAST_50kb_400000.insulation.bed".format(root),
        "{}F123_R1_all_as_3NP_S129_50kb_400000.insulation.bed".format(root),
        "/home/ibai/pombo_johnny5/Sasha/Projects/F123/TADs/F123_all_as_3NPs_Phased_at50Kb/raw_files/F123_all_as_3NPs_CAST_50kb_ins800K_800000.insulation.bed",
        "/home/ibai/pombo_johnny5/Sasha/Projects/F123/TADs/F123_all_as_3NPs_Phased_at50Kb/raw_files/F123_all_as_3NPs_S129_50kb_ins800K_800000.insulation.bed",
        "/home/ibai/pombo_johnny5/Sasha/Projects/F123/TADs/F123_all_as_3NPs_NonPhased/F123_all_as_3NPs_50kb_ins800K_800000.insulation.bed"
        #"/home/ibai/pombo_johnny5/Sasha/Projects/F123/TADs/F123_all_as_3NPs_Phased_at50Kb/raw_files/F123_all_as_3NPs_CAST_50kb_ins400K_400000.insulation.bed",
        #"/home/ibai/pombo_johnny5/Sasha/Projects/F123/TADs/F123_all_as_3NPs_Phased_at50Kb/raw_files/F123_all_as_3NPs_S129_50kb_ins400K_400000.insulation.bed",
        #"/home/ibai/pombo_johnny5/Sasha/Projects/F123/TADs/F123_all_as_3NPs_NonPhased/F123_all_as_3NPs_50kb_ins400K_400000.insulation.bed"
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

output = "/home/ibai/tads_corr.txt"
with open (output,"w") as stdout:
    stdout.write("chrom\tR1_R2\t1NP_3NP\tPAN_alleles\tR2_alleles\tR1_1NP_alleles\tR1_3NP_alleles\tR1_alleles\tF123_alleles\n")
    for chr_ in chrs_:
        corr1 = compare_two_datasets(dfs[0],dfs[1],chr_)
        corr2 = compare_two_datasets(dfs[2],dfs[3],chr_)
        corr3 = compare_two_datasets(dfs[4],dfs[5],chr_)
        corr4 = compare_two_datasets(dfs[6],dfs[7],chr_)
        corr5 = compare_two_datasets(dfs[8],dfs[9],chr_)
        corr6 = compare_two_datasets(dfs[10],dfs[11],chr_)
        corr7 = compare_two_datasets(dfs[12],dfs[13],chr_)
        corr8 = compare_two_datasets(dfs[14],dfs[15],chr_)
        stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr_,corr1,corr2,corr3,corr4,corr5,corr6,corr7,corr8))

fig = plt.figure(figsize=[18,2])
ax = fig.add_subplot(111)
df = pd.read_csv(output,sep="\t",index_col="chrom")
df = df.T
ax = df.plot(kind='bar',figsize=(18,12),legend=False)
ax.set_ylim([0,1])
plt.savefig("{}Insulation_score_correlation.pdf".format(root))
plt.close()

corr_values_1 = df.iloc[0].values
corr_values_2 = df.iloc[1].values
corr_values_3 = df.iloc[2].values
corr_values_4 = df.iloc[3].values
corr_values_5 = df.iloc[4].values
corr_values_6 = df.iloc[5].values
corr_values_7 = df.iloc[6].values
corr_values_8 = df.iloc[7].values
print (np.mean(corr_values_1), np.mean(corr_values_2), np.mean(corr_values_3), np.mean(corr_values_4), np.mean(corr_values_5), np.mean(corr_values_6), np.mean(corr_values_7),np.mean(corr_values_8))

df1 = pd.DataFrame({ 'group' : np.repeat("R1_R2",len(chrs_)), 'value': corr_values_1 }) 
df2 = pd.DataFrame({ 'group' : np.repeat("1NP_3NP",len(chrs_)), 'value': corr_values_2 }) 
df3 = pd.DataFrame({ 'group' : np.repeat("PAN_alleles",len(chrs_)), 'value': corr_values_3 }) 
df4 = pd.DataFrame({ 'group' : np.repeat("R2_alleles",len(chrs_)), 'value': corr_values_4 }) 
df5 = pd.DataFrame({ 'group' : np.repeat("R1_1NP_alleles",len(chrs_)), 'value': corr_values_5 }) 
df6 = pd.DataFrame({ 'group' : np.repeat("R1_3NP_alleles",len(chrs_)), 'value': corr_values_6 }) 
df7 = pd.DataFrame({ 'group' : np.repeat("R1_alleles",len(chrs_)), 'value': corr_values_7 }) 
df8 = pd.DataFrame({ 'group' : np.repeat("F123_alleles",len(chrs_)), 'value': corr_values_8 }) 
df=df1.append(df2).append(df3).append(df4).append(df5).append(df6).append(df7).append(df8)
f, ax = plt.subplots(figsize=(18, 12))
sns.violinplot( x='group', y='value', data=df)
plt.ylim(0,1)
plt.savefig("{}IS_score_correlation_violin.pdf".format(root))
# %%

# scatter plots between alleles and bulk

fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(18,12)) 
aux = dfs[16]
nonphased = aux[aux.chrom.isin(chrs_)].copy()
ax1.scatter(dfs[14]["score"].values,dfs[15]["score"].values)
ax2.scatter(dfs[14]["score"].values,nonphased["score"].values)
ax3.scatter(dfs[15]["score"].values,nonphased["score"].values)




# %%
import seaborn as sns
cast_df = dfs[14][["region","score"]]
cast_df = cast_df.rename(columns={"score":"score_cast"})
s129_df = dfs[15][["region","score"]]
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
plt.savefig("/home/ibai/pombo_johnny5/F123/figures/Figure1/rest/IS_scores_correlation_{}.pdf".format(date),dpi=300)


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
plt.savefig("/home/ibai/pombo_johnny5/F123/figures/Figure1/rest/IS_scores_distribution_{}.pdf".format(date),dpi=300)



# %%
