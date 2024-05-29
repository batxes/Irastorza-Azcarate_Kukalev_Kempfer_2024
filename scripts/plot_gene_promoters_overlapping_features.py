#%%
import pandas as pd
import matplotlib.pyplot as plt
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

# Figure 2h

master_df = pd.read_csv(root+"/300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed",sep="\t",index_col=None)

fig,ax = plt.subplots(figsize=(2,3))
# order  H3K4me3,ATAC, H3K27ax, CTCF, rad21 
x = ["H3K4me3","ATAC", "H3K27ac", "CTCF", "rad21"]
values1 = []
for n in ["n_cast_h3k4me3","n_cast_ATAC","n_cast_h3k27ac","n_cast_CTCF","n_cast_Rad21"]:
    values1.append(len(master_df.query('type == "cast" and {} >= 1'.format(n))))
print (values1)
ax.bar(x,values1,label="cast",color=cast_color)
values4 = []
for n in ["n_s129_h3k4me3","n_s129_ATAC","n_s129_h3k27ac","n_s129_CTCF","n_s129_Rad21"]:
    values4.append(len(master_df.query('type == "cast" and {} >= 1'.format(n))))
print (values4)
ax.bar(x,values4,label="s129",color=s129_color,bottom=values1)
values_aux = [x+y for (x,y) in zip(values1,values4)]
values2 = []
for n in ["n_common_h3k4me3","n_common_ATAC","n_common_h3k27ac","n_common_CTCF","n_common_Rad21"]:
    values2.append(len(master_df.query('type == "cast" and {} >= 1'.format(n))))
ax.bar(x,values2,label="common",color="silver",bottom=values_aux)
print (values2)
values3 = []
for n in ["n_nonphased_h3k4me3","n_nonphased_ATAC","n_nonphased_h3k27ac","n_nonphased_CTCF","n_nonphased_Rad21"]:
    values3.append(len(master_df.query('type == "cast" and {} >= 1'.format(n))))
print (values3)
values_aux2 = [x+y for (x,y) in zip(values_aux,values2)]
ax.bar(x,values3,label="nonphased",color=non_phased_color,bottom=values_aux2)
ax.axhline(1308,color=cast_color,linestyle="--",label="total cast genes")
ax.legend(loc="best")
ax.set_title("Features in CAST genes")
fig.savefig(save_path+"features_in_cast_genes.pdf",dpi=300)
print ("####")
fig,ax = plt.subplots(figsize=(2,3))
values1 = []
for n in ["n_cast_h3k4me3","n_cast_ATAC","n_cast_h3k27ac","n_cast_CTCF","n_cast_Rad21"]:
    values1.append(len(master_df.query('type == "s129" and {} >= 1'.format(n))))
ax.bar(x,values1,label="cast",color=cast_color)
print (values1)
values4 = []
for n in ["n_s129_h3k4me3","n_s129_ATAC","n_s129_h3k27ac","n_s129_CTCF","n_s129_Rad21"]:
    values4.append(len(master_df.query('type == "s129" and {} >= 1'.format(n))))
ax.bar(x,values4,label="s129",color=s129_color,bottom=values1)
print (values4)
values2 = []
for n in ["n_common_h3k4me3","n_common_ATAC","n_common_h3k27ac","n_common_CTCF","n_common_Rad21"]:
    values2.append(len(master_df.query('type == "s129" and {} >= 1'.format(n))))
values_aux2 = [x+y for (x,y) in zip(values1,values4)]
ax.bar(x,values2,label="common",color="silver",bottom=values_aux2)
print (values2)
values3 = []
for n in ["n_nonphased_h3k4me3","n_nonphased_ATAC","n_nonphased_h3k27ac","n_nonphased_CTCF","n_nonphased_Rad21"]:
    values3.append(len(master_df.query('type == "s129" and {} >= 1'.format(n))))
print (values3)
values_aux = [x+y for (x,y) in zip(values_aux2,values2)]
ax.bar(x,values3,label="nonphased",color=non_phased_color,bottom=values_aux)
ax.axhline(923,color=s129_color,linestyle="--",label="total s129 genes")
ax.legend(loc="best")
ax.set_title("Features in S129 genes")
fig.savefig(save_path+"features_in_s129_genes.pdf",dpi=300)



# %%

