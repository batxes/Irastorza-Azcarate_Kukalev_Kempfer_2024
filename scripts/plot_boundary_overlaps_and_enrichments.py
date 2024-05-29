#%%

"""
This script analyzes TAD boundary overlaps and the enrichment of different features in those overlaps. 
This plot was used in Figure 1g, 1h and SI figure 2e and 2f.

it plots upsetplots, enrichment profiles, heatmaps of those enrichments and gets significante between sets for insulation scores and features.
"""
# figure 1g, 1h, Suppl Figure 2e, 2f

import pandas as pd
import pybedtools as pb
import matplotlib.pyplot as plt
from upsetplot import UpSet
import numpy as np
import seaborn as sns
import os

import sys
import os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..' )
sys.path.append(mymodule_dir)
from variables import *

# Boundaries overlap. Intersection for any number of boundaries. Add 1 bp to all boundaries on each side

def read_beds (list_of_paths):
    import pandas as pd
    import pybedtools as pbt
    list_of_beds, list_of_extended_beds, list_of_pbts = [], [], []
    for file_path in list_of_paths:
        list_of_beds.append (pd.read_csv (file_path, sep='\t',names=["chrom","start","end"]))
    for bed in list_of_beds:
        bed_1bp =  extend_1bp (bed)
        list_of_extended_beds.append (bed_1bp)
        pbt_1bp = pbt.BedTool.from_dataframe (bed_1bp)
        list_of_pbts.append (pbt_1bp)
    concatenated_beds = pd.concat(list_of_extended_beds, axis=0)
    concatenated_beds_pbt = pbt.BedTool.from_dataframe (concatenated_beds)    
    return list_of_pbts, concatenated_beds_pbt

def extend_1bp (boundary_df, extra_distance=1): # add 1 bp to each side of each interval
    import pandas as pd
    import numpy as np
    extended_bed = boundary_df.copy()
    extended_bed.columns = ['chrom', 'start', 'stop']
    extended_bed.loc [:, 'new_start'] = extended_bed ['start'] - extra_distance
    extended_bed.loc [:, 'new_stop'] = extended_bed ['stop'] + extra_distance
    extended_bed = extended_bed [['chrom', 'new_start', 'new_stop']]
    extended_bed.columns = ['chrom', 'start', 'stop']
    return extended_bed

def save_temporal_beds (list_of_beds):
    list_of_bed_files = []
    for number, bed in enumerate (list_of_beds):
        file_name = 'bed_nr_' + str (number) + '.bed'
        bed.to_csv (file_name, sep='\t', index=False, header=False)
        list_of_bed_files.append (file_name)
    return list_of_bed_files

def save_beds (list_of_pbts, list_of_combinations, concatenated_beds_pbt, output_path):
    for pbt_file, combination in zip (list_of_pbts, list_of_combinations):
        combination_name = '_'.join (combination)
        overlap_pbt = concatenated_beds_pbt.intersect (pbt_file, wa=True)
        overlap_df = pb.BedTool.to_dataframe (overlap_pbt)
        if len (combination) == 1:
            overlap_df.to_csv (output_path + 'Intersect.unique.for.' + combination_name + '.bed', sep='\t', index=False, header=False)
        if len (combination) > 1:
            overlap_df.to_csv (output_path + 'Intersect.common.for.' + combination_name + '.bed', sep='\t', index=False, header=False)

def remove_temporal_beds (list_of_beds):
    import os
    for bed in list_of_beds:
        os.remove(bed)

def make_upset_plot (multiliner_boolean, filename=''):
    from upsetplot import plot
    boolean_df = (multiliner_boolean > 0)
    boolean_columns = ['_' + i for i in boolean_df.columns]
    boolean_df.columns = boolean_columns
    merged_df = boolean_df.join (multiliner_boolean)
    merged_df.set_index (boolean_columns, inplace=True)
    plt.figure (figsize=(15, 3))
    plot(merged_df, sort_by='cardinality', show_counts='%d') 
    if filename != '':
        plt.savefig (filename)
    plt.show()

def upset_plot_only_alleles (numbers_list, combinations_list, upset_file_name):
    from upsetplot import from_memberships
    upset_final_df = from_memberships(combinations_list ,data=numbers_list)
    facecolor="black"
    upset = UpSet(upset_final_df, 
                            #subset_size='count', 
                            intersection_plot_elements=3,
                            sort_by="cardinality",
                            show_counts='%d',
                            facecolor=facecolor)
    upset.style_subsets(present="CAST",
                        facecolor=cast_color,
                        label="cast")
    upset.style_subsets(present="S129",
                        facecolor=s129_color,
                        label="s129")
    upset.style_subsets(present=["CAST","S129"],
                        facecolor=both_color,
                        label="ase_only")
    fig = plt.figure(figsize=(40,80))
    upset.plot(fig=fig)
    fig.savefig("{}{}.pdf".format(save_path,upset_file_name),dpi=300)
    plt.show()

def upset_plot (numbers_list, combinations_list, upset_file_name):
    from upsetplot import from_memberships
    upset_final_df = from_memberships(combinations_list ,data=numbers_list)
    facecolor="black"
    upset = UpSet(upset_final_df, 
                            #subset_size='count', 
                            intersection_plot_elements=3,
                            sort_by="cardinality",
                            show_counts='%d',
                            facecolor=facecolor)
    upset.style_subsets(present="CAST",
                        facecolor=cast_color,
                        label="cast")
    upset.style_subsets(present="S129",
                        facecolor=s129_color,
                        label="s129")
    upset.style_subsets(present=["CAST","S129"],
                        facecolor=both_color,
                        label="ase_only")
    upset.style_subsets(present="Unphased",
                        facecolor=non_phased_color,
                        label="unphased")
    upset.style_subsets(present=["Unphased","CAST","S129"],
                        facecolor="black")
    upset.style_subsets(present=["Unphased","S129"],
                        facecolor="black")
    upset.style_subsets(present=["Unphased","CAST"],
                        facecolor="black")
    fig = plt.figure(figsize=(40,80))
    upset.plot(fig=fig)
    fig.savefig("{}{}.pdf".format(save_path,upset_file_name),dpi=300)
    plt.show()

def intersect_3way (list_of_paths, list_of_names, upset_file_name, output_path):
    name1, name2, name3 = list_of_names [0], list_of_names [1], list_of_names [2]
    list_of_pbts, concatenated_pbts = read_beds (list_of_paths)
    pbt1, pbt2, pbt3 = list_of_pbts [0], list_of_pbts [1], list_of_pbts [2]
    # Intersections
    pbt1_unique = (pbt1 - pbt2 - pbt3)
    pbt2_unique = (pbt2 - pbt1 - pbt3)
    pbt3_unique = (pbt3 - pbt1 - pbt2)
    pbt1_pbt2 = ((pbt1 + pbt2) - pbt3)
    pbt1_pbt3 = ((pbt1 + pbt3) - pbt2)
    pbt2_pbt3 = ((pbt2 + pbt3) - pbt1)
    pbt1_pbt2_pbt3 = (pbt1 + pbt2 + pbt3)
    # Upset plot and intersections
    numbers_list = [len (pbt1_pbt2_pbt3), len (pbt1_unique), len (pbt2_unique), len (pbt3_unique), len (pbt1_pbt2), len (pbt1_pbt3), len (pbt2_pbt3)]
    combinations_list = [[name1, name2, name3], [name1],[name2],[name3],[name1, name2],[name1, name3],[name2, name3]]
    upset_plot (numbers_list, combinations_list, upset_file_name)
    list_of_pbts = [pbt1_pbt2_pbt3, pbt1_unique, pbt2_unique, pbt3_unique, pbt1_pbt2, pbt1_pbt3, pbt2_pbt3]
    save_beds (list_of_pbts, combinations_list, concatenated_pbts, output_path)

def intersect_2way (list_of_paths, list_of_names, upset_file_name, output_path):
    name1, name2 = list_of_names [0], list_of_names [1]
    list_of_pbts, concatenated_pbts = read_beds (list_of_paths)
    pbt1, pbt2 = list_of_pbts [0], list_of_pbts [1]
    # Intersections
    pbt1_unique = pbt1 - pbt2
    pbt2_unique = pbt2 - pbt1
    pbt1_pbt2 = pbt1 + pbt2
    # Upset plot and intersections
    numbers_list = [len (pbt1_pbt2), len (pbt1_unique), len (pbt2_unique)]
    combinations_list = [[name1, name2], [name1],[name2]]
    upset_plot_only_alleles (numbers_list, combinations_list, upset_file_name)
    list_of_pbts = [pbt1_pbt2, pbt1_unique, pbt2_unique]
    save_beds (list_of_pbts, combinations_list, concatenated_pbts, output_path)

def multi_way_intersection (list_of_paths, list_of_names, output_path, plot_file_name=''):  
    import pybedtools as pbt
    if len (list_of_paths) == len (list_of_names):
        # files work
        list_of_beds = read_beds (list_of_paths)
        number_of_beds = len (list_of_beds)
        list_of_beds_1bp = extend_1bp (list_of_beds)
        temporal_bed_files = save_temporal_beds (list_of_beds_1bp)
        # intersections by multiliner
        x = pbt.BedTool()
        multiliner = x.multi_intersect(i=temporal_bed_files)
        multiliner_df = pbt.BedTool.to_dataframe (multiliner)
        multiliner_df_lenght = len (multiliner_df.columns)
        multiliner_boolean = multiliner_df.iloc [:, (multiliner_df_lenght - number_of_beds):multiliner_df_lenght]
        multiliner_boolean.columns = list_of_names
        remove_temporal_beds (temporal_bed_files)
        # show upset-plot
        make_upset_plot (multiliner_boolean, plot_file_name)
        # return bed-file coordinates
        concatenated_beds = pd.concat(list_of_beds_1bp, axis=0)
        concatenated_beds_pbt = pbt.BedTool.from_dataframe (concatenated_beds)
        # Rename combinations of intersections
        dict_keys = list(range(1, len (list_of_names) + 1))
        dict_keys = list(map(str, dict_keys))
        dict_beds_names = dict (zip (dict_keys, list_of_names))
        for key, value in dict_beds_names.items():
            multiliner_df['score'] = multiliner_df['score'].replace(key, value, regex=True)
        combinations = list (set (multiliner_df['score'].values))
        # Extract bed-files for each type of intersection
        for combination in combinations:
            extracted_coordinates = multiliner_df.query ('score == @combination')
            extracted_coordinates_pbt = pbt.BedTool.from_dataframe (extracted_coordinates)
            overlap_pbt = concatenated_beds_pbt.intersect (extracted_coordinates_pbt, wa=True)
            overlap_df = pbt.BedTool.to_dataframe (overlap_pbt)
            unique_or_common = combination.split (',')
            if len (unique_or_common) == 1:
                overlap_file_name = output_path + 'Unique_boundaries_for_' + str (combination) + '.bed'
            if len (unique_or_common) > 1:
                overlap_file_name = output_path + 'Common_boundaries_for_' + str (combination) + '.bed'
            overlap_df.to_csv (overlap_file_name, sep='\t', header=False, index=False)
    else:
        print ('Please check the number of arguments!')

def bed_intersect_merged (pbt1, pbt2):
    wa_common = pbt1.intersect(pbt2, wa=True, u=True)
    wb_common = pbt1.intersect(pbt2, wb=True, u=True)    
    pbt_merged = wa_common.cat(wb_common, postmerge=True)
    return pbt_merged

def bed_non_intersect (pbt1, pbt2):
    unique = pbt1.intersect(pbt2, wa=True, v=True)   
    return unique

def bedtools_intersect_two_way (df1, df2, df1_name, df2_name, filename=''):
    import pybedtools as pbt
    pbt1 = pbt.BedTool.from_dataframe(df1)
    pbt2 = pbt.BedTool.from_dataframe(df2)    
    common = bed_intersect_merged (pbt1, pbt2)
    pbt1_unique = bed_non_intersect (pbt1, pbt2)
    pbt2_unique = bed_non_intersect (pbt2, pbt1)
    numbers = [len (common), len (pbt1_unique), len (pbt2_unique)]
    from upsetplot import from_memberships
    upset = from_memberships ([[df1_name, df2_name], [df1_name],[df2_name]], data=numbers)
    plt.figure (figsize=(12, 8))
    upset.plot(upset, sort_by='cardinality', show_counts='%d')  
    if filename != '':
        plt.savefig (filename)
    plt.show()

def bedtools_intersect_three_way (df1, df2, df3, df1_name, df2_name, df3_name, filename=''):
    import pybedtools as pbt
    # Create pbts from dfs
    pbt1 = pbt.BedTool.from_dataframe(df1)
    pbt2 = pbt.BedTool.from_dataframe(df2)
    pbt3 = pbt.BedTool.from_dataframe(df3)
    pbt1_pbt2_common_temp = bed_intersect_merged (pbt1, pbt2) # Common for pbt1 and pbt2
    pbt1_pbt2_common = bed_non_intersect (pbt1_pbt2_common_temp, pbt3) # Only in pbt1 and pbt2
    pbt1_pbt3_common_temp = bed_intersect_merged (pbt1, pbt3)  # Common for pbt1 and pbt3
    pbt1_pbt3_common = bed_non_intersect (pbt1_pbt3_common_temp, pbt2) # Only for pbt1 and pbt3  
    pbt2_pbt3_common_temp = bed_intersect_merged (pbt2, pbt3) # Common for pbt2 and pbt3
    pbt2_pbt3_common = bed_non_intersect (pbt2_pbt3_common_temp, pbt1) # Only in pbt1 and pbt3   
    pbt1_pbt2_pbt3_common = bed_intersect_merged (pbt1_pbt2_common_temp, pbt3) # Common for all 3 
    print (len(pbt1_pbt2_pbt3_common.to_dataframe()))
    print ((pbt1_pbt2_pbt3_common.to_dataframe()))
    pbt1_pbt2_pbt3_common.to_dataframe().to_csv("~/test.tsv",sep="\t",index=False,header=None)
    pbt1_unique_temp = bed_non_intersect (pbt1, pbt2)
    pbt1_unique = bed_non_intersect (pbt1_unique_temp, pbt3)   # pbt1 unique
    pbt2_unique_temp = bed_non_intersect (pbt2, pbt3)
    pbt2_unique = bed_non_intersect (pbt2_unique_temp, pbt1)  # pbt2 unique
    pbt3_unique_temp = bed_non_intersect (pbt3, pbt1)
    pbt3_unique = bed_non_intersect (pbt3_unique_temp, pbt2)  # pbt3 unique  
    numbers = [len (pbt1_unique), len(pbt2_unique), len(pbt1_pbt2_common), len(pbt3_unique), len(pbt1_pbt3_common), len (pbt2_pbt3_common), len(pbt1_pbt2_pbt3_common)]

    from upsetplot import from_memberships
    upset_final_df = from_memberships ([[df1_name],[df2_name], [df1_name, df2_name], [df3_name], [df1_name, df3_name], [df2_name, df3_name], [df1_name, df2_name, df3_name]], data=numbers)
    print (upset_final_df)
    facecolor="black"
    upset = UpSet(upset_final_df, 
                            #subset_size='count', 
                            intersection_plot_elements=3,
                            sort_by="cardinality",
                            show_counts='%d',
                            facecolor=facecolor)
    upset.style_subsets(present="CAST",
                        facecolor=cast_color,
                        label="cast")
    upset.style_subsets(present="S129",
                        facecolor=s129_color,
                        label="s129")
    upset.style_subsets(present=["CAST","S129"],
                        facecolor=both_color,
                        label="ase_only")
    upset.style_subsets(present="Unphased",
                        facecolor=non_phased_color,
                        label="unphased")
    upset.style_subsets(present=["Unphased","CAST","S129"],
                        facecolor="black")
    upset.style_subsets(present=["Unphased","S129"],
                        facecolor="black")
    upset.style_subsets(present=["Unphased","CAST"],
                        facecolor="black")
    fig = plt.figure(figsize=(40,80))
    upset.plot(fig=fig)
    fig.savefig("{}upset_plot_TADs.pdf".format(save_path),dpi=300)
    plt.show()

    list_of_overlaps = [pbt1_unique, pbt2_unique, pbt1_pbt2_common, pbt3_unique, pbt1_pbt3_common, pbt2_pbt3_common, pbt1_pbt2_pbt3_common]
    return list_of_overlaps




#%%
# Figure 1g and SI Figure 2e
F123_all_as3NPs_TADs_winthoutX = '{}MiniIS.TAD.boundaries.F123.as3NPs.withoutX.bed'.format(root)
F123_all_as3NPs_CAST_TADs = '{}MiniIS.TAD.boundaries.F123.as3NPs.CAST.bed'.format(root)
F123_all_as3NPs_S129_TADs = '{}MiniIS.TAD.boundaries.F123.as3NPs.S129.bed'.format(root)

list_of_paths = [F123_all_as3NPs_S129_TADs,
                F123_all_as3NPs_CAST_TADs,
                F123_all_as3NPs_TADs_winthoutX]
list_of_names = ['S129', 'CAST', 'Unphased']

intersect_3way (list_of_paths, list_of_names, 'F123.boundaries.upset_supp', save_path)
intersect_2way ([F123_all_as3NPs_S129_TADs,F123_all_as3NPs_CAST_TADs], ['S129','CAST'], 'F123.boundaries.upset_main', save_path)
#multi_way_intersection (list_of_paths, list_of_names, save_path)

#%%

# Plot to get TAD size
cast = "{}TADs.F123.as3NPs.CAST.bed".format(root)
s129 = "{}TADs.F123.as3NPs.S129.bed".format(root)
unphased = "{}TADs.F123.as3NPs.bed".format(root)

cast_TADs = pd.read_csv(cast,sep="\t",names=["chrom","start","end"])
s129_TADs = pd.read_csv(s129,sep="\t",names=["chrom","start","end"])
unphased_TADs = pd.read_csv(unphased,sep="\t",names=["chrom","start","end"])
unphased_TADs = unphased_TADs[unphased_TADs.chrom.isin(chrs_)]
cast_TADs["size_"] = cast_TADs.end - cast_TADs.start
s129_TADs["size_"] = s129_TADs.end - s129_TADs.start
unphased_TADs["size_"] = unphased_TADs.end - unphased_TADs.start

fig,ax = plt.subplots(figsize=(4,4))
positions = np.arange(3) + 1
meanpointprops = dict(marker='o', markeredgecolor='black',
                      markerfacecolor='white',markersize=10)
bp = ax.boxplot([cast_TADs.size_, s129_TADs.size_, unphased_TADs.size_],
                showfliers=False,showmeans=True,patch_artist=True,positions=positions,notch=True,
                meanprops= meanpointprops)
colors = [cast_color,s129_color,non_phased_color]
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
plt.setp(bp["medians"],color="black")
ax.set_xticklabels(["CAST TADs","S129 TADs","unphased TADs"])
ax.set_ylabel("size in bps")
ax.set_title("Size of TADs")
plt.savefig("{}size_of_TADs.pdf".format(save_path),dpi=300)
print (np.mean(cast_TADs.size_))
print (np.mean(s129_TADs.size_))
print (np.mean(unphased_TADs.size_))

# %%
# Here we load all features that we want to overlap boundaries with.

# Load CTCF peaks
CTCF_all_peaks = pd.read_csv ('{}CTCF.R1R2.common.MACS2.broad.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
CTCF_CAST = pd.read_csv ('{}CTCF.CAST.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
CTCF_S129 = pd.read_csv ('{}CTCF.S129.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])

# Load Rad21 peaks
Rad21_all_peaks = pd.read_csv ('{}Rad21.R2R3.common.MACS2.broad.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
Rad21_CAST = pd.read_csv ('{}Rad21.CAST.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
Rad21_S129 = pd.read_csv ('{}Rad21.S129.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])

# Load H3K27Ac peaks
H3K27Ac_all_peaks = pd.read_csv ('{}H3K27Ac.R1R2.common.MACS2.broad.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
H3K27Ac_CAST = pd.read_csv ('{}H3K27Ac.CAST.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
H3K27Ac_S129 = pd.read_csv ('{}H3K27Ac.S129.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])

# Load H3K4me3 peaks
H3K4me3_all_peaks = pd.read_csv ('{}H3K4me3.R1R2.common.MACS2.broad.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
H3K4me3_CAST = pd.read_csv ('{}H3K4me3.CAST.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])
H3K4me3_S129 = pd.read_csv ('{}H3K4me3.S129.specific.peaks.bed'.format(root), skiprows=1, sep='\t', names=['chrom', 'start', 'stop'])

# Load ATAC-seq peaks
#ATAC_all_peaks = pd.read_csv ('/home/ibai/pombo_johnny5/F123/tracks/ATAC_seq/F123_ATAC_Bing_shared_peaks.bed', sep='\t', names=['chrom', 'start', 'stop'])
ATAC_all_peaks = pd.read_csv ('{}common_peaks.bed'.format(root), sep='\t', names=['chrom', 'start', 'stop'])
ATAC_CAST = pd.read_csv ('{}F123_bing_data030420_genome_CAST.bed'.format(root), sep='\t', names=['chrom', 'start', 'stop'])
ATAC_S129 = pd.read_csv ('{}F123_bing_data030420_genome_S129.bed'.format(root), sep='\t', names=['chrom', 'start', 'stop'])

# Load housekeeping genes
Housekeeping_genes = pd.read_csv ('{}ENSEMBLE.housekeeping.gene.list.bed'.format(root), sep='\t')
Housekeeping_genes_coordinates = Housekeeping_genes [['chrom', 'start', 'stop']]

# Load ENSEMBLE genes
ENSEMBLE_genes = pd.read_csv ('{}ENSEMBLE.GRCm38.99.from.gft.mm10.gene.list.bed'.format(root), sep='\t', names = ['chrom', 'start', 'stop', 'strand', 'ID', 'gene_name'])
ENSEMBLE_genes_coordinates = ENSEMBLE_genes [['chrom', 'start', 'stop']]

# Load all genes TSSs genes
all_genes = pd.read_csv ('{}300504.F123.most.expressed.isoform.promoter.states.TSS.allelic.reads.after.gene.filtering.transcript_biotype_updated_with_peaks.bed'.format(root), sep='\t')
all_genes_TSS_coordinates = all_genes [['chrom', 'TSS_start', 'TSS_end']]

# Load ASE genes
CAST_ASE_genes = pd.read_csv ('{}210504.F123.CAST.expressed.genes.bed'.format(root), sep='\t')
CAST_ASE_genes_coordinates = CAST_ASE_genes [['chrom', 'TSS_start', 'TSS_end']]

S129_ASE_genes = pd.read_csv ('{}210504.F123.S129.expressed.genes.bed'.format(root), sep='\t')
S129_ASE_genes_coordinates = S129_ASE_genes [['chrom', 'TSS_start', 'TSS_end']]

# Load SINEs
SINEs = pd.read_csv ('{}mm10.SINE.repeatmasker.txt.gz'.format(root), skiprows=0, sep='\t', names = ['chrom', 'start', 'stop', 'family', 'lenght', 'strand'])
SINEs_coordinates = SINEs [['chrom', 'start', 'stop']]

# Load LINEs
LINEs = pd.read_csv ('{}mm10.LINE.repeatmasker.txt.gz'.format(root), skiprows=0, sep='\t', names = ['chrom', 'start', 'stop', 'family', 'lenght', 'strand'])
LINEs_coordinates = LINEs [['chrom', 'start', 'stop']]

# Load S5p RNAPII
RNAPII_S5p = pd.read_csv ('{}F123.S5p.BCP_peaks.HM_mode.bed'.format(root), skiprows=0, sep='\t', names = ['chrom', 'start', 'stop', 'lenght', 'score'])
RNAPII_S5p_coordinates = RNAPII_S5p [['chrom', 'start', 'stop']]

# Load S7p RNAPII
RNAPII_S7p = pd.read_csv ('{}F123.S7p.BCP_peaks.HM_mode.bed'.format(root), skiprows=0, sep='\t', names = ['chrom', 'start', 'stop', 'lenght', 'score'])
RNAPII_S7p_coordinates = RNAPII_S7p [['chrom', 'start', 'stop']]

# Load K27me3 
K27me3 = pd.read_csv ('{}F123.K27me3.BCP_peaks.HM_mode.bed'.format(root), skiprows=0, sep='\t', names = ['chrom', 'start', 'stop', 'lenght', 'score'])
K27me3_coordinates = K27me3 [['chrom', 'start', 'stop']]

K27me3 = pd.read_csv ('{}F123.K27me3.BCP_peaks.HM_mode.bed'.format(root), skiprows=0, sep='\t', names = ['chrom', 'start', 'stop', 'lenght', 'score'])
K27me3_coordinates = K27me3 [['chrom', 'start', 'stop']]

lads_path = "{}mESC.mm10.liftover.LADs.bed".format(root)
lads = pd.read_csv(lads_path,sep="\t")[["chrom","start","end"]]
lads = lads[lads["chrom"].isin(chrs_)]

# Average profile function

def middle_boundary_point (boundaries_pbt):
    boundaries_df = pb.BedTool.to_dataframe (boundaries_pbt)
    chrom_list, mid_boundary_list = [], []
    for index, row in boundaries_df.iterrows ():
        chrom, start, stop = row [0], row [1], row [2]
        distance = stop - start
        mid_point = int (start + (distance / 2))  
        chrom_list.append (chrom)
        mid_boundary_list.append (mid_point)
    mid_boundary_df = pd.DataFrame (data=zip (chrom_list, mid_boundary_list), columns=['chrom', 'position'])
    return mid_boundary_df

def average_profile (boundaries_df, features_df, extra_distance, bin_size, merge_intervals=False):
    # Specify merge_intervals=True if you want to merge overlaping boundaries and 
    # use their central point as central area of average profiles 
    pbt_input = pb.BedTool.from_dataframe (boundaries_df)
    pbt_sorted = pbt_input.sort()
    if merge_intervals == True: # Merge overlaping intervals
        boundaries_pbt = pbt_sorted.merge()
    else:
        boundaries_pbt = pbt_sorted
    # Reformat and extend the coordinates in the original bed
    coordinates_df = middle_boundary_point (boundaries_pbt)
    coordinates_df ['start'] = coordinates_df ['position'] - extra_distance
    coordinates_df ['stop'] = coordinates_df ['position'] + extra_distance
    coordinates_df.drop (columns= ['position'], inplace=True)
    coordinates_df.reset_index (inplace=True)
    coordinates_df = coordinates_df.reindex(columns=['chrom', 'start', 'stop', 'index'])
    # Generate windows for each interval from bed
    coordinates_pbt = pb.BedTool.from_dataframe (coordinates_df)
    window_maker_pbt = pb.BedTool().window_maker(w=bin_size, i="src", b=coordinates_pbt)
    # Get the coverage of features per bin
    feature_pbt = pb.BedTool.from_dataframe (features_df)
    coverage_in_windows_pbt = window_maker_pbt.coverage (feature_pbt)
    coverage_in_windows_df = pb.BedTool.to_dataframe (coverage_in_windows_pbt)
    # Iterate over every interval from the original bed and collect all bins matching that interval
    start, stop = int (coordinates_df.head (n=1)['start']), int (coordinates_df.head (n=1)['stop'])
    heatmap_size = int ((stop - start) / bin_size)
    heatmap_df = pd.DataFrame (columns = range (heatmap_size), index=coordinates_df.index) # Empty dataframe
    intervals = coverage_in_windows_df['name'].unique() # array of intervals in the original bed 
    for coordinate in intervals:
        coordinate_subset = coverage_in_windows_df.query ('name==@coordinate')['score'].tolist()
        heatmap_df.loc [coordinate] = coordinate_subset
    # average_profiles
    sum_heatmap = heatmap_df.sum(axis = 0, skipna = True) 
    average_profile = sum_heatmap / len (heatmap_df)
    x_scale = np.arange (0, len (average_profile), 1)
    return x_scale, average_profile, heatmap_df

def plot_enrichment(ax1,ax2,hm,values_y,vmin,vmax):
    hm = hm.apply(pd.to_numeric)
    #CAST_unique_heatmap.columns = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    hm["sum"] = hm.sum(numeric_only=True, axis=1)
    hm = hm.sort_values("sum",ascending=False)
    hm = hm.drop(columns=["sum"])
    ax1.plot(values_y,color="black")
    #ax2.pcolor(hm,cmap="afmhot",vmin=0,vmax=10)
    ax2.imshow(hm,aspect="auto",cmap="viridis_r",vmin=vmin,vmax=vmax)
    #ax2.imshow(hm,aspect="auto",cmap="afmhot")
    plt.subplots_adjust(wspace=0.02, hspace=0.02)

def plot_enrichment_LAD(ax1,ax2,hm,values_y,vmin,vmax):
    hm = hm.apply(pd.to_numeric)
    hm["sum"] = hm.sum(numeric_only=True, axis=1)
    hm = hm.sort_values([40,"sum"],ascending=False)
    hm = hm.drop(columns=[40])
    hm = hm.drop(columns=["sum"])
    ax1.plot(values_y,color="black")
    #ax2.pcolor(hm,cmap="afmhot",vmin=0,vmax=10)
    ax2.imshow(hm,aspect="auto",cmap="viridis_r",vmin=vmin,vmax=vmax)
    #ax2.imshow(hm,aspect="auto",cmap="afmhot")
    plt.subplots_adjust(wspace=0.02, hspace=0.02)
# %%

# Script for SI Figure 2e
# This bit of the script plots enrichment plots for all boundaries from SI Figure 2e. 
# It also plots heatmaps for the each feature.

features = [CTCF_all_peaks,Rad21_all_peaks,H3K27Ac_all_peaks,H3K4me3_all_peaks,ATAC_all_peaks,Housekeeping_genes_coordinates,all_genes_TSS_coordinates,CAST_ASE_genes_coordinates,S129_ASE_genes_coordinates,SINEs_coordinates,LINEs_coordinates,RNAPII_S5p_coordinates,RNAPII_S7p_coordinates,K27me3_coordinates,lads, 
            CTCF_CAST,CTCF_S129,Rad21_CAST,Rad21_S129,H3K27Ac_CAST,H3K27Ac_S129,H3K4me3_CAST,H3K4me3_S129,ATAC_CAST,ATAC_S129,ENSEMBLE_genes_coordinates]
features_title = ["CTCF_all_peaks","Rad21_all_peaks","H3K27Ac_all_peaks","H3K4me3_all_peaks","ATAC_all_peaks","Housekeeping_genes_coordinates","all_genes_TSS_coordinates","CAST_ASE_genes_coordinates","S129_ASE_genes_coordinates","SINEs_coordinates","LINEs_coordinates","RNAPII_S5p_coordinates","RNAPII_S7p_coordinates","K27me3_coordinates","LADs", 
            "CTCF_CAST","CTCF_S129","Rad21_CAST","Rad21_S129","H3K27Ac_CAST","H3K27Ac_S129","H3K4me3_CAST","H3K4me3_S129","ATAC_CAST","ATAC_S129","ENSEMBLE_genes_coordinates"]
fig_big,ax_big = plt.subplots(len(features),7,figsize=(12,30),sharex=True)
colors=[non_phased_color,cast_color,s129_color,"black","black","black","black"]
paths = ['{}Intersect.unique.for.Unphased.bed'.format(root),
        '{}Intersect.unique.for.CAST.bed'.format(root),
        '{}Intersect.unique.for.S129.bed'.format(root),
        '{}Intersect.common.for.Unphased_CAST_S129.bed'.format(root),
        '{}Intersect.common.for.Unphased_CAST.bed'.format(root),
        '{}Intersect.common.for.Unphased_S129.bed'.format(root),
        '{}Intersect.common.for.CAST_S129.bed'.format(root)]
titles = ["Unphased","CAST","S129","shared","Unphased_CAST","Unphased_S129","Shared_CAST_S129"]

for j,feat in enumerate(features):
    plot_hm_data = []
    plot_y_data = []
    y_max = []
    y_min = []
    hm_max = []
    hm_min = []
    for n,x in enumerate(paths):
        max_aux_hm = []
        min_aux_hm = []
        if n < 3:
            df = pd.read_csv (x, sep='\t', names = ['chrom', 'start', 'stop'])
            #x_, y_, heatmap = average_profile (df, feat, 200000, 20000, True)
            x_, y_, heatmap = average_profile (df, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
        else:
            shared_pb = pb.BedTool(x)
            shared = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"])  # merge repeated boundaries 
            #x_, y_, heatmap = average_profile (shared, feat, 200000, 20000, True)
            x_, y_, heatmap = average_profile (shared, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
        max_aux_hm.append(heatmap.max())
        min_aux_hm.append(heatmap.min())

    flat = [item for sublist in plot_y_data for item in sublist]
    y_min = min(flat)
    y_max = max(flat)
    flat = [item for sublist in max_aux_hm for item in sublist]
    hm_max = max(flat)
    hm_max = np.percentile(flat,95)
    flat = [item for sublist in min_aux_hm for item in sublist]
    hm_min = min(flat)
    hm_min = np.percentile(flat,5)
    fig,axs = plt.subplots(2,len(plot_hm_data),figsize=(10,15), gridspec_kw={'height_ratios': [1, 7]},sharex=True)
    for n,x in enumerate(plot_hm_data):
        plot_enrichment(axs[0][n],axs[1][n],x,plot_y_data[n],hm_min,hm_max)
        axs[0][n].set_ylim([y_min,y_max])
        if n > 0:
            axs[0][n].set_yticks([])
        axs[0][n].set_title(titles[n])
        (row,col) = x.shape
        axs[0][n].set_xticks([0,int(col/2),col])
        axs[0][n].set_xticklabels([-1000,0,1000])
    fig.subplots_adjust(wspace=0.4, hspace=0.05)
    fig.suptitle(features_title[j])
    fig.savefig("{}boundary_enrichment_{}.pdf".format(save_path,features_title[j]),dpi=300)
    for n,x in enumerate(plot_y_data):
        ax_big[j][n].plot(x,color=colors[n])
        ax_big[j][n].set_ylim([y_min,y_max])
        ax_big[j][n].set_xticks([0,int(col/2),col])
        ax_big[j][n].set_xticklabels([-1000,0,1000])
        ax_big[j][0].set_ylabel(features_title[j],rotation=0)
        if n != 0:
            ax_big[j][n].set_yticks([])
        ax_big[j][n].set_xticks([])
        if j == len(features)-1:
            ax_big[j][n].set_xlabel(titles[n])
fig_big.savefig("{}boundary_enrichment_all_plots.pdf".format(save_path),dpi=300)



# %%

# This script bit is only focused on LADs. Plots the heatmap of the enrichment plots and clusters them. 
# Used in SI Figure 2f

features = [lads] # we can add more if needed.
features_title = ["LADs"]
colors=[cast_color,s129_color,non_phased_color]
changes_list = []
paths = ['{}Intersect.unique.for.CAST.main.bed'.format(root),
        '{}Intersect.unique.for.S129.main.bed'.format(root),
        '{}Intersect.common.for.CAST_S129.main.bed'.format(root)]
titles = ["CAST","S129","shared"]

for j,feat in enumerate(features):
    plot_hm_data = []
    plot_y_data = []
    plot_hm_data_sorted = []
    y_max = []
    y_min = []
    hm_max = []
    hm_min = []
    for n,x in enumerate(paths):
        max_aux_hm = []
        min_aux_hm = []
        if n < 3:
            df = pd.read_csv (x, sep='\t', names = ['chrom', 'start', 'stop'])
            x_, y_, heatmap = average_profile (df, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
        else:
            shared_pb = pb.BedTool(x)
            shared = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"])  # merge repeated boundaries 
            x_, y_, heatmap = average_profile (shared, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
        if features_title[j] == "LADs":
             # count number of LAD / interLAD changes in borders
            #changes = len(heatmap.query("(@heatmap[17] == 0 and @heatmap[16] == 0 and @heatmap[15] == 0 and @heatmap[24] == 1 and @heatmap[23] == 1 and @heatmap[22] == 1) or (@heatmap[17] == 1 and @heatmap[16] == 1 and @heatmap[15] == 1 and @heatmap[24] == 0 and @heatmap[23] == 0 and @heatmap[22] == 0)"))
            #less conservative
            changes = len(heatmap.query("(@heatmap[17] == 0 and @heatmap[16] == 0 and @heatmap[23] == 1 and @heatmap[22] == 1) or (@heatmap[17] == 1 and @heatmap[16] == 1 and @heatmap[23] == 0 and @heatmap[22] == 0)"))
            changes_index = list(heatmap.query("(@heatmap[17] == 0 and @heatmap[16] == 0 and @heatmap[23] == 1 and @heatmap[22] == 1) or (@heatmap[17] == 1 and @heatmap[16] == 1 and @heatmap[23] == 0 and @heatmap[22] == 0)").index)
            changes_list.append(list(heatmap.query("(@heatmap[17] == 0 and @heatmap[16] == 0 and @heatmap[23] == 1 and @heatmap[22] == 1) or (@heatmap[17] == 1 and @heatmap[16] == 1 and @heatmap[23] == 0 and @heatmap[22] == 0)").index))
            print ("Number of changes: {} from {}, {}%".format(changes,len(heatmap),changes/len(heatmap)))

        max_aux_hm.append(heatmap.max())
        min_aux_hm.append(heatmap.min())

        for i in range(len(heatmap)):
            if i in changes_index:
                heatmap.at[i, 40] = 1
            else: 
                heatmap.at[i, 40] = 0
        plot_hm_data_sorted.append(heatmap)

    flat = [item for sublist in plot_y_data for item in sublist]
    y_min = min(flat)
    y_max = max(flat)
    flat = [item for sublist in max_aux_hm for item in sublist]
    hm_max = max(flat)
    hm_max = np.percentile(flat,95)
    flat = [item for sublist in min_aux_hm for item in sublist]
    hm_min = min(flat)
    hm_min = np.percentile(flat,5)
    fig,axs = plt.subplots(2,len(plot_hm_data_sorted),figsize=(10,15), gridspec_kw={'height_ratios': [1, 7]},sharex=True)
    for n,x in enumerate(plot_hm_data_sorted):
        plot_enrichment_LAD(axs[0][n],axs[1][n],x,plot_y_data[n],hm_min,hm_max)
        axs[0][n].set_ylim([y_min,y_max])
        if n > 0:
            axs[0][n].set_yticks([])
        axs[0][n].set_title(titles[n])
        (row,col) = x.shape
        axs[0][n].set_xticks([0,int(col/2),col])
        axs[0][n].set_xticklabels([-1000,0,1000])
    fig.subplots_adjust(wspace=0.4, hspace=0.05)
    fig.suptitle(features_title[j])
    fig.savefig("{}boundary_enrichment_{}_main.pdf".format(save_path,features_title[j]),dpi=300)

#%%

# IS scores average profiles for already defined lowest centered boundaries

import math 
import scipy as scp

def get_significance_mann_whitney(x_,y_,text):
    x = [x for x in x_ if str(x) != 'nan']
    y = [y for y in y_ if str(y) != 'nan']
    stats, pvalue = scp.stats.mannwhitneyu(x,y,method="asymptotic")
    print  ("{}: {} {}".format(text,stats,pvalue))

def average_profile_of_IS_scores (TAD_boundary_df, IS_scores_df, square_size = "ins400K",bin_distance_from_lowest_score = 20):
    mm10_chrom_length = {'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116, 'chr5': 151834684, 'chr6':149736546,
                        'chr7': 145441459, 'chr8': 129401213, 'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022, 'chr13': 120421639,
                        'chr14': 124902244,'chr15': 104043685, 'chr16': 98207768, 'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566}
    chrom_list_mouse = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19']
    middle_IS_scores = []
    IS_scores_low_to_high_region_df = pd.DataFrame()
    #start_values = np.array (IS_scores_df.index.get_level_values('start')) # Find resolution
    #rint (IS_scores_df)
    start_values = IS_scores_df.start.values.tolist() # Find resolution
    resolution = start_values [1] - start_values [0] # Find resolution
    IS_scores_df = IS_scores_df.set_index(["chrom","start","stop"])
    IS_scores_df.sort_index(inplace=True) # Sort IS scores
    for row in TAD_boundary_df.itertuples(): # Iterate over TAD boundaries in supplied dataframe
        chromosome, start, stop = row [1], row [2], row [3]
        if chromosome in chrom_list_mouse: # Skip chromosomes except autosomal
            subset_IS_scores = IS_scores_df.loc [(chromosome, start):(chromosome, stop - resolution)]
            low_border = start - (bin_distance_from_lowest_score*resolution)
            high_border = stop + (bin_distance_from_lowest_score*resolution)

            if low_border > 0 and high_border < mm10_chrom_length.get (chromosome): # Remove TAD boundaries too close to start or end of the chromosome
                IS_scores_from_lowborder_to_highborder = IS_scores_df.loc [(chromosome, low_border):(chromosome, high_border - resolution)][square_size].tolist()
                #minimal_IS_score_value = min(IS_scores_from_lowborder_to_highborder) # Get lowest IS score value in the boundary
                minimal_IS_score_value = IS_scores_from_lowborder_to_highborder[21] # Take middle always
                middle_IS_scores.append(minimal_IS_score_value)
                #print (IS_scores_from_lowborder_to_highborder.index(minimal_IS_score_value),minimal_IS_score_value)
                #plt.plot(IS_scores_from_lowborder_to_highborder)
                #plt.show()
                column_name = chromosome + '-' + str(start) + '-' + str (stop)
                if not math.isnan(minimal_IS_score_value):
                    if len (IS_scores_from_lowborder_to_highborder) == (len (subset_IS_scores) + bin_distance_from_lowest_score*2):
                        IS_scores_from_lowborder_to_highborder = [x - minimal_IS_score_value for x in IS_scores_from_lowborder_to_highborder]
                        IS_scores_from_lowborder_to_highborder=IS_scores_from_lowborder_to_highborder[:43] #only for 20 left and righ
                        IS_scores_low_to_high_region_df [column_name] = IS_scores_from_lowborder_to_highborder
                    else:
                        empty_list = [np.nan] * (len (subset_IS_scores) + (bin_distance_from_lowest_score*2))
                        IS_scores_low_to_high_region_df [column_name] = empty_list  
     
    return IS_scores_low_to_high_region_df.T, middle_IS_scores

#def average_profile_of_IS_scores_normalized (TAD_boundary_df, IS_scores_df,square_size = "ins400K", bin_distance_from_lowest_score = 20):
#    mm10_chrom_lenght = {'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116, 'chr5': 151834684, 'chr6':149736546,
#                        'chr7': 145441459, 'chr8': 129401213, 'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022, 'chr13': 120421639,
#                        'chr14': 124902244,'chr15': 104043685, 'chr16': 98207768, 'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566}
#    chrom_list_mouse = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19']
#    IS_scores_low_to_high_region_df = pd.DataFrame()
#    #start_values = np.array (IS_scores_df.index.get_level_values('start')) # Find resolution
#    #resolution = start_values [1] - start_values [0] # Find resolution
#    start_values = IS_scores_df.start.values.tolist() # Find resolution
#    resolution = start_values [1] - start_values [0] # Find resolution
#    IS_scores_df = IS_scores_df.set_index(["chrom","start","stop"])
#    IS_scores_df.sort_index(inplace=True) # Sort IS scores
#    for row in TAD_boundary_df.itertuples(): # Iterate over TAD boundaries in supplied dataframe
#        chromosome, start, stop = row [1], row [2], row [3]
#        if chromosome in chrom_list_mouse: # Skip chromosomes except autosomal
#            subset_IS_scores = IS_scores_df.loc [(chromosome, start):(chromosome, stop - resolution)]
#            minimal_IS_score = subset_IS_scores[square_size].argmin() # Get the bin with lowest IS score in the boundary
#            minimal_IS_score_value = subset_IS_scores[square_size].min() # Get lowest IS score value in the boundary
#            if minimal_IS_score == minimal_IS_score: # Check for NAN values
#                print (minimal_IS_score)
#                chromosome_minIS, start_minIS, stop_minIS = minimal_IS_score [0], minimal_IS_score [1], minimal_IS_score [2]
#                low_border = start_minIS - (bin_distance_from_lowest_score*resolution)
#                high_border = stop_minIS + (bin_distance_from_lowest_score*resolution) 
#                if low_border > 0 and high_border < mm10_chrom_lenght.get (chromosome_minIS): # Remove TAD boundaries too close to start or end of the chromosome
#                    IS_scores_from_lowborder_to_highborder = IS_scores_df.loc [(chromosome, low_border):(chromosome, high_border - resolution)]['IS_score'].tolist()
#                    column_name = minimal_IS_score[0] + '-' + str(minimal_IS_score[1]) + '-' + str (minimal_IS_score[2])
#                    normalized_IS_Scores_from_lowborder_to_highborder = IS_scores_from_lowborder_to_highborder - minimal_IS_score_value # Normalize IS scores by setting the minimal value to zero
#                    IS_scores_low_to_high_region_df [column_name] = normalized_IS_Scores_from_lowborder_to_highborder
#    average_profile_values = IS_scores_low_to_high_region_df.mean(axis = 1, skipna = True).tolist()
#    return average_profile_values
#
def plot_is_scores(cast_df,s129_df,un_df,title,savename,ax):
    #fig,ax = plt.subplots(figsize=(2,2))
    ax.plot(cast_df.mean(),color=cast_color)
    ax.plot(s129_df.mean(),color=s129_color)
    ax.plot(un_df.mean(),color=non_phased_color)
    ax.set_xticks([1,11,21,31,41])
    ax.set_xticklabels([-1000,-500,0,500,1000])
    #ax.set_xlabel("Distance in Kbp")
    ax.set_title(title)
    ax.set_ylim([-0.03,0.20])
    #plt.savefig("/home/ibai/pombo_johnny5/F123/figures/Figure1/rest/{}.pdf".format(savename),dpi=300)

#%%

# Insulation aggregate plots for SI Figure2e

#this is supplementary. Main figure would be later.

paths = ['{}Intersect.unique.for.Nonphased.bed'.format(root),
            '{}Intersect.unique.for.CAST.bed'.format(root),
            '{}Intersect.unique.for.S129.bed'.format(root),
            '{}Intersect.common.for.S129_CAST_Nonphased.bed'.format(root),
            '{}Intersect.common.for.CAST_Nonphased.bed'.format(root),
            '{}Intersect.common.for.S129_Nonphased.bed'.format(root),
            '{}Intersect.common.for.S129_CAST.bed'.format(root)]

unique_unphased_df =  pd.read_csv (paths[0], sep='\t', names = ['chrom', 'start', 'stop'])
unique_cast_df =  pd.read_csv (paths[1], sep='\t', names = ['chrom', 'start', 'stop'])
unique_s129_df =  pd.read_csv (paths[2], sep='\t', names = ['chrom', 'start', 'stop'])
shared_pb = pb.BedTool(paths[3])
shared_unphased_cast_s129_df = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"])  
shared_pb = pb.BedTool(paths[4])
shared_unphased_cast_df = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"])  
shared_pb = pb.BedTool(paths[5])
shared_unphased_s129_df = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"])  
shared_pb = pb.BedTool(paths[6])
shared_cast_s129_df = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"]) 

fig,([[ax1,ax2],[ax3,ax4],[ax5,ax6],[ax7,ax8]]) = plt.subplots(4,2,figsize=(3,6),sharex=True,sharey=True)
cast_IS_scores = pd.read_csv("{}210910.F123.as.3NPs.curated.CAST.insulation.scores.at50Kb.table".format(root),sep="\t")
s129_IS_scores = pd.read_csv("{}210910.F123.as.3NPs.curated.S129.insulation.scores.at50Kb.table".format(root),sep="\t")
unphased_IS_scores = pd.read_csv("{}F123.as3NPs.IS.scores.csv".format(root),sep="\t")
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(unique_unphased_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(unique_unphased_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(unique_unphased_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Unique unphased","IS_plots_Unique_unphased",ax1)
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(unique_cast_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(unique_cast_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(unique_cast_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Unique cast","IS_plots_unique_CAST",ax2)
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(unique_s129_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(unique_s129_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(unique_s129_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Unique s129","IS_plots_unique_S129",ax3)
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(shared_unphased_cast_s129_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(shared_unphased_cast_s129_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(shared_unphased_cast_s129_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Shared all","IS_plots_shared_all",ax4)
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(shared_unphased_cast_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(shared_unphased_cast_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(shared_unphased_cast_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Shared unphased cast","IS_plots_shared_unphased_cast",ax5)
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(shared_unphased_s129_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(shared_unphased_s129_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(shared_unphased_s129_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Shared unphased s129","IS_plots_shared_unphased_s129",ax6)
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(shared_cast_s129_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(shared_cast_s129_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(shared_unphased_cast_s129_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Shared cast s129","IS_plots_shared_cast_s129",ax7)
ax7.set_xlabel("Distance to boundary in Kbp")
ax8.set_xlabel("Distance to boundary in Kbp")
ax1.set_ylabel("Normalized Insulation Score")
ax3.set_ylabel("Normalized Insulation Score")
ax5.set_ylabel("Normalized Insulation Score")
ax7.set_ylabel("Normalized Insulation Score")

plt.savefig("{}{}.pdf".format(save_path,"ALL_IS_plots.pdf"),dpi=300)

# %%

# Insulation aggregate plots for Figure1h

paths = ['{}Intersect.unique.for.CAST.main.bed'.format(root),
        '{}Intersect.unique.for.S129.main.bed'.format(root),
        '{}Intersect.common.for.S129_CAST.main.bed'.format(root)]

shared_df =  pd.read_csv (paths[2], sep='\t', names = ['chrom', 'start', 'stop'])
unique_cast_df =  pd.read_csv (paths[0], sep='\t', names = ['chrom', 'start', 'stop'])
unique_s129_df =  pd.read_csv (paths[1], sep='\t', names = ['chrom', 'start', 'stop'])


fig,([ax1,ax2,ax3]) = plt.subplots(1,3,figsize=(6,2),sharex=True,sharey=True)
cast_IS_scores = pd.read_csv("{}210910.F123.as.3NPs.curated.CAST.insulation.scores.at50Kb.table".format(root),sep="\t")
s129_IS_scores = pd.read_csv("{}210910.F123.as.3NPs.curated.S129.insulation.scores.at50Kb.table".format(root),sep="\t")
unphased_IS_scores = pd.read_csv("{}F123.as3NPs.IS.scores.csv".format(root),sep="\t")

df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(shared_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(shared_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(shared_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Shared","IS_plots_shared",ax1)
get_significance_mann_whitney(cast_middle_scores,s129_middle_scores,"cast_vs_s129")
get_significance_mann_whitney(cast_middle_scores,unphased_middle_scores,"cast_vs_unphased")
get_significance_mann_whitney(unphased_middle_scores,s129_middle_scores,"s129_vs_unphased")

df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(unique_cast_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(unique_cast_df,s129_IS_scores)
df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(unique_cast_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
get_significance_mann_whitney(cast_middle_scores,s129_middle_scores,"cast_vs_s129")
get_significance_mann_whitney(cast_middle_scores,unphased_middle_scores,"cast_vs_unphased")
get_significance_mann_whitney(unphased_middle_scores,s129_middle_scores,"s129_vs_unphased")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Unique cast","IS_plots_unique_CAST",ax2)
df_cast_scores,cast_middle_scores = average_profile_of_IS_scores(unique_s129_df,cast_IS_scores)
df_s129_scores,s129_middle_scores = average_profile_of_IS_scores(unique_s129_df,s129_IS_scores)

df_un_scores,unphased_middle_scores = average_profile_of_IS_scores(unique_s129_df,unphased_IS_scores,square_size = "as3NPs_ins400K")
plot_is_scores(df_cast_scores,df_s129_scores,df_un_scores,"Unique s129","IS_plots_unique_S129",ax3)
get_significance_mann_whitney(cast_middle_scores,s129_middle_scores,"cast_vs_s129")
get_significance_mann_whitney(cast_middle_scores,unphased_middle_scores,"cast_vs_unphased")
get_significance_mann_whitney(unphased_middle_scores,s129_middle_scores,"s129_vs_unphased")

ax2.set_xlabel("Distance to boundary in Kbp")
ax1.set_ylabel("Normalized Insulation Score")
ax3.set_ylabel("Normalized Insulation Score")
ax2.set_ylabel("Normalized Insulation Score")

plt.savefig("{}{}.pdf".format(save_path,"Main_IS_plots.pdf"),dpi=300)
#%%

# Here we do random permutations to see if the features in the boundaries are significantly enriched.


# Some changes to standarize the data
CTCF_all_peaks = CTCF_all_peaks.rename(columns={"stop":"end"})
Rad21_all_peaks = Rad21_all_peaks.rename(columns={"stop":"end"})
ATAC_all_peaks = ATAC_all_peaks.rename(columns={"stop":"end"})
Housekeeping_genes_coordinates = Housekeeping_genes_coordinates.rename(columns={"stop":"end"})
lads = lads.rename(columns={"stop":"end"})

import random 
random.seed(2)
chrs_ = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19']
CHROM_SIZES="{}mm10.chrom.sizes".format(root)
chrom_size = {}
with open(CHROM_SIZES) as f:
    for line in f:
        (key, val) = line.split()
        chrom_size[key] = int(val)

def permute(df,start="start",end="end"):
    s = start
    e = end
    df = df[df["chrom"].isin(chrs_)]
    shift = random.randint(1,20000000) #shift between 1 and 20 mega
    print (shift)
    df[s] = df[s] + shift
    df[e] = df[e] + shift
    new_start_values = []
    new_end_values = []
    new_chroms = []
    for chr_ in chrs_:
        size = chrom_size[chr_]
        aux_df = df[df["chrom"] == chr_].copy()
        start_values = aux_df[s].values.tolist()
        end_values = aux_df[e].values.tolist()
        for start,end in zip(start_values,end_values): #if they are bigger than chr, 
            if start > size:
                start = start - size + 3000000
                end = end - size + 3000000
            elif end > size:
                start = start - size + 3000000
                end = end - size + 3000000
            if start < 0:
                new_start_values.append(0)
            else:
                new_start_values.append(start)
            new_end_values.append(end)
            new_chroms.append(chr_)
    df["chrom"] = new_chroms 
    df[s] = new_start_values
    df[e] = new_end_values
    return df

number_of_averages = 10 #for testing. Add more.
permuted_CTCF_list,permuted_RAD21_list,permuted_ATAC_list, permuted_HK_list, permuted_LAD_list = [],[],[],[], []
for i in range(number_of_averages):
    permuted_CTCF = permute(CTCF_all_peaks,"start","end")
    permuted_RAD21 = permute(Rad21_all_peaks,"start","end")
    permuted_ATAC = permute(ATAC_all_peaks,"start","end")
    permuted_HK = permute(Housekeeping_genes_coordinates,"start","end")
    permuted_LAD = permute(lads,"start","end")
    permuted_CTCF_list.append(permuted_CTCF)
    permuted_RAD21_list.append(permuted_RAD21)
    permuted_ATAC_list.append(permuted_ATAC)
    permuted_HK_list.append(permuted_HK)
    permuted_LAD_list.append(permuted_LAD)

features_permuted = [permuted_CTCF_list,permuted_RAD21_list,permuted_ATAC_list,permuted_HK_list,permuted_LAD_list]
features = [CTCF_all_peaks,Rad21_all_peaks,ATAC_all_peaks,Housekeeping_genes_coordinates,lads]
features_title = ["permuted_CTCF","permuted_Rad21","permuted_ATAC","permuted_HK","permuted_LADs"]

fig_big,ax_big = plt.subplots(len(features),7,figsize=(12,4),sharex=True)
colors=[non_phased_color,cast_color,s129_color,"black","black","black","black"]
for j,feat in enumerate(features):
    print (features_title[j])
    paths = ['{}Intersect.unique.for.Nonphased.bed'.format(root),
             '{}Intersect.unique.for.CAST.bed'.format(root),
             '{}Intersect.unique.for.S129.bed'.format(root),
             '{}Intersect.common.for.S129_CAST_Nonphased.bed'.format(root),
             '{}Intersect.common.for.CAST_Nonphased.bed'.format(root),
             '{}Intersect.common.for.S129_Nonphased.bed'.format(root),
             '{}Intersect.common.for.S129_CAST.bed'.format(root)]

    titles = ["Unphased","CAST","S129","shared","Unphased_CAST","Unphased_S129","Shared_CAST_S129"]

    plot_hm_data = []
    plot_y_data = []
    plot_y_data_permuted = []
    y_max = []
    y_min = []
    hm_max = []
    hm_min = []
    for n,x in enumerate(paths):
        print ("Boundary number {}".format(n))
        max_aux_hm = []
        min_aux_hm = []
        aux_list = []
        if n < 3:
            df = pd.read_csv (x, sep='\t', names = ['chrom', 'start', 'stop'])
            x_, y_, heatmap = average_profile (df, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
            for i in range(number_of_averages):
                x_permuted, y_permuted, heatmap_permuted = average_profile (df, features_permuted[j][i], 1000000, 50000, True)
                aux_list.append(y_permuted)
            df_aux = pd.concat(aux_list,axis=1)
            df_aux["avg"] = df_aux.mean(axis=1)
            plot_y_data_permuted.append(df_aux["avg"])
            
        else:
            shared_pb = pb.BedTool(x)
            shared = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"])  # merge repeated boundaries 
            x_, y_, heatmap = average_profile (shared, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
            for i in range(number_of_averages):
                x_permuted, y_permuted, heatmap_permuted = average_profile (df, features_permuted[j][i], 1000000, 50000, True)
                aux_list.append(y_permuted)
            df_aux = pd.concat(aux_list,axis=1)
            df_aux["avg"] = df_aux.mean(axis=1)
            plot_y_data_permuted.append(df_aux["avg"])
        max_aux_hm.append(heatmap.max())
        min_aux_hm.append(heatmap.min())

    flat = [item for sublist in plot_y_data for item in sublist]
    y_min = min(flat)
    y_max = max(flat)
    flat = [item for sublist in max_aux_hm for item in sublist]
    hm_max = max(flat)
    hm_max = np.percentile(flat,95)
    flat = [item for sublist in min_aux_hm for item in sublist]
    hm_min = min(flat)
    hm_min = np.percentile(flat,5)
    for n,x in enumerate(plot_y_data):
        print ("Plotting {}".format(n))

        ax_big[j][n].plot(plot_y_data_permuted[n],color="silver",linestyle='dashed')
        ax_big[j][n].plot(x,color=colors[n])
        ax_big[j][n].set_ylim([y_min,y_max])
        ax_big[j][n].set_xticks([0,int(col/2),col])
        ax_big[j][n].set_xticklabels([-1000,0,1000])
        ax_big[j][0].set_ylabel(features_title[j],rotation=0)
        if n != 0:
            ax_big[j][n].set_yticks([])
        ax_big[j][n].set_xticks([])
        if j == len(features)-1:
            ax_big[j][n].set_xlabel(titles[n])
fig_big.savefig("{}boundary_enrichment_permuted_all_plots.pdf".format(save_path),dpi=300)

# %%

# Now the same but for main Figure 1h

fig_big,ax_big = plt.subplots(len(features),7,figsize=(12,4),sharex=True)
colors=[cast_color,s129_color,non_phased_color]
for j,feat in enumerate(features):
    print (features_title[j])
    paths = ['{}Intersect.unique.for.CAST.main.bed'.format(root),
             '{}Intersect.unique.for.S129.main.bed'.format(root),
             '{}Intersect.common.for.S129_CAST.main.bed'.format(root)]

    titles = ["CAST","S129","shared"]

    plot_hm_data = []
    plot_y_data = []
    plot_y_data_permuted = []
    y_max = []
    y_min = []
    hm_max = []
    hm_min = []
    for n,x in enumerate(paths):
        print ("Boundary number {}".format(n))
        max_aux_hm = []
        min_aux_hm = []
        aux_list = []
        if n < 3:
            df = pd.read_csv (x, sep='\t', names = ['chrom', 'start', 'stop'])
            x_, y_, heatmap = average_profile (df, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
            for i in range(number_of_averages):
                x_permuted, y_permuted, heatmap_permuted = average_profile (df, features_permuted[j][i], 1000000, 50000, True)
                aux_list.append(y_permuted)
            df_aux = pd.concat(aux_list,axis=1)
            df_aux["avg"] = df_aux.mean(axis=1)
            plot_y_data_permuted.append(df_aux["avg"])
            
        else:
            shared_pb = pb.BedTool(x)
            shared = shared_pb.sort().merge().to_dataframe(names=["chrom","start","stop"])  # merge repeated boundaries 
            x_, y_, heatmap = average_profile (shared, feat, 1000000, 50000, True)
            plot_y_data.append(y_)
            plot_hm_data.append(heatmap)
            for i in range(number_of_averages):
                x_permuted, y_permuted, heatmap_permuted = average_profile (df, features_permuted[j][i], 1000000, 50000, True)
                aux_list.append(y_permuted)
            df_aux = pd.concat(aux_list,axis=1)
            df_aux["avg"] = df_aux.mean(axis=1)
            plot_y_data_permuted.append(df_aux["avg"])
        max_aux_hm.append(heatmap.max())
        min_aux_hm.append(heatmap.min())

    flat = [item for sublist in plot_y_data for item in sublist]
    y_min = min(flat)
    y_max = max(flat)
    flat = [item for sublist in max_aux_hm for item in sublist]
    hm_max = max(flat)
    hm_max = np.percentile(flat,95)
    flat = [item for sublist in min_aux_hm for item in sublist]
    hm_min = min(flat)
    hm_min = np.percentile(flat,5)
    for n,x in enumerate(plot_y_data):
        print ("Plotting {}".format(n))

        ax_big[j][n].plot(plot_y_data_permuted[n],color="silver",linestyle='dashed')
        ax_big[j][n].plot(x,color=colors[n])
        ax_big[j][n].set_ylim([y_min,y_max])
        # uncomment this if we run the whole script from the beginning
        #ax_big[j][n].set_xticks([0,int(col/2),col])
        #ax_big[j][n].set_xticklabels([-1000,0,1000])
        ax_big[j][0].set_ylabel(features_title[j],rotation=0)
        if n != 0:
            ax_big[j][n].set_yticks([])
        ax_big[j][n].set_xticks([])
        if j == len(features)-1:
            ax_big[j][n].set_xlabel(titles[n])
fig_big.savefig("{}boundary_enrichment_permuted_all_plots_main.pdf".format(save_path),dpi=300)



# %%
