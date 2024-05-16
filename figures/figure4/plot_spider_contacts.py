#%%

from scipy import interpolate
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

df = pd.read_csv("{}full_contacts_with_features_and_genes.tsv".format(root),sep="\t")

#%%
def coordinates_from_location_string (region):
    chrom = region.split(':')[0]
    if len (region.split(':')) == 1:
        start, stop = '', ''
    else:
        coordinates = region.split(':')[1]
        start = coordinates.split('-') [0]
        stop = coordinates.split ('-') [1]
        start = int(start.replace(',', ''))
        stop = int(stop.replace(',', ''))
    if start > stop:
        raise Exception ('Check coordinates of your region! It can not end before the start')
    return chrom, start, stop  

def get_diff_contacts(with_PRC=False):
    if with_PRC:
        cast_diff_df = df.query("delta >= 1 and bin1_PRC >= 1 and bin2_PRC >= 1")
        s129_diff_df = df.query("delta <= -1 and bin1_PRC >= 1 and bin2_PRC >= 1")
    else:
        cast_diff_df = df.query("delta >= 1")
        s129_diff_df = df.query("delta <= -1")

    return cast_diff_df,s129_diff_df

def get_contact_list(region,contacts_df,clean_region = False,max_dist=100000000):
    chrom, start, end = coordinates_from_location_string (region)
    contacts_df = contacts_df.loc[contacts_df.chrom1 == chrom]
    starts1 = contacts_df['start1'].values.tolist()
    ends1 = contacts_df['end1'].values.tolist()
    starts2 = contacts_df['start2'].values.tolist()
    ends2 = contacts_df['end2'].values.tolist()
    starts = []
    ends = []
    for n,bin_ in enumerate(zip(starts1,ends1,starts2,ends2)):
        bin_1 = bin_[0]
        bin_2 = bin_[1]
        bin_3 = bin_[2]
        bin_4 = bin_[3]
        if (start <= bin_1) and (start <= bin_2) and (end >= bin_3) and (end >= bin_4):
            if clean_region != False:
                bin_x, bin_y = clean_region[0],clean_region[1]
                if ((bin_3 - bin_1) < max_dist) and not((bin_1 > bin_x) and (bin_1 < bin_y)) and not((bin_3 > bin_x) and (bin_3 < bin_y)):
                    starts.append(bin_1)
                    ends.append(bin_3)
            else:
                if (bin_3 - bin_1) < max_dist:
                    starts.append(bin_1)
                    ends.append(bin_3)
    return starts,ends


## filters are filter the positions we want to show only. 
## highlight region is to color some regions with another color
## highluigh contacts are to color the contacts from the filters
## homotypic if we want to only connect same type of peaks, for example

def connect_points(ax,points1,points2,c,filter_starts = 0,filter_ends = 0,homotipic=False,highlight_regions=[],highlight_color="darkred",highlight_contacts=False):

    def is_in_region(pos):
        for reg in highlight_regions:
            start,end = reg.split("-")
            if (int(start) <= pos) and (int(end) >= pos):
                return True
        return False 
    

    #first filter the connections if needed
    if filter_starts != 0:
        filter_df = pd.DataFrame({"chrom":["chrZ"]*len(filter_starts),"start":filter_starts,"end":filter_ends})
        filter_pb = pb.BedTool.from_dataframe(filter_df)
        points1_ends = [x+30000 for x in points1]
        points1_df = pd.DataFrame({"chrom":["chrZ"]*len(points1),"start":points1,"end":points1_ends})
        points1_pb = pb.BedTool.from_dataframe(points1_df)
        points2_ends = [x+30000 for x in points2]
        points2_df = pd.DataFrame({"chrom":["chrZ"]*len(points2),"start":points2,"end":points2_ends})
        points2_pb = pb.BedTool.from_dataframe(points2_df)

        filtered_p1 = points1_pb.intersect(filter_pb,wa=True,wb=True).to_dataframe().start.values.tolist()
        filtered_p2 = points2_pb.intersect(filter_pb,wa=True,wb=True).to_dataframe().start.values.tolist()
        print ("Connections filtered")
        print (len(filtered_p1))
        print (len(filtered_p2))

    contacting_peaks = []
    for p1_aux,p2_aux in zip(points1,points2):
        p1 = 0
        p2 = 0
        rest_p1 = 0
        rest_p2 = 0
        if filter_starts != 0:
            if homotipic:
                if (p1_aux in filtered_p1) and (p2_aux in filtered_p2):
                    p1 = p1_aux
                    p2 = p2_aux
                else:
                    rest_p1 = p1_aux
                    rest_p2 = p2_aux
            else:
                if (p1_aux in filtered_p1):
                    p1 = p1_aux
                    p2 = p2_aux
                    contacting_peaks.append(p1)
                elif (p2_aux in filtered_p2):
                    p1 = p1_aux
                    p2 = p2_aux
                    contacting_peaks.append(p2)
                else:
                    rest_p1 = p1_aux
                    rest_p2 = p2_aux
        else:
            p1 = p1_aux
            p2 = p2_aux
        if p1 != 0:
            midpoint =(p2-p1)/2+p1 
            x = [p1,midpoint,p2]
            max = (p2-p1)/30000
            #if we want a fixed y
            #max = 100
            y= [0,max,0]
            x2 = np.linspace(x[0], x[-1], 100)
            y2 = interpolate.pchip_interpolate(x, y, x2)
            if len(highlight_regions) != 0:
                if is_in_region(p1) or is_in_region (p2):
                    ax.plot(x2, y2,highlight_color,alpha=1,linewidth=2)
                else:
                    ax.plot(x2, y2,c,alpha=1,linewidth=2)
            else:
                ax.plot(x2, y2,c,alpha=1,linewidth=2)
        if highlight_contacts:
            if rest_p1 != 0:
                midpoint =(rest_p2-rest_p1)/2+rest_p1 
                x = [rest_p1,midpoint,rest_p2]
                max = (rest_p2-rest_p1)/30000
                y= [0,max,0]
                x2 = np.linspace(x[0], x[-1], 100)
                y2 = interpolate.pchip_interpolate(x, y, x2)
                ax.plot(x2, y2,"black",alpha=0.0,linewidth=2)
    
    if filter_starts != 0:
        for f in contacting_peaks:
            ax.plot(f,0,color=c,marker="o")

# %%

# region for Figure 4
region = "chr13:21000000-24000000" 

fig,(ax1,ax2) = plt.subplots(2,1,figsize=(12,6),sharex=True)

cast_contacts_df, s129_contacts_df = get_diff_contacts(with_PRC=True)
starts_c,ends_c = get_contact_list(region,cast_contacts_df,clean_region=[22800000,23700000])
starts_s,ends_s,= get_contact_list(region,s129_contacts_df,clean_region=[22800000,23700000])


connect_points(ax1,starts_c,ends_c,cast_color)
connect_points(ax2,starts_s,ends_s,s129_color)
plt.savefig(save_path+"spider_plot_{}.pdf".format(region),dpi=100)
# %%

