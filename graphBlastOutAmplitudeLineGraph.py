#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:31:06 2019

@author: benlitterer
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 06:51:32 2019

@author: benlitterer
"""

import glob
import os
#import pandas as pd
import time
import argparse
import statistics
from math import pi
from bokeh.io import output_file, show
from bokeh.plotting import figure 
from bokeh.palettes import Viridis
from collections import OrderedDict
start_time = time.time()
parser = argparse.ArgumentParser(description="path to directory containing your BLAST output")
parser.add_argument("directory")
args = parser.parse_args()

#db_hits_list=[set(['a','b','z','q']), set(['a','b','c','f']), set(['a','c']), set(['a','c','f',]), set(['a','d'])]
#iterate over all .out files in current working directory
print("here are the files that contain output and the output length")
matrix_dict = {}
for filename in glob.glob(args.directory + '/*.out'):
    hit_dict = {}
    if 'slurm' not in filename and os.stat(filename).st_size > 0:
        matrix_name = filename.split("/")[1].strip(".out")
        """
        FIRST VERSION: this is 4x slower but loads into a pandas dataframe for easier data exploration
        print("FILE: " + filename)
        hit_set = set(pd.read_csv(filename, '\t', header=None).loc[:, 1])
        """
        #matrix_dict[matrix_name] = {}
        with open(filename) as file: 
            content = file.readlines()
        hit_list = []
        prev_qname = content[0].split('\t')[0]
        last_val = content[len(content)-1].split('\t')[0]
        for i in range(len(content)): 
            line = content[i]
            curr_hit = line.split('\t')[1]
            curr_qname = line.split('\t')[0]
            if curr_qname != prev_qname or i == len(content)-1:
                # the last element needs to be appended before the hit_set is added
                if i == len(content)-1:
                    hit_list.append(curr_qname)
                    
                hit_dict[prev_qname] = hit_list
                hit_list=[]
            hit_list.append(curr_hit)
            prev_qname = curr_qname
        matrix_dict[matrix_name] = hit_dict
"""
matrix_qlens = {key: len(value) for (key, value) in matrix_dict.items()}
matrix_qlens = sorted(matrix_qlens.items(), key = lambda x: -x[1])
"""
max_val_name = max(matrix_dict, key=lambda k: len(matrix_dict[k]))
#print(f"MAX VAL IS : {max_val_name}")
x_ticks = list(matrix_dict[max_val_name].keys())
#print(x_ticks)
output_file("graphBlastAmplitudeWDoubles.html")
fig = figure(title = "Hit count with doubles included", x_range = x_ticks, x_axis_label = "Blast query" , y_axis_label = "# of hits (including duplicates)", plot_width = 1000, plot_height=800)
set_fig = figure(title = "Hit count unique proteins only", x_range = x_ticks, x_axis_label = "Blast_query", y_axis_label = "# of hits (excluding duplicates)", plot_width = 1000, plot_height=800)
#^^^works 
#print(matrix_dict["MC22"].keys())
#print(Viridis[11])

#important to remember that this sorting is being done while using the data that keeps duplicate values 
matrix_list = sorted(matrix_dict.keys(), key=lambda x: statistics.median([len(matrix_dict[x][query]) for query in matrix_dict[x]]))

color_array = Viridis[11]
color_count = 0
for matrix in matrix_list: 
    len_dict = OrderedDict({query : len(matrix_dict[matrix][query]) for query in matrix_dict[matrix]})
    len_set_dict = OrderedDict({query : len(set(matrix_dict[matrix][query])) for query in matrix_dict[matrix]})
    fig.line(x=[key for key in len_dict], y=[item[1] for item in len_dict.items()], line_width=2, color=str(color_array[color_count]), alpha=.6, legend=matrix)
    set_fig.line(x=[key for key in len_set_dict], y=[item[1] for item in len_set_dict.items()], line_width=2, color=str(color_array[color_count]), alpha=.6, legend=matrix)
    color_count+=1
    
fig.legend.click_policy="hide"
fig.legend.title = "Matrices"
fig.xaxis.major_label_orientation = pi/3
fig.yaxis.major_label_orientation = pi/3

set_fig.legend.click_policy="hide"
set_fig.legend.title = "Matrices" 
set_fig.xaxis.major_label_orientation = pi/3
set_fig.yaxis.major_label_orientation = pi/3

show(fig)

output_file("graphBlastAmplitudeNoDoubles.html")
show(set_fig)

        



