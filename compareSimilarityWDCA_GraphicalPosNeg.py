#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 09:58:37 2019

@author: benlitterer
"""

#TODO think about how to present gos that overlap with the blosum output's go's. Keep them in algo or not? 
#TODO also maybe consider not letting B62 match with itself and use that as a referance?? 
from goatools import obo_parser
import argparse
import os 
from collections import Counter
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import RdYlBu, Category20

go_obo = "/home/benlitterer/Academic/Research/ProjectMatrices"
go = obo_parser.GODag(go_obo + "/go-basic.obo")

parser = argparse.ArgumentParser(description="path to directory containing your Gene Annotation Output")
parser.add_argument("blast_directory")
parser.add_argument("go_directory")
args = parser.parse_args()
#cwd = "/home/benlitterer/Academic/Research/Summer2019/referance" 
goannots = {}

BLAST_OUT_DIR = args.blast_directory 

GO_dict = {}
db_path = args.go_directory #"/home/benlitterer/Academic/Research/Summer2019/testAllMatricesMediumInput/MappedGOs01Cutoff/"
def make_new_GOdict(): 
    for dirname in os.listdir(db_path): 
        curr_matrix = dirname[:-3] + ".out"
        if curr_matrix not in GO_dict.keys(): 
            GO_dict[curr_matrix] = {}
        for file in os.listdir(db_path + dirname): 
            if file not in GO_dict[curr_matrix].keys():
                GO_dict[curr_matrix][file] = []
            with open(db_path + dirname + "/" + file, "r") as curr_query: 
                query_lines = curr_query.readlines()
                for line in query_lines: 
                    if line[0] == "G": 
                        GO = line.strip("\n")
                        GO_dict[curr_matrix][file].append(GO)
            
make_new_GOdict()

print(GO_dict["MC20.out"]["XP_016885595.1GOs.out"])

def common_parent_go_ids(terms, go):
    '''
        This function finds the common ancestors in the GO 
        tree of the list of terms in the input.
    '''
    # Find candidates from first
    rec = go[terms[0]]
    candidates = rec.get_all_parents()
    candidates.update({terms[0]})
    
    # Find intersection with second to nth term
    for term in terms[1:]:
        rec = go[term]
        parents = rec.get_all_parents()
        parents.update({term})
        
        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)
    return candidates

def deepest_common_ancestor(terms, go):
    '''
        This function gets the nearest common ancestor 
        using the above function.
        Only returns single most specific - assumes unique exists.
    '''
    # Take the element at maximum depth. 
    try: 
        dca = max(common_parent_go_ids(terms, go), key=lambda t: go[t].depth)
        return dca
    except ValueError: 
        print("didn't work")

def min_branch_length(bl_go_id, kj_go_id, go):
    '''
        Finds the minimum branch length between two terms in the GO DAG.
    '''
    # First get the deepest common ancestor
    dca = deepest_common_ancestor([bl_go_id, kj_go_id], go)
    
    # Then get the distance from the DCA to each term
    dca_depth = go[dca].depth
    d1 = go[bl_go_id].depth - dca_depth
    d2 = go[kj_go_id].depth - dca_depth
    
    # get the total distance - i.e., to the deepest common ancestor and back.
    dist = d1 + d2
    
    # make the sign equal to the kj go - the blosum go 
    #print(f"{go[bl_go_id].depth}, {go[kj_go_id].depth}")
    if (go[bl_go_id].depth - go[kj_go_id].depth) <= 0: 
        dist = -dist
        
    return dist
    

print(min_branch_length("GO:0016021","GO:0016021",go))

#set output file
output_file("GO_dist_freqs_PosNeg.html")

#create figure 
p = figure(title="Distance of New GOs to BLOSUM62 GOs")

#set axis labels
p.xaxis.axis_label = "closest distance between new matrix annotation and B62's annotation"
p.yaxis.axis_label = '# of GO annotations'

#as many colors as there are matrices 
#color_palette = RdYlBu[len(GO_dict.keys())]
color_palette = Category20[len(GO_dict.keys())]
color_counter = 0

for matrix, queries in GO_dict.items(): 
    print(f"\n\t {matrix}")
    distances = []
    for query, GOs in queries.items(): 
        print(f"\t {query}")
        temp_distances = []
        #make sure BL62 returned something for this query 
        if query in GO_dict["BLOSUM62.out"].keys(): 
            if GOs: 
                for GO in GOs: 
                    #print("----------------------")
                    which_tree = go[GO].namespace
                    
                        
                    #the idea here is that blosum62 is being referanced against taking out the identity matches 
                    if matrix == "BLOSUM62.out":
                        distance_options = [min_branch_length(BL_GO, GO, go) for BL_GO in GO_dict["BLOSUM62.out"][query] if GO != BL_GO and go[BL_GO].namespace == which_tree and go[GO].name and go[BL_GO].name]
                        #print([min_branch_length(BL_GO, GO, go) for BL_GO in GO_dict["BLOSUM62.out"][query] if go[BL_GO].namespace == which_tree and GO != BL_GO and go[GO].name and go[BL_GO].name])
                    
                    #try and find all dca's for this go vs blosum62gos
                    else: 
                        #dca_options = [min_branch_length(BL_GO, GO ,go) for BL_GO in GO_dict["BLOSUM62.out"][query] if go[BL_GO].namespace == which_tree and go[GO].name and go[BL_GO].name] 
                        distance_options = [min_branch_length(BL_GO, GO, go) for BL_GO in GO_dict["BLOSUM62.out"][query] if go[BL_GO].namespace == which_tree and go[GO].name and go[BL_GO].name]
                    
                    if distance_options: 
                        closest_GO_dist = min(distance_options, key=abs)
                        #print(closest_GO_dist)
                        
                        distances.append(closest_GO_dist)
                        temp_distances.append(closest_GO_dist)
                
                        
                #this is for sorting by the vals 
                print("\t \t" + str(sorted(dict(Counter(temp_distances)).items(), key=lambda x: x[1], reverse=True)))
                    
    dist_distro = sorted(dict(Counter(distances)).items(), key=lambda x: x[0], reverse=True)
    
    #get all the counts for which the distance is positive
    pos_half = dict([pair for pair in dist_distro if pair[0] >= 0 ])
    print(pos_half)
    
    #get all the counts for which the distance is negative, flip the sign but make the count 
    neg_half = {abs(pair[0]): -pair[1] for pair in dist_distro if pair[0] < 0}
    print(neg_half)
    
    #p.line(x=[pair[0] for pair in dist_distro], y=[pair[1] for pair in dist_distro], color=color_palette[color_counter],muted_color=color_palette[color_counter], line_width=2, alpha=1, muted_alpha=.3, legend=matrix)
    p.line(x=[pair[0] for pair in pos_half.items()], y=[pair[1] for pair in pos_half.items()], color=color_palette[color_counter],muted_color=color_palette[color_counter], line_width=2, alpha=1, muted_alpha=.3, legend=matrix)
    p.line(x=[pair[0] for pair in neg_half.items()], y=[pair[1] for pair in neg_half.items()], color=color_palette[color_counter],muted_color=color_palette[color_counter], line_width=2, alpha=1, muted_alpha=.3, legend=matrix)
    color_counter += 1
    print(f"\t \t----------------------------------------\n\t \tTotal: {dist_distro}") 

#it says mute but actually the mute setting on the glyph increases the alpha value 
p.legend.click_policy = "hide"

show(p)
