#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:34:54 2019

@author: benlitterer
"""

import argparse
from collections import Counter 
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from goatools import obo_parser


go_obo = "/home/benlitterer/Academic/Research/ProjectMatrices"
go = obo_parser.GODag(go_obo + "/go-basic.obo")


parser = argparse.ArgumentParser(description="path to directory containing your Gene Annotation Output")
parser.add_argument("directory")
args = parser.parse_args()

"""
with open(args.directory, 'r') as file:
    data = file.read().replace('\n','')
"""

GO_dict = {}

"""
for db in data.split('^'): 
    if len(db) > 0: 
        entries = db.split("#")
        db_name = entries[0]
        keys=[entry.split("[")[0] for entry in entries]
        keys=keys[1:]
        vals=[]
        for entry in entries: 
            key_val_pairs = entry.split('[')
            if len(key_val_pairs)>1:
                vals.append([val.strip(" '") for val in key_val_pairs[1].strip().replace("]", "").split(',')])
        
        GO_dict[db_name] = dict(zip(keys,vals))

print(GO_dict)
"""

db_path = args.directory

GO_dict = {}

def make_new_GOdict(): 
    for dirname in os.listdir(db_path): 
        curr_matrix = dirname[:-3] + ".out"
        if curr_matrix not in GO_dict.keys(): 
            GO_dict[curr_matrix] = {}
        for file in os.listdir(db_path + dirname):
            curr_file = file[:-7]
            if file not in GO_dict[curr_matrix].keys():
                GO_dict[curr_matrix][curr_file] = {}
            with open(db_path + dirname + "/" + file, "r") as curr_query: 
                query_lines = curr_query.readlines()
                for i in range(len(query_lines)):
                    line = query_lines[i]
                    if line[0] == "G": 
                        GO = line.strip("\n")
                        assoc_prot_count = len(query_lines[i+1].split(","))
                        GO_dict[curr_matrix][curr_file][GO] = assoc_prot_count
            
make_new_GOdict()
        
LEVELS_TO_SEARCH = 8
def get_parents(GO_id, counter, parents): 
    if counter >= LEVELS_TO_SEARCH or parents == None: 
        return parents
    else: 
        parent_list = [parent.id for parent in go[GO_id].parents]
        if len(parent_list)!=0:
            for second_GO_id in parent_list: 
                #print(second_GO_id)
                parents.append(second_GO_id)
                #print("parents is currently equal to " + str(parents))
                my_parents = get_parents(second_GO_id, counter+1, parents)
            return parents+my_parents
                
        else: 
            return parents
        counter += 1
        
#ok to only iterate over the keys in one matrix because if the other one doesn't have the key then theres no overlap at all and vice versa
plt.rcParams.update({'font.size': 6})

print(GO_dict["MC30.out"])
# dict of dicts with outer keys being the query, inner key being the database, and innermost values being the children in question 

container_children_of_blast = {}
for key in GO_dict["BLOSUM62.out"]:
    print("QUERY : " + key)
    #list of counters corresponding to each db 
    count_list = []
    
    #column number changes based on whether the current database has the current key 
    cols = []
    for db, query in GO_dict.items(): 
        if key in query.keys(): 
            
            #IMPORTANT: shows that this is working correctly
            #print("DB: " + db + " KEY: " + str(key))
            cols.append(db)
            count_list.append(query[key])
    #print(len(count_list))
    print(cols)
    
    df = pd.DataFrame(count_list).T.fillna(value=0)
    df.columns = cols
    
    #force output to show all rows in dfo
    pd.set_option('display.max_rows', df.shape[0])
    
    #print first term_limit number of hits from BLOSUM
    term_limit = 6
    print("leading BLOSUM62 Hits")
    
    #should be a series I believe 
    leading_BL_terms = df[["BLOSUM62.out"]].sort_values(by=["BLOSUM62.out"],ascending=False).head(term_limit)
    print(leading_BL_terms)
    
    fig, axes = plt.subplots(nrows=math.ceil(len(cols)/4.0), ncols=3,squeeze=False)
    #which column and row to place graph in 
    col_count = 0
    row_count = 0
    print(count_list)
    
    children_of_blast = {}
    for col in df.columns: 
        #print("ROW: " + str(which_row) + " COL: " + str(which_col))
        if col != "BLOSUM62.out":
            #if you want to avoid seeing all of the empty space where neither column has hits 
            #df=df[(df["BLOSUM62.out"] != 0 ) | (df[col] != 0)]
            #print("the current df :")
            #print(df[["BLOSUM62.out", col]])
            
            #print first term_limit number of hits unique to kejues matrix 
            print("leading terms unique to " + col)
            
            leading_comparison_terms = df[(df["BLOSUM62.out"] == 0) & (df[col] != 0)].sort_values(by=[col], ascending=False)[[col]].head(term_limit)
            print(leading_comparison_terms)
            
            

            if not df.empty: 
                child_list = []
                for annot in leading_comparison_terms.index: 
                    parents = [annot]
                    parents = set(get_parents(annot, 0, parents))
                    
                    child_list = []
                    for parent in parents: 
                        if parent in leading_BL_terms.index: 
                            cp_statement = str(annot) + " @ level " + str(go[annot].level) + ":" + str(go[annot].name) +  " \nis a child of \n" + str(parent) + ":" + str(go[parent].name) + " @ level " + str(go[parent].level)
                            print(cp_statement)
                            child_list.append(cp_statement)
                children_of_blast[col] = child_list
                    
                splot = df[["BLOSUM62.out",col]].sort_values(by=["BLOSUM62.out",col],ascending=True)
                print(row_count, col_count)
                # use this param to tile graphs, but it seems redundant given the legend title = "BLOSUM62 vs " + str(col)
                l = splot.plot(kind="line", ax=axes[row_count,col_count], stacked=False,alpha=.75)
                l.set(xlabel="GO terms", ylabel="frequency")
        
                #moves plotting coordinates 
                if col_count == 2: 
                    col_count = 0
                    row_count +=1
                else: 
                    col_count+=1
    
    fig.suptitle("QUERY: " + str(key))
    plt.show()
    
    #add dictionary based on db keys to dictionary with query based keys 
    container_children_of_blast[key] = children_of_blast
    
# if at a point in time you need to access the output from kjs matrix that are children of the other matrix 
print(container_children_of_blast)  
