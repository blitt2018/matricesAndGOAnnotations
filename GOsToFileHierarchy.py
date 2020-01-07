#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:00:36 2019

@author: benlitterer
"""

import glob
import os
from multiprocessing import Pool 
import time
import argparse
import requests

start_time = time.time()
parser = argparse.ArgumentParser(description="path to directory containing your BLAST output")
parser.add_argument("directory")
parser.add_argument("output")
args = parser.parse_args()

""" ADJUST THESE FOR YOUR PIPELINE"""

EVAL_CUTOFF =  .01                                                     #keep evals smaller or equal to this number 

EVAL_COLUMN_NUMBER = 3                                                 #which column your evals are in the first one being 1 rather than zero 

OUTPUT_DIRECTORY = args.output

#create directory to put mapped protein db in 
new_dir_name = args.directory + OUTPUT_DIRECTORY
if not os.path.exists(new_dir_name):
    os.mkdir(new_dir_name)
    
#db_hits_list=[set(['a','b','z','q']), set(['a','b','c','f']), set(['a','c']), set(['a','c','f',]), set(['a','d'])]
db_hit_list=[]
#iterate over all .out files in current working directory
#print("here are the files that contain output and the output length")
file_list = []

for filename in glob.glob(args.directory + '/*.out'):
    #print("\n \n" + filename + "\n \n")
    hit_dict = {}
    if 'slurm' not in filename and os.stat(filename).st_size > 0:
        #print(filename)
        """
        FIRST VERSION: this is 4x slower but loads into a pandas dataframe for easier data exploration
        print("FILE: " + filename)
        hit_set = set(pd.read_csv(filename, '\t', header=None).loc[:, 1])
        """
        file_list.append(filename)
        with open(filename) as file: 
            content = file.readlines()
        """
        Here it is shown that there can be duplicates in the blast output. The algorithm only takes off duplicates within hits from one query, so duplicates that correspond to
        another query are maintained. 
        print(len([line.split('\t')[1] for line in content]))
        print(len(set([line.split('\t')[1] for line in content])))
        """
        hit_list = []
        prev_qname = content[0].split('\t')[0]
        last_val = content[len(content)-1].split('\t')[0]
        for i in range(len(content)): 
            line = content[i]
            curr_hit = line.split('\t')[1]
            curr_qname = line.split('\t')[0]
            curr_eval = float(line.split('\t')[EVAL_COLUMN_NUMBER-1])
            if curr_qname != prev_qname or i == len(content)-1:
                # the last element needs to be appended before the hit_set is added
                if i == len(content)-1:
                    hit_list.append(curr_qname)
                    
                hit_set = set(hit_list)
                """
                Duplicates that are responses to the same query should be alright to count as one hit. This code counts the amount of duplicates per set
                print(str(len(hit_list) - len(hit_set)) + " duplicates in this set" )
                """
                hit_dict[prev_qname] = hit_set
                hit_list=[]
            if curr_eval <= EVAL_CUTOFF: 
                #print("\t HIT: " + curr_hit + " EVAL: " + str(curr_eval))
                hit_list.append(curr_hit)
            prev_qname = curr_qname
        db_hit_list.append(hit_dict)
#a couple random tests 
#print(db_hit_list[3]['XP_006526966.1'])
#print(db_hit_list[2]['XP_017173608.1'])

#generate all values for each query result 
#seems to be unneccesary
"""
all_vals_dict = {}
for hit_dict in db_hit_list:
    for key in hit_dict.keys():
        all_vals_dict[key] = hit_dict[key]
"""

def call_uni_api(query_str): 
    payload = {
            # TESTING "query":"YP_220563.1 OR YP_220550.1"
            "query":query_str,
            "format":"tab",
            "columns":"genes,id,go-id, entry",
            "include":"yes",
            "compress":"no",
            "format":"tab"
            }
    
    r =requests.get("https://www.uniprot.org/uniprot/", params=payload)
    return r

print(db_hit_list)
print(file_list)
                  

for i in range(len(db_hit_list)):
    hit_dict = db_hit_list[i]
    associated_file = file_list[i]
    dont_include_hit_list = []

    #generate a list of sets that doesn't include set being currently iterated over 
    # try using dont_include_hit_list = [d for d in db_hit_list if d != hit_dict]
    for second_hit_dict in db_hit_list: 
        if hit_dict!=second_hit_dict: 
            dont_include_hit_list.append(second_hit_dict)
    
    # name of directory to output to 
    out_dir = file_list[i].split("/")[-1].split(".")[0]
    
    print("^" + out_dir)

    #gets all unique values in this dictionary and adds them to a list 
    for key in hit_dict: 
        #name of file to write output to 
        out_file = str(key) 
        
        print("QUERY: " + key)
        go_key_dict = {}
        corresponding_sets = [d[key] for d in dont_include_hit_list if key in d.keys()]
        unique_vals = list(hit_dict[key])
        #print("amt of queries " + str(len(unique_vals)))
        print(unique_vals)

        if len(unique_vals) > 0: 
            pool_args=[]
            for i in range(len(unique_vals)):
                curr_str = unique_vals[i]
                #print(curr_str)
                if len(pool_args) == 4 or i == len(unique_vals)-1: 
                    if i == len(unique_vals)-1: 
                        pool_args.append(curr_str)
                    p = Pool(4)
                    pool_out_list = p.map(call_uni_api, pool_args)
                    p.close()
                    p.join()
                    print("call being done on: " + str(pool_args))
                    for j in range(len(pool_out_list)): 
                        r = pool_out_list[j]
                        hit_str = pool_args[j]
                        line_list = [hit_str]
                        #extend doesn't return a value because lists are mutable 
                        lines = r.text.split('\n')
                        
                        #always want to get second line because the first is simply the column labels 
                        if len(lines)>1:
                            
                            line_list.extend(lines[1].split('\t'))
                            #print(line_list)
                            #print(line_list)
                            GOs = [item.strip(' ') for item in line_list[-1].split(';')]
                            #print(hit_str)
                            #print(GOs)
                            for GO in GOs: 
                                if GO in go_key_dict.keys(): 
                                    go_key_dict[GO].append(hit_str)
                                else: 
                                    go_key_dict[GO] = [hit_str]
                    pool_args=[curr_str]
                else: 
                    pool_args.append(curr_str)

            output_str = ""
            for key, val in go_key_dict.items(): 
                output_str += key + "\n"
                output_str += str(val) + "\n" 
                
                new_dir_name = args.directory + OUTPUT_DIRECTORY + "/" + out_dir + "GOs/"
                if not os.path.exists(new_dir_name):
                    os.mkdir(new_dir_name)
            print(output_str)

            out_file_name = args.directory  + OUTPUT_DIRECTORY + "/" + out_dir + "GOs/" + out_file + "GOs.out"
            print(out_file_name)
            outfile = open(out_file_name, "w")
            outfile.write(output_str)
            outfile.close()

#print("--- %s seconds ---" % (time.time() - start_time))
