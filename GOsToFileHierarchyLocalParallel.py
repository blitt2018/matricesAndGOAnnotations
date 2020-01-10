#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:00:09 2019

@author: benlitterer

This file takes blast output, retrieves the corresponding GO ids for the results, and then maps these results to a file structure. 
If the output directory specified exists, files are placed there. If not, the output directory specified is created. Processing is done with parallel threads.
"""


import glob
import os
from multiprocessing import Pool 
import time
import argparse
import pandas as pd

start_time = time.time()
parser = argparse.ArgumentParser(description="directory: path to directory containing your BLAST output. output_dir: name of directory to be head of output file tree")
parser.add_argument("directory")
parser.add_argument("output_dir") 

args = parser.parse_args()

""" ADJUST THESE FOR YOUR PIPELINE"""

EVAL_CUTOFF =  .01                                                     #keep evals smaller or equal to this number 

EVAL_COLUMN_NUMBER = 3                                                 #which column your evals are in the first column being 1 rather than zero 

#this line is redundant but exists because this code was modified based off of already working code 
OUTPUT_DIRECTORY = args.output_dir

# path to file containing GO ids for every uniprot protein id
GO_REF_FILENAME = "/home/benlitterer/Academic/Research/ProjectMatrices/referance/uniprotKB_unfiltered_human.gaf"

#create directory to put mapped protein db in 
new_dir_name = args.directory + OUTPUT_DIRECTORY
if not os.path.exists(new_dir_name):
    os.mkdir(new_dir_name)

#db_hits_list=[set(['a','b','z','q']), set(['a','b','c','f']), set(['a','c']), set(['a','c','f',]), set(['a','d'])]
db_hit_list=[]
#iterate over all .out files in current working directory
#print("here are the files that contain output and the output length")
file_list = []

"""
This for loop goes through each output file in the blast output directory and creates a list of dictionaries called db_hit_list. One particular dictionary represents output for a particular blast version. 
Each dictionary contains key value pairs where the key is the query given to blast, and the value is a set (unique items only) of hits given from that version of blast. 

[{query1: {hit1, hit2, hit3}},{query1: {hit1, hit4, hit5}},{query1: {hit6, hit7, hit8}}]

To track which dictionary goes with which version of blast (and by extension which blast matrix), a list called file_list is created where the order corresponds with the order 
of the dictionaries in db_hit_list
"""
print("files read: ")
for filename in glob.glob(args.directory + '/*.out'):
    #print("\n \n" + filename + "\n \n")
    hit_dict = {}
    #don't open slurm output files, and don't open empty files 
    if 'slurm' not in filename and os.stat(filename).st_size > 0:
        print(filename)
        """
        FIRST VERSION: this is 4x slower but loads into a pandas dataframe for easier data exploration
        print("FILE: " + filename)
        hit_set = set(pd.read_csv(filename, '\t', header=None).loc[:, 1])
        """
        file_list.append(filename)
        with open(filename) as file: 
            content = file.readlines()
        """
        IMPORTANT NOTE
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
                
            """
            IMPORTANT: this if statement keeps hits that don't satisfy eval requirement from being considered. This means that no GO ids will be generated for these hits. 
            """
            if curr_eval <= EVAL_CUTOFF: 
                #print("\t HIT: " + curr_hit + " EVAL: " + str(curr_eval))
                hit_list.append(curr_hit)
            prev_qname = curr_qname
        db_hit_list.append(hit_dict)
#a couple random tests 
#print(file_list)
print("the 6th output file in the dict is " + str(file_list[6]))

#enter queries that you used when blasting originally
#print("this file should contain the following queries: YP_003024026.1, NP_001335218.1, XP_024308266.1" )
#print("this statement has been programatically evaluated to be: " + str(db_hit_list[6]["YP_003024026.1"] != None and db_hit_list[6]["NP_001335218.1"] != None and db_hit_list[6]["XP_024308266.1"] != None))

"""
don't think this is necessary 
#generate all values for each query result 
all_vals_dict = {}
for hit_dict in db_hit_list:
    for key in hit_dict.keys():
        all_vals_dict[key] = hit_dict[key]
"""


cols = ["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence Code", "With (or) From", "Aspect", "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type", "Taxon and Interacting taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"]
relevant_cols = ["DB", "DB_Object_ID", "GO_ID"]

#read in the data linking uniprot ids to go ids. cut down the columns to only the ones we need. 
GO_db = pd.read_csv(GO_REF_FILENAME, "\t", names = cols, usecols = relevant_cols)

def getGOsDBMethod(GO_id): 
    return list(set(GO_db[GO_db["DB_Object_ID"] == GO_id]["GO_ID"]))


for i in range(len(db_hit_list)):
    hit_dict = db_hit_list[i]
    associated_file = file_list[i]
    dont_include_hit_list = []

    #create list of all dictionaries that doesn't include the one currently being iterated over
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
                    pool_out_list = p.map(getGOsDBMethod, pool_args)
                    p.close()
                    #p.join()
                    #print("call being done on: " + str(pool_args))
                    for j in range(len(pool_out_list)): 
                        GOs = pool_out_list[j]
                        hit_str = pool_args[j]
                        print(GOs)
                        #always want to get second line because the first is simply the column labels 
                        if GOs:
                            for GO in GOs: 
                                if GO in go_key_dict.keys(): 
                                    go_key_dict[GO].append(hit_str)
                                else: 
                                    go_key_dict[GO] = [hit_str]
                    pool_args=[curr_str]
                else: 
                    pool_args.append(curr_str)
            print(go_key_dict)
            output_str = ""
            
            for key, val in go_key_dict.items(): 
                output_str += key + "\n"
                output_str += str(val) + "\n" 
                
                new_dir_name = args.directory + OUTPUT_DIRECTORY + "/" + out_dir + "GOs/"
                if not os.path.exists(new_dir_name):
                    os.mkdir(new_dir_name)
            
            out_file_name = args.directory  + OUTPUT_DIRECTORY + "/" + out_dir + "GOs/" + out_file + "GOs.out"
            print(out_file_name)
            outfile = open(out_file_name, "w")
            outfile.write(output_str)
            outfile.close()
            
#print("--- %s seconds ---" % (time.time() - start_time))
