"""
Created on 04/11/21 

@author: Maxim Maron
"""
import os
import json
import csv
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-i","--inputFile",help="File containing a column with Accession and another column with amino acid position")
parser.add_argument("-o","--fileName",help="Name of output file")
args = parser.parse_args()

#set your working director
#os.chdir('/mnt/c/Users/maxim/Dropbox (EinsteinMed)/CloudStation/PRMT/Proteomics/')

#open and read the JSON file from mobiDB
with open('disorder_UP000005640.mjson', 'r') as f:
    lines = f.readlines()

#remove both leading and trailing charachters to clean up the format
data = [json.loads(l.strip()) for l in lines]

#you can write the file out to a new json, if desired
#json.dump(data, open("disorder_proteome.json", "w"))

#load in the json file to continue the analysis
#with open('disorder_proteome.json', 'r') as f:
#    data = json.load(f)

# Convert the list of dicts to a dictionary of lists   
data_dict = { 
    k: [l[k] for l in data] 
    for k in ['sequence', 'acc', 'mobidb_consensus']
    }

#create empty dictionary and load in the accessions as the keys and the list of disordered regions as values
disorder_dict = {}
for i in range(len(data_dict['acc'])):
    value = len(data_dict['mobidb_consensus'][i]['disorder']['predictors'])-1
    disorder_dict[data_dict['acc'][i]] = data_dict['mobidb_consensus'][i]['disorder']['predictors'][value]['regions']

#Loop through dictionary and remove the string 'D' from the list of disordered regions
for v1 in disorder_dict.values():
    for v2 in v1:
        v2.remove('D')

#load in your csv that contains a column with an accession and an integer for the amino acid position in question
data_file = open(args.inputFile)
reader=csv.DictReader(data_file)
df_dict = {}
#loop through the csv and create a dictionary with the accession as the key and a list of the amino acids as the values
for row in reader:
    if row['Accession'] in df_dict:
        df_dict[row['Accession']].append(row['position'])
    else:
        df_dict[row['Accession']] = [row['position']]
data_file.close()

#open an output file for writing
outF = open(args.fileName, "w")
#write the column headers
outF.write('\t'.join(['Accession', 'position', 'result']) +'\n')
#loop through the dictionary of accessions and amino acids in question and compare them to the mobiDB disordered amino acid ranges
for key in df_dict.keys():
    if key in disorder_dict.keys():
        for value in df_dict[key]:
            for i in range(len(disorder_dict[key])):
                if int(disorder_dict[key][i][0]) <= int(value) <=  int(disorder_dict[key][i][1]):
                    outF.write('\t'.join([str(key), str(value), 'TRUE']) +'\n')
outF.close()
