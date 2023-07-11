from itertools import permutations,combinations_with_replacement
import pandas as pd
from collections import defaultdict
import numpy as np
import sys
import os
import re



def connections_retrieving(fl_name,sample,dir_way,thr):

   dt = pd.read_table(fl_name,header = None)
   #print(dt[2])

   #use info about contig source (-5 column)
   dt_new = dt.iloc[:,[0,-5,-1]].sort_values(by = [0])
   #print("new",dt_new)

   dct = defaultdict(list)
   dct_to_compare = defaultdict(list)
   for n in range(dt_new.shape[0]):
       l = list(dt_new.iloc[n,[0,1,2]])
       dct_to_compare[l[0]].append(l[1])
       dct[l[0]].append(l[2])
   #print(dct)
   
   network =  defaultdict(list)
   for n in dct:
      #check if the plasmid contig source is between at least one fitted host (eg.: Klebsiella`s contig - Klebsiella MAG)
      complement = any("g__"+taxa.split("_")[0] in dct[n] for taxa in dct_to_compare[n])
      print(dct_to_compare[n])
      print(complement)

      if len(dct[n])>1:
          var = list(permutations(set(dct[n]), 2))
      else:
          var = list(set(list(combinations_with_replacement(set(dct[n]), 2))))

      for k in var:
            print(k,n)
            network["node"].append(n)
            network["from"].append(k[0])
            network["to"].append(k[-1])
            if complement:
               network["compliment"].append("yes")
            else:
               network["compliment"].append("no")
          
      
       #var = list(combinations_with_replacement(set(dct[n]), 2))

       #print(list(combinations_with_replacement(set(dct[n]), 2)), dct[n],n)
       #print(var, dct[n],n ) 
       #var = list(set(list(combinations_with_replacement(set(dct[n]), 2))+ list(permutations(set(dct[n]), 2))))
       #print(list(permutations(set(dct[n]), 2)))

      
   network_df = pd.DataFrame(network)
   network_df["sample"] = sample
   network_df.to_csv(f"{dir_way}/{sample}_data.csv",index=None)
   return network_df



content = [feature.strip() for feature in open(sys.argv[1]).readlines()]
samples = content[::2]
thresholds = [float(thr) for thr in content[1::2]]
dict = dict(zip(samples, thresholds))     
#print(d)


for r, d, f in os.walk(sys.argv[2]):
    for file in f:
        if "graph" in file:
            sample = file[:file.find("_")]
            thr = dict[sample]
            connections_retrieving(f"{sys.argv[2]}/{file}",sample,sys.argv[2],thr)



def concat(dir_way):
  all_files = list()
  def merging(ar,data):
    #print(ar,data)
    ar = pd.read_csv(ar)
    ar= ar[["Contig","Best_Hit_ARO","Drug Class"]]
    data = pd.read_csv(data)

    l = [re.search("[^_]*_[^_]*_[^_]*_[^_]*_[^_]*_[^_]*", i).group(0) for i in ar["Contig"]]
    ar["Contig"] = l
   
    merged_data = data.merge(ar,how='left',left_on = "node",right_on = "Contig")
    return merged_data

  for r, d, f in os.walk(dir_way):
    for file in f:
        if "data" in file:
          sample = file[:file.find("_")]
          ar = f'../find_amr/output/{sample}_amr.csv'
          data = f'{dir_way}/{sample}_data.csv'

          mrg = merging(ar,data)
          all_files.append(mrg)
          #print(pd.concat(all_files))
  final_file = pd.concat(all_files)
  final_file.to_csv(f"{sys.argv[3]}",index=None)


  return final_file


concat(sys.argv[2]) 
          

