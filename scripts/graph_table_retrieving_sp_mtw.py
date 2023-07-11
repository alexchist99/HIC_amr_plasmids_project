import sys
import pandas as pd


#define input and output files
infile1 = open(sys.argv[1])  
infile2 = open(sys.argv[2])
infile3 = open(sys.argv[3])

mag_info_file1 = pd.read_table(sys.argv[4],header=None)
mag_info_file2 = pd.read_table(sys.argv[5],header=None)
mag_info_file = mag_info_file1.merge(mag_info_file2,on=0)

threhold_file = open(sys.argv[6])
outfile = open(sys.argv[7], "w")
out_bins = open(sys.argv[8], "w")


#print(mag_info_file.head())


#for no inter MAG contatcts thresholds is 0.2


thr = float(threhold_file.readlines()[-1])

print(thr)



#check mag composition
cont_mag = dict()
for line in infile1:
    cont_mag[line.strip().split(":")[1]] = int(line.strip().split(":")[0].strip("bin."))
mag = set(cont_mag.keys())



#check palsmid composition
plas_names = set()
plas_dict = dict()
plasmid_dir = dict()
for line in infile3:
    words = line.strip().split()
    cont_name = words[0]
    plas_dict[cont_name] = list([words[1], words[2], words[3:]])
    plas_names.add(cont_name)



#check which contig to with has the contact
link_bin = dict()
n = 0
for line in infile2:
    words = line.strip().split(",")
    if words[1].strip('"') in mag:

        if words[2].strip('"') in plas_names:
            if words[2].strip('"') not in link_bin.keys():
                link_bin[words[2].strip('"')] = dict()
                link_bin[words[2].strip('"')][cont_mag[words[1].strip('"')]] = float(words[0])
            elif cont_mag[words[1].strip('"')] not in link_bin[words[2].strip('"')].keys():
                link_bin[words[2].strip('"')][cont_mag[words[1].strip('"')]] = float(words[0])
            else:
                link_bin[words[2].strip('"')][cont_mag[words[1].strip('"')]] = max(float(words[0]), link_bin[words[2].strip('"')][cont_mag[words[1].strip('"')]])

    if words[2].strip('"') in mag:
        if words[1].strip('"') in plas_names:
            if words[1].strip('"') not in link_bin.keys():
                link_bin[words[1].strip('"')] = dict()
                link_bin[words[1].strip('"')][cont_mag[words[2].strip('"')]] = float(words[0])
            elif cont_mag[words[2].strip('"')] not in link_bin[words[1].strip('"')].keys():
                link_bin[words[1].strip('"')][cont_mag[words[2].strip('"')]] = float(words[0])
            else:
                link_bin[words[1].strip('"')][cont_mag[words[2].strip('"')]] = max(float(words[0]), link_bin[words[1].strip('"')][cont_mag[words[2].strip('"')]])


#get high-quality mags and add their clasiification
bin_info = dict()
good_bin = set()
for index,line in mag_info_file.iterrows():
    words = list(line)
    if float(words[1]) >= 80 and float(words[2]) <= 5:
        #print(words)
        tax = words[3].split(";")
        family = tax[4]
        gener = tax[5]
        num_bin = int(words[0].strip("bin."))
        bin_info[num_bin] = (float(words[1]), float(words[2]), family, gener)
        good_bin.add(num_bin)



#get the threshold for inter and intra mags 

num_bin = 0
num_plas_connect = 0
num_plas_many = 0
plas_f_mag_f = dict()


left_bins = set()
for key in link_bin.keys():
    num_bin = 0
    for bin in link_bin[key].keys():
        if link_bin[key][bin] >= thr and bin in good_bin:
            if cont_mag.get(key,"None")=="None" or cont_mag.get(key,"None") in good_bin:
                left_bins.add(bin)
                num_bin += 1
                print(key,bin,
                  link_bin[key][bin],
                  cont_mag.get(key,"None") == bin,
                  plas_dict[key][0],
                  plas_dict[key][1],
                  *plas_dict[key][2],
                  *bin_info[bin],
                  file=outfile, sep="\t")
                
for bin in sorted(list(left_bins)):
    out_bins.write(f"bin.{bin}\n")

outfile.close()

