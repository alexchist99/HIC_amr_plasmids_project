
import numpy as np
import pandas as pd
import sys
from scipy import stats
import matplotlib.pyplot as plt
from collections import defaultdict as dfd

def draw(plt, f1,plasmids,h,dnst):

    #add function to get dataframe with link density of inter/intra chromosome density, based on metawrap and viralverify:
    
    def get_dataframe(contigs,plasmids):
        c = False
        plasmid = []
        d = dfd(list)
        for line in open(plasmids):
            #print(line.split())
            node= line.split()[0].strip('"')
            Type = line.split()[4]
            which_mag =line.split()[3]
            mag_name = line.split()[-1].split("__")[1]
            weihgt = float(line.split()[2])
            node_in_bin = contigs.get(node, -1)
            #print(node_in_bin)
            #print(line.split()[0].strip('"'))


            if node_in_bin!=-1:
                d["Pl_w"].append(np.log(weihgt))
                if Type=="Chromosome":
                    if int(which_mag=="True")!=0:
                            d["Type"].append("intra_chromosome")
                    else: 
                            c=True
                            d["Type"].append("inter_chromosome")
                else:
                     d["Type"].append(Type)
                d["Mag"].append(mag_name)

        data_bondary = pd.DataFrame(d)
        return data_bondary,c
    
    # contig -> MAG for contigs in high-quality MAGs
    contigs = dict() 
    for line in open(f1):
        contigs[line.strip().split(":")[1]] = int(line.strip().split(":")[0].strip("bin."))

    data_bondary,is_inter = get_dataframe(contigs,plasmids)
    
    #if we have more then one MAG in metawrap then we have inter chromosome contacts and can define threshold, 
    # otherwise just draw a plot 

    if is_inter:

        inter = data_bondary.query('Type=="inter_chromosome"')["Pl_w"]
        intra = data_bondary.query('Type=="intra_chromosome"')["Pl_w"]

        #get the threshold

        def p_val(inr,outr,name):
            thr_list = []
            inr, outr = np.round(np.exp(inr),2), np.round(np.exp(outr),2)
            for k in np.arange(0.2,int(max(outr)),0.1):
                p = np.sum(outr >= k) / len(outr)
                if p<= 0.05:
                    thr_list.append(k)
                    break
            for j in np.arange(0.2,int(max(inr)),0.1):
                p = np.sum(inr >= j) / len(inr)
                if p>= 0.05:
                    thr_list.append(j)
                    break
            mean_thr = np.mean(thr_list)

                    #print(f"p-value = {p:.4f},sample = {name},thr = {k:.4f}")
            ks = stats.kstest(inr,outr)

            with open(sys.argv[5],"w") as fl:
                fl.write(f"{name}\n{mean_thr:.4f}\n")
                #fl.write(f"sample = {name}\nintra_thr = {k:.4f}\ninter_thr = {j:.4f}\nmean_thr = {mean_thr:.4f}\nks_statistic = {ks[0]:.4f}\nks_pvalue = {ks[1]:.4f}\n")   
            return mean_thr
        thr = p_val(intra,inter,sys.argv[3])
        
            
        # for pl_line in open(pl_file):
        #     plasmid_w = np.log(float(pl_line.split()[2]))
        #     plasmid.append(plasmid_w)
            
        
        plt.hist(list(inter), color='red', density=dnst, bins=100, label="Inter_chromosome", alpha=0.5)
        plt.hist(list(intra), color='green', density=dnst, bins=100, label="Intra_chromosome", alpha=0.5)
        #plt.hist(plasmid, color='blue', density=dnst, bins=100, label="Plasmid", alpha=0.5)
        plt.legend()
        
         #      #set custom ticks
        lbls = [0.1, 0.2, 0.3, 0.6, 1, 2, 5, 10, 30]
        #0.470
        
        plt.axvline(np.log(thr), color='k', linestyle='dashed', linewidth=1)
        plt.set_xticks(np.log(lbls))
        plt.set_xticklabels(lbls)
        #plt.grid(axis='x', color='k', linestyle='--')
        plt.set_title(h,fontsize = 10)
        plt.set_xlim(np.log(0.1), np.log(30))
        #plt.set_ylabel("Link count")
        plt.set_ylabel("Relative reachness")
        plt.set_xlabel("Normalized HI-C contact intensity")
        plt.figure.savefig(sys.argv[4])

    else:
        intra = data_bondary.query('Type=="intra_chromosome"')["Pl_w"]

        plt.hist(list(intra), color='green', density=dnst, bins=100, label="Intra_chromosome", alpha=0.5)
        plt.legend()
        lbls = [0.1, 0.2, 0.3, 0.6, 1, 2, 5, 10, 30]

        #plt.axvline(np.log(thr), color='k', linestyle='dashed', linewidth=1)
        plt.set_xticks(np.log(lbls))
        plt.set_xticklabels(lbls)
        #plt.grid(axis='x', color='k', linestyle='--')
        plt.set_title(h,fontsize = 10)
        plt.set_xlim(np.log(0.1), np.log(30))
        #plt.set_ylabel("Link count")
        plt.set_ylabel("Relative reachness")
        plt.set_xlabel("Normalized HI-C contact intensity")
        plt.figure.savefig(sys.argv[4])
        
        with open(sys.argv[5],"w") as fl1:
             fl1.write(f"{sys.argv[3]}\n{0.2}\n")
                #fl1.write(f"sample = {sys.argv[3]}\nintra_thr = None\ninter_thr = None\nmean_thr = None\nks_statistic = None\nks_pvalue = None\n")  

fig, axs = plt.subplots()
fig.tight_layout(pad=5.0)
draw(axs, sys.argv[1], sys.argv[2],sys.argv[3],True)
        


