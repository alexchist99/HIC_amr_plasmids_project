import pandas as pd
import numpy as np
import sys 

def without_assembly_info(plasflow,PlDB_all,VirV,MobR,MobT,thr):

    MobT = MobT.loc[MobT["rep_type(s)"]!="-"]
    MobT["sample_id"] = [i.split()[0] for i in MobT["sample_id"]]
    MobT = MobT [["sample_id","rep_type(s)","predicted_mobility"]]
    MobT = MobT.rename(columns = {"sample_id":"contig_name"})
    

    MobR["NODE"] = ["".join(i.split(":")[1:]) for i in MobR["contig_id"]]
    MobR["molecule_type"] = ["plasmid" if i!="-" else "-" for i in MobR["mash_neighbor_identification"]]
    MobR = MobR.rename(columns = {'contig_id':'contig_name'})
    MobR = MobR[["contig_name","molecule_type","NODE","mash_neighbor_identification"]]
    MobR = MobR.rename(columns = {'molecule_type':'Mob_recon'})

    
    PlDB_all = PlDB_all.loc[(PlDB_all[5]>=80) & (PlDB_all[6]>=60)]
    PlDB_all=PlDB_all.drop_duplicates(subset=[0],inplace=False)
    PlDB_all = PlDB_all.rename(columns = {0:'contig_name'})
    PlDB_all["NCBI_id0.8cov0.6"]=PlDB_all[8]
    PlDB_all = PlDB_all[["contig_name","NCBI_id0.8cov0.6"]]
   
    
    plasflow_pl = plasflow.loc[plasflow.label.str.contains("plasmid")][["contig_name","contig_length","label"]]
    plasflow_pl = plasflow_pl[["contig_name","label"]]
    plasflow_pl = plasflow_pl.rename(columns = {'label':'Plasflow'})

    
    Vir=VirV.rename(columns = {'Contig name':'contig_name'})
    Vir = Vir[["contig_name","Prediction","Length","Circular"]]
    Vir = Vir.rename(columns = {'Prediction':'ViralVerify'})

    
    
    a= pd.merge(pd.merge(pd.merge(PlDB_all,Vir,how='outer')
                         ,MobR,how='outer'),plasflow_pl,how='outer')
    merged_dt = pd.merge(a,MobT,how='outer')
    
    #merged_dt = a[["contig_name","NODE","NCBI_id0.8cov0.6","Plasflow","Mob_recon","ViralVerify","Circular","Length"]]
    merged_dt = merged_dt[["contig_name","NODE","NCBI_id0.8cov0.6","Plasflow","Mob_recon","ViralVerify","Circular","Length","mash_neighbor_identification","rep_type(s)","predicted_mobility"]]
    merged_dt["rep_type(s)"] = merged_dt["rep_type(s)"].fillna("Nope")

    nn = []

    for n in range(merged_dt.shape[0]):
        l = list(merged_dt.iloc[n,[1,2,3,4,5]])
        #print(l)
        pl = " ".join([str(k) for k in l]).count("plasmid")
        Pl = " ".join([str(k) for k in l]).count("Plasmid")
        pl_gen = Pl+pl
        #or list(merged_dt.iloc[n,[-2]])[0]!="Nope"
        if pl_gen > int(thr):
            nn.append(n)
    
    filtered_dt = merged_dt.iloc[nn]
    filtered_dt = filtered_dt.fillna("Nope")
   
    
    pl= []
    for i in range(filtered_dt.shape[0]):
        n1,n2 = filtered_dt.iloc[i,[0,6]]
        #print(n1,n2)
        if "+" in n2:
            pl.append("Plasmid_full_assembly")
        else:
            pl.append("Plamid_contig")
    filtered_dt["Plasmid_status"] = pl

    filtered_dt["Plasmid_annotation"] = ["_".join(k.split()[1:3]) for k in filtered_dt["NCBI_id0.8cov0.6"]]
    filtered_dt["ViralVerify"] = ["_".join(k.split()) for k in filtered_dt["ViralVerify"]] 

    #filtered_dt = filtered_dt[["contig_name","NODE","Plasflow","Mob_recon","ViralVerify","Length","Plasmid_status","Name"]]
    filtered_dt_final = filtered_dt.drop(columns = ["Plasflow","Circular","NCBI_id0.8cov0.6","mash_neighbor_identification","NODE"])    
    return filtered_dt_final
    
    
plasflow = pd.read_table(sys.argv[1])
PlDB_all = pd.read_table(sys.argv[2],header=None)
VirV= pd.read_csv(sys.argv[3],sep =",")
MobR = pd.read_csv(sys.argv[4],sep="\t")
MobT = pd.read_table(sys.argv[5],sep="\t")
thr = sys.argv[6]
res = without_assembly_info(plasflow,PlDB_all,VirV,MobR,MobT,thr)
res.to_csv(sys.argv[7], header=None, index=None, sep='\t', mode='w')
    

