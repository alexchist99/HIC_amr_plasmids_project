

###WARNING!
#1) DATABASES YOU SHOULD DOWNLOAD SEPARETELY VIA THE FOLLOWING LINKS: 
#NCBI NUCLEOTIDE DATABASE OF PALSMIDS: https://ccb-microbe.cs.uni-saarland.de/plsdb INTO HIC_amr_plasmids_project/databases/PLSDB/

#PFAM DATABASE: http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz INTO HIC_amr_plasmids_project/databases/Pfam/
#2) YOU CAN KEEP ALL OF OUTPUT FILES FROM EACH TOOL JUST BY REMOVING OPTIONAL CODE IN EACH RULE 
#3) SOME TOOLS (E.G MOBTYPER AND VIRALVERIFY) WORK VERY SLOW EVEN WITH 40 THREADS BY DEFAULT. SO,I HIGHLY RECOMEND YOU TO RUN IT ON 
#"SCREEN" OR "TMUX" AND GRAB SOME TEA OR SOMETHING STRONGER


import glob
import os


path="/home/alexclear/HIC_amr_plasmids_project/data/raw_data"

files= [os.path.basename(f) for f in glob.glob(f'{path}/*')]


list_files = sorted(list(set([i.split("_")[0] for i in set(files)])))

#check working directory and samples
#print(f"working directory: {os. getcwd()}")
print(f"samples: {list_files}")


rule all1:
    input: expand("../find_plasmids/output/samples_output/{id}_plasmids.txt", id=list_files)


rule Plasflow:
    input: "/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/{id}/spades_out/contigs.fasta"
    output: "../find_plasmids/output/{id}/plsflow_out/plasflow_predictions_{id}.tsv"
    conda:"plasflow"
    shell: '''
    python3 ../tools/PlasFlow/PlasFlow.py --input {input} --output {output} --threshold 0.95
    '''

rule Mob_T:
    input: "/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/{id}/spades_out/contigs.fasta"
    output: "../find_plasmids/output/{id}/mobt_out/MOB_TYP_{id}.txt"
    conda:"mobsuite"
    shell: '''
    mob_typer --multi --infile {input} --out_file ../find_plasmids/output/{wildcards.id}/mobt_out/mobtyper_results.txt -n 40 
    cat ../find_plasmids/output/{wildcards.id}/mobt_out/mobtyper_results.txt |grep -v "non-mobilizable"|cut -f14,1,17,6|awk '{{print $4" "$5"\t"$1"\t"$2"\t"$3}}'> {output}
    '''


rule Mob_R:
    input: "/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/{id}/spades_out/contigs.fasta"
    output: "../find_plasmids/output/{id}/mobr_out/MOB_REC_{id}.txt"
    conda:"mobsuite"
    shell: '''
    mob_recon --infile {input} --outdir ../find_plasmids/output/{wildcards.id}/mobr_out -c -f -n 40
    cat ../find_plasmids/output/{wildcards.id}/mobr_out/contig_report.txt|cut -f5,21 > {output}
    '''

rule Virv:
    input: "/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/{id}/spades_out/contigs.fasta"
    output: "../find_plasmids/output/{id}/viral_v/viralverify_{id}.csv"
    conda:"viralverify"
    shell: '''
    viralverify -f {input} --hmm ../databases/Pfam/Pfam-A.hmm -o ../find_plasmids/output/{wildcards.id}/viral_v -t 40
    cat ../find_plasmids/output/{wildcards.id}/viral_v/contigs_result_table.csv|awk -F"," '{{print $1","$2","$3","$4}}'> {output}
    '''


rule NCBI:
    input: "/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/{id}/spades_out/contigs.fasta"
    output: "../find_plasmids/output/{id}/NCBI/PLSDB_{id}.txt"
    conda:"hicmag_py37"
    shell: '''
   blastn -query {input} -db ../databases/PLSDB/plsdb.fna -out {output} -outfmt "6 qseqid  sseqid qlen slen evalue pident qcovs length salltitles" -num_threads 40 -evalue 10E-5
'''

rule concatenate:
   #`params:2` means that we consider contig as plasmid if more or equal then two tools give plasmid identification output  
   input:
      plasflow = rules.Plasflow.output,
      mob_typer = rules.Mob_T.output,
      mob_recon = rules.Mob_R.output,
      viralverify =  rules.Virv.output,
      #viralverify = "/home/alexclear/HIC_amr_plasmids_project/find_plasmids/output/SD002/viral_v/viralverify_SD002.csv",
      plsdb = rules.NCBI.output
      #plsdb = "/home/alexclear/HIC_amr_plasmids_project/find_plasmids/output/SD002/NCBI/PLSDB_SD002.txt",
      #mob_typer = "/home/alexclear/HIC_amr_plasmids_project/find_plasmids/output/SD002/mobt_out/MOB_TYP_SD002.txt"
   output: "../find_plasmids/output/samples_output/{id}_plasmids.txt"
   conda:"exp"
   params: 2
   shell: '''
   python3 ../scripts/concatenation2_mtw.py {input.plasflow} {input.plsdb} {input.viralverify} {input.mob_recon} {input.mob_typer} {params} {output}
       
   #optional,to save space (delete fasta files):
   rm -rf ../find_plasmids/output/{wildcards.id}/mobt_out/*.fasta
   rm -rf ../find_plasmids/output/{wildcards.id}/mobr_out/*.fasta
   rm ../find_plasmids/output/{wildcards.id}/plsflow_out/*.fasta
 
  
   rm -rf ../find_plasmids/output/{wildcards.id}/viral_v/*.fa
   rm -rf ../find_plasmids/output/{wildcards.id}/viral_v/*.fasta
   rm -rf ../find_plasmids/output/{wildcards.id}/viral_v/Prediction_results_fasta
   
   '''
    