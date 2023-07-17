import glob
import os

path="/home/alexclear/HIC_amr_plasmids_project/data/raw_data"

files= [os.path.basename(f) for f in glob.glob(f'{path}/*')]

#WARNING: Change paths to relevant


list_files = sorted(list(set([i.split("_")[0] for i in set(files)])))
#print(list_files)


rule all3:
    input: expand("/home/alexclear/HIC_amr_plasmids_project/find_amr/output/{id}_amr.csv", id=list_files)


rule get_plasmid_sequences:
   input: 
       pl_id = "/home/alexclear/HIC_amr_plasmids_project/find_plasmids/output/samples_output/{id}_plasmids.txt",
       assembly = "../data/spades_contigs_out/{id}_contigs.fasta"
   output: "/home/alexclear/HIC_amr_plasmids_project/find_amr/input/{id}_plasmids.faa"
   conda: "exp"
   shell: '''
   cat {input.pl_id}|cut -f1 > ../find_amr/input/{wildcards.id}_pl.txt
   seqtk subseq {input.assembly} {input.pl_id} > {output}
   rm ../find_amr/input/{wildcards.id}_pl.txt
   '''
   

rule get_amr:
   input: 
       pl_seq = rules.get_plasmid_sequences.output
   output: "/home/alexclear/HIC_amr_plasmids_project/find_amr/output/{id}_amr.csv"
   conda: "rgi"
   shell: '''
      cd ../tools/rgi
      rgi main --input_sequence {input.pl_seq} --output_file /home/alexclear/HIC_amr_plasmids_project/find_amr/output/{wildcards.id}_tmp --input_type contig --local --clean -n 40 
      cd ../../try
      cat ../find_amr/output/{wildcards.id}_tmp.txt |cut -f2,6,9,15,16,17|sed "s/\t/,/g"> {output}
      rm ../find_amr/output/{wildcards.id}_tmp.txt
      rm ../find_amr/output/{wildcards.id}_tmp.json
      '''
