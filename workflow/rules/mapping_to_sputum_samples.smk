import glob
import os

path="/home/alexclear/HIC_amr_plasmids_project/data/raw_data"

files= [os.path.basename(f) for f in glob.glob(f'{path}/*')]

list_files = sorted(list(set([i.split("_")[0] for i in set(files)])))

#get samples corresponding to sputum isolates samples
list_files = [list_files[1],list_files[-2]]


rule all4:
   input: expand("/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/{id}_classif_intersected.txt", id=list_files)

"/home/alexclear/HIC_amr_plasmids_project/data/raw_data_wgs"
rule get_sd_plasmid_sequences:
   input: 
      spades = "/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/{id}/spades_out/contigs.fasta",
      nodes = "/home/alexclear/HIC_amr_plasmids_project/find_plasmids/output/samples_output/{id}_plasmids.txt"
   output:"/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/plasmids_seq/plasmids_{id}.fa"
   conda: "exp"
   shell: '''
   cat {input.nodes} |cut -f1|sort|uniq > /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/temp_{wildcards.id}.txt
   seqtk subseq {input.spades} /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/temp_{wildcards.id}.txt > {output}
   rm /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/temp_{wildcards.id}.txt
   '''

rule mapping:
  input: 
      reference = rules.get_sd_plasmid_sequences.output,
      read_sputum_r1 = "/home/alexclear/HIC_amr_plasmids_project/data/raw_data_wgs/{id}_B_R1.fastq.gz",
      read_sputum_r2 = "/home/alexclear/HIC_amr_plasmids_project/data/raw_data_wgs/{id}_B_R2.fastq.gz"
  output: 
      bam =  "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/sample_{id}.bam",
      stats = "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/log/sample_{id}.stats.txt"
  conda: "exp"
  shell:'''
  #renamed B-9782 into SD002_B; B-9817 into SD008_B as they match to each other

  #index:
  bwa index {input.reference}
 
  # map paired-end sample
  bwa mem -t 4 {input.reference} {input.read_sputum_r1} {input.read_sputum_r1} > /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/sample_{wildcards.id}.sam
  

  # convert sam to sorted bam, index and stats
	samtools view -bS /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/sample_{wildcards.id}.sam|samtools sort -o /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/sample_{wildcards.id}.bam
	samtools index {output.bam}
	samtools stats {output.bam} | grep ^SN > {output.stats}
  rm /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/sample_{wildcards.id}.sam
  '''

rule get_mapped_only:
  input:"/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/sample_{id}.bam"
  output: "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/mapped_{id}.bam"
  conda:"exp"
  shell:'''
  # get mapped reads which mate is always mapped 
samtools view -u -F12 {input} > temp_{wildcards.id}.bam
# name-sort bam file because we need read pairs to go one by one 
samtools sort temp_{wildcards.id}.bam -o {output} 
rm temp_{wildcards.id}.bam
  '''

rule get_high_covered_contigs:
   input: 
         mapped = rules.get_mapped_only.output,
         all_cont = rules.get_sd_plasmid_sequences.input.nodes
   output: 
        intersected_contigs = '/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/{id}_intersected_nodes_list.txt',
        all_contigs = '/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/{id}_all_nodes_list.txt'
   conda:"exp"
   shell:'''
   #covergae 90(%) is affordable for getting reliable result (checked via samtools covergae --hist):
   samtools coverage {input.mapped}|cut -f1,6|awk '{{if ($2>=90){{print $1}}}}'|sort|uniq|grep "NODE" > {output.intersected_contigs}

   #get all nodes from one of the previous step:
   cat {input.all_cont}|cut -f1|sort|uniq> {output.all_contigs}
   '''


rule get_amr_contigs:
      input: 
          target = rules.get_high_covered_contigs.output.intersected_contigs,
          amr_info = "/home/alexclear/HIC_amr_plasmids_project/find_amr/output/{id}_amr.csv"
      output:  
          amr_of_intersected_contigs = '/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/{id}_amr_intersected.txt'
      conda: "exp"
      shell:'''

      #grep -f /home/alexclear/HIC_amr_plasmids_project/sputum_analysis/SD002_intersected_nodes_list.txt /home/alexclear/HIC_amr_plasmids_project/make_graphs/out/SD002/graph/SD002_graph_table.txt|cut -f1,11
      
      #get intersected nodes with amr:
      grep -f {input.target} {input.amr_info} > {output.amr_of_intersected_contigs}
      '''
rule get_classification:
   input: 
       mapped_inters = rules.get_high_covered_contigs.output.intersected_contigs,
       all_cont = rules.get_sd_plasmid_sequences.input.nodes
   output: "/home/alexclear/HIC_amr_plasmids_project/sputum_analysis/{id}_classif_intersected.txt"
   shell: '''
   grep -f {input.mapped_inters} {input.all_cont}|cut -f1,8 > {output}
   '''

