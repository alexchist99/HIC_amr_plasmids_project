import glob
import os

path="/home/alexclear/HIC_amr_plasmids_project/data/raw_data"

files= [os.path.basename(f) for f in glob.glob(f'{path}/*')]




list_files = sorted(list(set([i.split("_")[0] for i in set(files)])))
list_files = [list_files[4],list_files[9]]
#list_files1 = [*list_files[:4],list_files[-1]]

print(list_files)

"../make_graphs/out/{id}/graph/{id}_graph_table.txt"
rule all2:
    input: expand("../make_graphs/out/{id}/graph/{id}_graph_table.txt", id=list_files)


rule spades_metawrap:
   input:  
            spades_contigs = "../data/spades_contigs_out/{id}_contigs.fasta",
            metawrap_list = "../data/metawrap/{id}_mtw"  
            
   output:  
            marked_spades = "../data/spades_contigs_out/{id}_replaced.fa",
            genomes = "../make_graphs/out/{id}/genomes/marked_mtw_mags.faa"
   conda:"exp"
   shell:'''
       bash ../scripts/spades_vs_metawrap_preprocess.sh {input.spades_contigs} {input.metawrap_list} {wildcards.id} {output.marked_spades} {output.genomes}
       less {output}|grep ">"|wc -l 
       less {input.spades_contigs}|grep ">"|wc -l 
   '''


# if there are indexes in ../data/spades_contigs_out/{id}_file then  bamm make -d return ERROR:
'''
ERROR:
You didn't specify that index files have been kept but there appears to be bwa index files present.
I'm cowardly refusing to run so as not to risk overwriting.
Force overwriting to create new indices
'''

rule hic_mapping:
  input: 
       fastq_r1 = '../data/raw_data/{id}_HiC_R1.fastq.gz',
       fastq_r2 = '../data/raw_data/{id}_HiC_R2.fastq.gz',
       spades_contigs = '../data/spades_contigs_out/{id}_replaced.fa'
  output: 
       link = '../make_graphs/out/{id}/bamm/links.tsv' 
  conda: "for_graphs"

  shell: '''
       bamm make -d {input.spades_contigs} -c {input.fastq_r1}  {input.fastq_r2} -o ../make_graphs/out/{wildcards.id}/bamm -t 65
       bamm parse -l {output.link} -b ../make_graphs/out/{wildcards.id}/bamm/*.bam -v
       '''

rule calculate_coverage:
  input: 
      #link_file = rules.hic_mapping.output.link,
      link_file  = "../make_graphs/out/{id}/bamm/links.tsv",
      depth_info = '../data/metabat/{id}_depth_metabat.tsv'
  output: 
       contig_info = '../make_graphs/out/{id}/contig_info.txt'
  conda: 'hicmag_py37'
  
  shell: '''
     jgi_summarize_bam_contig_depths --outputDepth {input.depth_info} ../make_graphs/out/{wildcards.id}/bamm/*.bam -v
     cut -f1,2,3 {input.depth_info} | tail -n +2 > {output.contig_info}
         '''

rule reprocess_counts:
  input:
      link_file = "../make_graphs/out/{id}/bamm/links.tsv"
      #rules.hic_mapping.output.link

  output: 
       links_count_strip = '../make_graphs/out/{id}/links_count_strip.txt',
       links_good_strip = '../make_graphs/out/{id}/links_good_strip.txt'
  params:"../data/metawrap/{id}_mtw_stats"
      
  shell: '''
     cut -f1,2 {input.link_file}| sort | uniq -c > ../make_graphs/out/{wildcards.id}/bamm/links_count.tsv
     
     awk 'substr($2,1,5)==substr($3,1,5)' ../make_graphs/out/{wildcards.id}/bamm/links_count.tsv > ../make_graphs/out/{wildcards.id}/bamm/links_intra.txt
     awk '$2>80 && $3<5' {params} | awk '{{print $1}}'|sort|uniq > ../make_graphs/out/{wildcards.id}/good_mags.txt
     
     #use metawrap results to get good MGSs

     grep -f ../make_graphs/out/{wildcards.id}/good_mags.txt ../make_graphs/out/{wildcards.id}/bamm/links_intra.txt|awk '{{print $1"\t"$2"\t"$3}}' | awk 'sub(/:bin.[0-9]*/,"",$2)'| awk 'sub(/:bin.[0-9]*/,"",$3)'> {output.links_good_strip}
     awk '{{$1=$1; print}}' ../make_graphs/out/{wildcards.id}/bamm/links_count.tsv | awk 'sub(/:bin.[0-9]*/,"",$2)'| awk 'sub(/:bin.[0-9]*/,"",$3)'> {output.links_count_strip}
         '''

rule hiczin_normalization:
  input: 
     links_count_strip = rules.reprocess_counts.output.links_count_strip,
     links_good_strip = rules.reprocess_counts.output.links_good_strip,
     contig_info = rules.calculate_coverage.output.contig_info

  output: 
       good_file = '../make_graphs/out/{id}/good.csv'
  conda:"exp"
  shell: '''
  Rscript ../scripts/hiczin.R -i {input.links_count_strip} -s {input.links_good_strip} -f {input.contig_info} -o {output.good_file} -m nb
  '''


rule find_threshold:
     input: 
          bins_info = rules.spades_metawrap.input.metawrap_list,
          good_file = "../make_graphs/out/{id}/good.csv",
          virv = "../find_plasmids/output/{id}/viral_v/viralverify_{id}.csv",
          checkm = "../data/metawrap/{id}_mtw_stats",
          gtdb = "../data/gtdbtk/{id}_gtdbtk_sorted.txt"
     output: 
           pic = "../make_graphs/out/{id}/{id}_hist.jpg",
           bondary = "../make_graphs/out/{id}/fileout.txt"
     conda:"exp"
     shell:'''
      cat {input.bins_info}|awk '{{print $2":"$1}}'> ../data/metawrap/{wildcards.id}_bins_info.txt
      cat {input.checkm}| cut -f1,2,3 > ../make_graphs/out/{wildcards.id}/checkm_modified.txt
      cat {input.virv}| tr  ',' '\t' > ../make_graphs/out/{wildcards.id}/temp_virv.txt
      printf "{wildcards.id}\n0.2\n" > ../make_graphs/out/{wildcards.id}/tmp_fileout.txt
      python3 ../scripts/graph_table_retrieving_sp_mtw.py ../data/metawrap/{wildcards.id}_bins_info.txt {input.good_file} ../make_graphs/out/{wildcards.id}/temp_virv.txt ../make_graphs/out/{wildcards.id}/checkm_modified.txt {input.gtdb} ../make_graphs/out/{wildcards.id}/tmp_fileout.txt ../make_graphs/out/{wildcards.id}/temp_virv_graph.txt ../make_graphs/out/{wildcards.id}/temp_bins.txt
      python3 ../scripts/plots_hic_intensity_spades.py ../data/metawrap/{wildcards.id}_bins_info.txt ../make_graphs/out/{wildcards.id}/temp_virv_graph.txt {wildcards.id} {output.pic} {output.bondary} 

      rm -rf ../make_graphs/out/{wildcards.id}/temp_virv.txt
      rm -rf ../make_graphs/out/{wildcards.id}/tmp_fileout.txt
      rm -rf ../make_graphs/out/{wildcards.id}/temp_virv_graph.txt
      rm -rf ../make_graphs/out/{wildcards.id}/temp_bins.txt
     '''

rule get_graph:
    input: 
        bins_info = "../data/metawrap/{id}_mtw",
        good_file = "../make_graphs/out/{id}/good.csv",
        plasmid_file = "../find_plasmids/output/samples_output/{id}_plasmids.txt",
        checkm = "../data/metawrap/{id}_mtw_stats",
        gtdb = "../data/gtdbtk/{id}_gtdbtk_sorted.txt",
        thr = rules.find_threshold.output.bondary
    output: 
         graph = "../make_graphs/out/{id}/graph/{id}_graph_table.txt",
         bins_thr_filtering = "../make_graphs/out/{id}/{id}_bins.txt"
    conda:"exp"
    shell:'''
         cat {input.bins_info}|awk '{{print $2":"$1}}'> ../data/metawrap/{wildcards.id}_bins_info.txt
         cat {input.checkm}| cut -f1,2,3 > ../make_graphs/out/{wildcards.id}/{wildcards.id}_checkm_modified.txt
         python3 ../scripts/graph_table_retrieving_sp_mtw.py ../data/metawrap/{wildcards.id}_bins_info.txt {input.good_file} {input.plasmid_file} ../make_graphs/out/{wildcards.id}/{wildcards.id}_checkm_modified.txt {input.gtdb} {input.thr} {output.graph} {output.bins_thr_filtering}
 '''
