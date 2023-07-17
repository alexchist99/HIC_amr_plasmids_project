import glob
import os

path="/home/alexclear/HIC_amr_plasmids_project/data/raw_data"

files= [os.path.basename(f) for f in glob.glob(f'{path}/*')]


list_files = sorted(list(set([i.split("_")[0] for i in set(files)])))
print(list_files)



rule all:
  input: "../tree_vs_plasmidome_resistome/models_article.pdf"


rule preprocessing:
  input:  expand("../make_graphs/out/{id}/fileout.txt", id=list_files)
  output: "../tree_vs_plasmidome_resistome/fileout_thr.txt"
  shell: '''
  cp ../make_graphs/out/*/graph/*.txt ../tree_vs_plasmidome_resistome/
  cat {input} > {output}
        '''
  

rule get_adj_matrix:
   input:  rules.preprocessing.output
   output:"../tree_vs_plasmidome_resistome/amr_merged.csv"
   conda:"exp"
   message: "get adjacent table... "
   shell: '''
   python3 ../scripts/get_adjacency_matrix.py {input} ../tree_vs_plasmidome_resistome {output}
   rm ../tree_vs_plasmidome_resistome/*_data.csv 
   rm ../tree_vs_plasmidome_resistome/*_graph_table.txt 
   '''


rule get_amr_adj_matrix:
   input: 
      matrix = rules.get_adj_matrix.output,
      tree = "../make_tree/output/phylogenomic-tree.txt"
   output: "../tree_vs_plasmidome_resistome/models_article.pdf"
   conda:"exp"
   shell: ''' 
   Rscript ../scripts/AR_distribution_hic_SPADES.R -m {input.matrix} -t {input.tree} -o {output}
   rm Rplots.pdf
   '''

