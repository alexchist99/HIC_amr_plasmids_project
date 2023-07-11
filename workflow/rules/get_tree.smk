import glob
import os



#get the tree

path="/home/alexclear/HIC_amr_plasmids_project/data/raw_data"

files= [os.path.basename(f) for f in glob.glob(f'{path}/*')]

list_files = sorted(list(set([i.split("_")[0] for i in set(files)])))
#print(list_files)



rule all6:
    input: "../make_tree/output/phylogenomic-tree.txt"


rule get_good_MAGS:
    output: "a.txt"
    conda:"exp"
    shell:''' 
#mkdir ../make_tree/input_pangenome
for i in {list_files}
do 
while read -r line
do
grep -w $line ../make_graphs/out/$i/genomes/marked_mtw_mags.faa|cut -f2- -d">" > ../make_tree/${{i}}_temp_file
NAME=$(cat /store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/${{i}}/metawrap/gtdbtk_out/classify/gtdbtk_sorted.txt|awk -F";" '{{for(j=1; j<=NF-1; ++j) printf "%s ", $j; print ""}}'|awk '{{print $1"\t"$NF}}'|grep -w $line|cut -f2|awk -F "__" '{{print $2}}')
seqtk subseq ../make_graphs/out/${{i}}/genomes/marked_mtw_mags.faa ../make_tree/${{i}}_temp_file > ../make_tree/input_pangenome/${{NAME}}.fna
rm ../make_tree/${{i}}_temp_file

done < ../make_graphs/out/$i/${{i}}_bins.txt
done
touch a.txt
     '''


#rules.get_good_MAGS.output
rule make_contig_database:
   input: "a.txt"
   output: "../make_tree/input_pangenome/genomes.txt"
   conda:"anvio-7.1"
   shell: '''

for i in `ls ../make_tree/input_pangenome/*.fna | awk 'BEGIN{{FS=".fna"}}{{print $1}}'`
do  
    anvi-script-reformat-fasta ${{i}}.fna -o ${{i}}.fa -l 0 --simplify-names
    anvi-gen-contigs-database -f ${{i}}.fa -o ${{i}}.db -T 10
    anvi-run-hmms -c ${{i}}.db
done

#get external-genomes file with .fa and .db

(cd ../make_tree/input_pangenome; ls *.fa) > ../make_tree/input_pangenome/file_fa
(cd ../make_tree/input_pangenome; ls *.db) > ../make_tree/input_pangenome/file_db
 
paste ../make_tree/input_pangenome/file_fa ../make_tree/input_pangenome/file_db > ../make_tree/input_pangenome/temp_genomes.txt
echo -e "name\tcontigs_db_path\n$(cat ../make_tree/input_pangenome/temp_genomes.txt)" > {output}

#delete temporary files
rm a.txt
rm ../make_tree/input_pangenome/*.fna
rm ../make_tree/input_pangenome/file_fa
rm ../make_tree/input_pangenome/file_db
rm ../make_tree/input_pangenome/temp_genomes.txt
   '''
   

rule alignment:
  input: rules.make_contig_database.output
  output: "../make_tree/output/concatenated-proteins.fa",
  conda:"anvio-7.1"
  shell: '''

# get aa sequences (aligned)
anvi-get-sequences-for-hmm-hits --external-genomes {input} \
                                -o {output} \
                                --hmm-source Bacteria_71 \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate                 
'''



rule tree:
  input: rules.alignment.output
  output: "../make_tree/output/phylogenomic-tree.txt",
  conda:"anvio-7.1"
  shell: '''
#get tree
anvi-gen-phylogenomic-tree -f {input} \
                           -o {output}     
  '''