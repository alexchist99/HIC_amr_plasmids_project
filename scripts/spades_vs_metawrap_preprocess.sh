#!/bin/bash 

#get and rename fasta with metawrap nodes only
seqkit seq -w 0 $1 | awk '/>/ {getline seq} {sub (">","",$0);print $0, seq}' | sort -k1 | join -1 1 -2 1 - <(sort -k1 $2 ) -o 1.1,2.2,1.2| awk '{print ">"$1":"$2"\n"$3}' | seqkit seq -w 60 > $3_mtawrap_only.fa


#get list of headers 
cat $2|cut -f1  > $3_metawrap_cont_first

cat $1|grep ">"|cut -d ">" -f2 > $3_contig_name

#get fasta with spades nodes only
grep -Fxv -f $3_metawrap_cont_first $3_contig_name > $3_only_spades_nodes  #get headers only in sd002_contig_name 
seqtk subseq $1 $3_only_spades_nodes > $3_only_spades_nodes.fa #get fasta with only spades nodes
sed '/^>/s/$/:bin.0/' $3_only_spades_nodes.fa > $3_spades_only.fa  #add node.0:

#concatenate
cat $3_mtawrap_only.fa $3_spades_only.fa > $4

mv $3_mtawrap_only.fa $5

#delete temp files
rm $3_metawrap_cont_first
rm $3_contig_name
rm $3_only_spades_nodes
rm $3_only_spades_nodes.fa
rm $3_spades_only.fa

