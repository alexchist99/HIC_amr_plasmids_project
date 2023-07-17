mkdir ../data/metabat
mkdir ../data/metawrap
mkdir ../data/spades_contigs_out
mkdir ../data/gtdbtk


for i in "SD001" "SD002" "SD003" "SD005" "SD006" "SD008" "SD009";
do 
ln -s /store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/$i/spades_out/contigs.fasta ../data/spades_contigs_out/${i}_contigs.fasta
ln -s /store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/$i/depth_metabat.tsv ../data/metabat/${i}_depth_metabat.tsv
ln -s /store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/$i/metawrap/metawrap_50_20_bins.stats ../data/metawrap/${i}_mtw_stats
ln -s /store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/$i/metawrap/gtdbtk_out/classify/gtdbtk_sorted.txt ../data/gtdbtk/${i}_gtdbtk_sorted.txt
ln -s /store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/$i/metawrap/metawrap_50_20_bins.contigs ../data/metawrap/${i}_mtw
done


for i in "IL" "Kas" "Sys" "Mez" "Roj";
do 
ln -s /store/bioinf/analysis/hicmicrobiome/5_covid_samples_new/$i/spades_out/contigs.fasta ../data/spades_contigs_out/${i}_contigs.fasta
ln -s /store/bioinf/analysis/hicmicrobiome/5_covid_samples_new/$i/depth_metabat.tsv  ../data/metabat/${i}_depth_metabat.tsv
ln -s /store/bioinf/analysis/hicmicrobiome/5_covid_samples_new/$i/metawrap/metawrap_50_20_bins.stats ../data/metawrap/${i}_mtw_stats
ln -s /store/bioinf/analysis/hicmicrobiome/5_covid_samples_new/$i/metawrap/gtdbtk_out/classify/gtdbtk_sorted.txt ../data/gtdbtk/${i}_gtdbtk_sorted.txt
ln -s /store/bioinf/analysis/hicmicrobiome/5_covid_samples_new/$i/metawrap/metawrap_50_20_bins.contigs ../data/metawrap/${i}_mtw
done