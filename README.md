# HIC_amr_plasmids_project

## Review:
### Analysis on plasmidome and antibiotic resistance genes possible transfer using classic shotgun WGS and HIC data 

## Usage:
* Clone the project repo:
```bash
git clone https://github.com/alexchist99/HIC_amr_plasmids_project.git
```
* Import environments from envs:
```bash
cd HIC_amr_plasmids_project/envs
for env in *; do conda env create -f $env; done
cd ..
```
* Make folders with input files (the pathes to source files might be different, so you can edit script):
```bash
cd data/raw_data
ln -s /store/bioinf/data/own/run_2022_01/SD00*_HiC_R[12].fastq.gz .
ln -s /store/bioinf/data/own/run_2021_04_filtered/*_HiC_R[12].fastq.gz .
cd ../..

cd workflow
bash make_folders.sh 
```

* Rename sputum files to simplify the code:
```bash
mv data/raw_data_wgs/B-9782_R[12].fastq.gz data/raw_data_wgs/SD002_B_R[12].fastq.gz
mv data/raw_data_wgs/B-9817_R[12].fastq.gz data/raw_data_wgs/SD008_B_R[12].fastq.gz
```
* Download tools [Plasflow](https://github.com/smaegol/PlasFlow) and [rgi](https://github.com/arpcard/rgi) into *tools/* :
* Download databases [Pfam 34v.](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz) and [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb) into *databases/Pfam/*  and *databases/PLSDB/*  accordingly 

* Then use snakemake (you can parallel the workflow by changing "-c" option (-c10: 10 files are processed at the same time)):
```bash
cd workflow
conda activate snakemake
snakemake --use-conda -c10 
```

## Workflow:

![alt text](make_graphs/out/SD005/SD005_hist.jpg)
