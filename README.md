# HIC_amr_plasmids_project

## Review:
### Analysis on plasmidome and antibiotic resistance genes possible transfer using classic shotgun WGS and HIC data 

#### WARNING!:
Initial files are available only on ibg server 

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
* Rename sputum files to simplify the code:
```bash
mv data/raw_data_wgs/B-9782_R[12].fastq.gz data/raw_data_wgs/SD002_B_R[12].fastq.gz
mv data/raw_data_wgs/B-9817_R[12].fastq.gz data/raw_data_wgs/SD008_B_R[12].fastq.gz
```
* Download tools [Plasflow](https://github.com/smaegol/PlasFlow) and [rgi](https://github.com/arpcard/rgi) into * *tools/* *:
* Download databases [Pfam 34v.](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz) and [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb) into * *databases/Pfam/* * and * *databases/PLSDB/* * accordingly 

* Then use:
```bash
cd workflow
conda activate snakemake
snakemake --use-conda -c1
```
