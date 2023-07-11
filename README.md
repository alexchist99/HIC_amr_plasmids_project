### WARNING!:
Initial files are available only on ibg server 

### Folders:


### Usage:
* Clone the project repo:
```bash
git clone 
```
* Import environments from envs:
```bash
cd envs
for env in *; do conda env create -f $env; done
```
* Then use:
```bash
cd workflow
conda activate snakemake
snakemake --use-conda -c1
```
