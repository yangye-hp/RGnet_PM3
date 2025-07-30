
# Overview





# Installation

1. PM3 Scripts

Install PM3 from github.
```
git clone https://github.com/yangye-hp/RGnet_PM3.git
```

2. Conda

The recommended way to configure PM3 in a clean environment is by using Conda.
```
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install. Enter "yes" to initialize conda
bash Miniconda3-latest-Linux-x86_64.sh

# Refresh the terminal
source ~/.bashrc

# Configure channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Install PM3 in a new environment named 'PM3':
```
conda env create -f PM3.yml
conda activate PM3
```

3. VEP and VEP plugins
Annotation of gene variants using Ensembl Variant Effect Predictor(VEP) and its plugins. 

```
# Install VEP by conda
conda install ensembl-vep=113

# Download and unpack VEP's offline cache dataset for GRCh38. Make sure that the versions of VEP and VEP cache dataset are consistent
mkdir $HOME/.vep
cd $HOME/.vep
wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz
wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz
tar xzf homo_sapiens_vep_113_GRCh37.tar.gz

# Download VEP's plugins
cd $HOME/.vep
git clone https://github.com/Ensembl/VEP_plugins.git
mv VEP_plugins Plugins
```

4. Reference genome
```
# Download the reference FASTA for GRCh38
mkdir -p $/HOME/RGnet_PM3/database
cd $/HOME/RGnet_PM3/database
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Download the reference FASTA for GRCh37
wget http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
```

5. Allele frequency annotation
```
# Download the Joint Frequency in gnomAD
cd $/HOME/RGnet_PM3/database
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr3.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr3.vcf.bgz.tbi
```
The database and Scripts directory contain the data and conf for a full example of the GJB2 gene. To run with a complete genome, users needs download the appropriate databases and reverse the file name in the vcfanno.conf file.


# Usage


```
cd $/HOME/RGnet_PM3
nohup python ./Scripts/pm3.py \
--vcf ./test/annotated_raw.AF_stats.vcf.gz \
--bed ./inputfile/AR_genelist_hg19_nochr.bed \
--sam ./inputfile/case.txt \
--gv GRCh37 \
--fa ./database/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
--conf ./Scripts/vcfanno.conf \
--plp ./inputfile/PUB_P.txt \
--scorePy ./Scripts/score.py \
--incis ./inputfile/PM3_remove.txt \
--ncd No \
--output output.tsv &
```

# Note
Users can enter the following command to view the use of this method and define the parameters themselves.
```
python ./Scripts/pm3.py --help
```


