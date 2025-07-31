
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
conda env create -f RGnet_PM3/PM3.yml
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
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr13.vcf.bgz.tbi
```
The database and Scripts directory contain the data and conf for a full example of the GJB2 gene. To run with a complete genome, users needs download the appropriate databases and reverse the file name in the vcfanno.conf file.


# Usage
Users can enter the following command to view the use of this method and define the parameters themselves.
```
python ./Scripts/pm3.py --help
```
```
usage: pm3.py [-h] --vcf VCF --bed BED --sam SAM [--cache CACHE] --gv GV --fa FA --conf CONF [--af AF] --plp PLP [--blb BLB] [--incis INCIS]
                        [--transcript TRANSCRIPT] [--ncd {Yes,No}] [--intrans INTRANS] --scorePy SCOREPY --output OUTPUT

Extract AR genes and case/family samples, VEP annotation , gnomAD_joint_grpmax_AF annotation, variants filtration, RGnet construction, PM3 tagging.

options:
  -h, --help            show this help message and exit
  --vcf VCF             The vcf file.
  --bed BED             The bed file of AR genes.
  --sam SAM             The case or family sample file.
  --cache CACHE         The path of VEP cache files.
  --gv GV               The version of genome.
  --fa FA               The reference fasta file.
  --conf CONF           The configure file of vcfanno.
  --af AF               The allele frequency of PM2.
  --plp PLP             The P/LP file.
  --blb BLB             The B/LB file.
  --incis INCIS         The in cis file.
  --transcript TRANSCRIPT
                        ENST IDs file of the transcript of the target gene.
  --ncd {Yes,No}        Whether to include non-coding region variants.
  --intrans INTRANS     The in trans file.
  --scorePy SCOREPY     The scoring script.
  --output OUTPUT       The name of output file.
```



# Example
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

```
# Example of the output file
CHROM   POS     REF     ALT     Point_Hom       Symbol  Transcript      VEP_IMPACT      VEP_Consequence AF_Max  VEP_HGVSc       VEP_HGVSp       HomozygousNum      Score   PM3_Level
11      67286437        G       A       0       CABP2   ENST00000294288 MODIFIER        3_prime_UTR_variant     0.00013659      ENST00000294288.4:c.*173C>T0       0.0     NA
11      67287257        G       C       0       CABP2   ENST00000294288 LOW     splice_region_variant&intron_variant    0.0009362       ENST00000294288.4:c.637+7C>G               0       0.0     NA
11      67287263        C       A       0.5     CABP2   ENST00000294288 HIGH    splice_donor_variant    0.0013616       ENST00000294288.4:c.637+1G>T      10.5     PM3_Supporting
11      67287517        C       T       0       CABP2   ENST00000294288 LOW     splice_donor_5th_base_variant&intron_variant    0.00013373      ENST00000294288.4:c.489+5G>A               0       0.0     NA
11      67290852        G       A       0       CABP2   ENST00000294288 MODIFIER        5_prime_UTR_variant     5.7657e-05      ENST00000294288.4:c.-55C>T00.0     NA
11      67287523        TC      T       1       CABP2   ENST00000294288 HIGH    frameshift_variant&splice_region_variant        .       ENST00000294288.4:c.487del ENSP00000294288.4:p.Glu163SerfsTer55    2       1.5     PM3_Moderate
11      67288601        G       A       0       CABP2   ENST00000294288 HIGH    stop_gained     1.3359e-05      ENST00000294288.4:c.274C>T      ENSP00000294288.4:p.Arg92Ter       0       0.5     PM3_Supporting
11      67288539        G       A       0       CABP2   ENST00000294288 LOW     synonymous_variant      0.00031359      ENST00000294288.4:c.336C>T      ENSP00000294288.4:p.Tyr112%3D      0       0.0     NA
11      67290766        C       T       0       CABP2   ENST00000294288 MODERATE        missense_variant        0.0020017       ENST00000294288.4:c.32G>A ENSP00000294288.4:p.Arg11Gln     0       0.0     NA
```


