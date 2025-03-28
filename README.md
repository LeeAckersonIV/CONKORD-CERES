## **CONKORD on USDA CERES**
**This file contains the pigT2T methodology, code, and logistic reasoning for employing the *CONKORD* pipeline on *USDA-CERES*. The CONKORD pipeline was developed to determine CNV of rDNA morphs by the primateT2T group. Note: In the *Ape Genome* paper, CONKORD Version 7 was employed, here we use the most recent release - version 8.**
* This pipeline is described in the supplementary notes of this paper: https://doi.org/10.1101/2024.07.31.605654
* The GitHub repository containing the source code is located here: https://github.com/borcherm/CONKORD (**Matthew Borchers**)

File Author: *Lee Ackerson*; ackers24@msu.edu

___

### **Description of the CONKORD pipeline**

A snakemake pipeline for counting/estimating the copy number of genomic features. The python script `conkord.py` creates a config.yml file for snakemake, which then runs a "Snakemake all" command that executes several subprocesses.

**Example Usage:**

```
python conkord.py \
    --no_uniq \
    -k 31 \
    -bed feature_coordinates.bed \
    -f feature.fa \
    -r sequencing_data/ \
    -t 15 \
    --cluster \
    -g genome.fa \
    --gzip 
```

**Paramter Description:**

* `--no_uniq` Use when evaluating rDNA; rarely used for other genomic features

* `-k` kmer size; {default=31}

* `-f` The fasta for your feature of interest

* `-bed` A bed file with coordinates for your feature of interest in the reference genome.

* `-r` Directory with only the sequencing reads in fastq or fastq.gz format

* `-g` The reference genome

* `-gzip` If the reads are gzipped you must pass this flag

* `-t` Number of threads

* `-w_size` The window size for finding G/C matched windows to normalize to (default 2000)

* `--cluster` If you are using a cluster with slurm this will cause conkord.py to submit an sbatch command for snakemake with the number of cores provided

___

### **Prepare Inputs for the CONKORD pipeline [Pig HxYL T2T Assembly]**

#### **1. Generate rDNA morphs for each haplotype & chromosome via ribotin**
```
# run ribotin to determine consensus.fa for each rDNA morph
ribotin-verkko \
  --approx-morphsize 15000 \ 
  --guess-tangles-using-reference Pig_rDNA.fasta \
  -i Pig/assembly/verkko2.2_hifi-duplex_trio \
  -o Pig/assembly/ribotin-verkko_pig-rDNA_tangle
```

#### **2. Manually generate BED files for each haplotype assembly**
```
# --------------------------------------------------- #
# run repeat masker to identify the location of rDNA  #
# --------------------------------------------------- #

## ==== Repeat Masker ==== ##

RepeatMasker -species "pig" \
    -libdir /project/cattle_genome_assemblies/config_files_scripts/RepeatMasker_4.0.6_lib \
    -no_is \
    -pa 96 \
    -gff \
    -s \
    -dir RM_hap1 \
    assembly.haplotype1.fasta

RepeatMasker -species "pig" \
    -libdir /project/cattle_genome_assemblies/config_files_scripts/RepeatMasker_4.0.6_lib \
    -no_is \
    -pa 96 \
    -gff \
    -s \
    -dir RM_hap2 \
    assembly.haplotype2.fasta
```

**Leverage output files to extract the hapmers containing/flanking the rDNA region:**

```
# ======================================== #
# PigT2T Assembly Example                  #
## verkko output files needed:             #
###     assembly.paths.tsv                 #
###     assembly.scfmap.tsv                #
###     assembly.haplotype{1,2}.fasta.gaps #
# ======================================== #

# one of the rDNA morph flanking utigs in the pig assembly was: utig4-58+ (maternal)
# ----------------------------------------------------------------------------------
grep utig4-58+ assembly.paths.tsv | awk '{print $1}'
    # > dam_compressed.k31.hapmer_from_utig4-1449

grep dam_compressed.k31.hapmer_from_utig4-1449 assembly.scfmap.tsv
    # > path	dam_compressed.k31.hapmer-0000007	dam_compressed.k31.hapmer_from_utig4-1449

grep dam_compressed.k31.hapmer-0000007 assembly.haplotype2.fasta.gaps
    # > dam_compressed.k31.hapmer-0000007	55395822	55403322
```

The range [55395822, 55403322] represents the gap in the verkko assembly caused by the rDNA tangle. However, we want the entire range of the rDNA region for CONKORD. We can get this by manually inspecting the RepeatMasker output and taking the **first** and **last** coorindates having rDNA annotation for the query hapmer.

```
# change to repeat masker outout dir; for haplotpye of interest
cd RepeatMasker/RM_Hap1/

# find contig comprising rDNA tangle; this will print entire table (manually investigate)
grep dam_compressed.k31.hapmer-0000007 assembly.haplotype1.fasta.out.rRNA.filtered 

# get the first coordinate for the forward strand (+):
grep dam_compressed.k31.hapmer-0000007 assembly.haplotype1.fasta.out.rRNA.filtered | awk '$9 == "+" {print $6}' | head -n 1
    # > 55333083
	
# get the last coordinate for the forward strand (+):
grep dam_compressed.k31.hapmer-0000007 assembly.haplotype1.fasta.out.rRNA.filtered | awk '$9 == "+" {print $7}' | tail -n 1
    # > 55433905
```
We now know that for the rDNA tangle on `dam_compressed.k31.hapmer-0000007` (Chr8), rDNA can be found between coordinates 55333083 and 55433905. In the interest or robustness, and to ensure that the entire region is encapsulated - I rounded these coordinates to a nearby number within ~250 bp:
* Start coordinate rounded down: 55333083 -> 55333000
* End coordinate rounded up: 55433905 -> 55434000

This was done such that some sequence *NOT* containing rDNA was included in the start and stop coordinate region, and all information was placed into a bed file:
```
echo -e 'dam_compressed.k31.hapmer-0000007\t55333000\t55434000' > assembly.haplotype1.chr8rDNA.bed
```

The above process is repeated for each tangle in each haplotype such that every hapmer has its own line in the BED file. The BED file should be formatted as such:

    Chr/Hapmer# | Start Coordinate rDNA region | End Coordinate rDNA region

Note: If there is a gap/break in your assembly at the point of an rDNA tangle - meaning the two rDNA flanking utigs belong to seperate hapmers/contigs (which was the case for the pigT2T assembly), then both hapmers should be included in the BED file - ex. the paternal Chr8 rDNA morph exists between/on `sire_compressed.k31.hapmer-0000432` & `sire_compressed.k31.hapmer-0000446`

So, even though there is 2 rDNA arrays/tangles for the sire, the `assembly.haplotype2.chr8rDNA.bed` file will have 2 lines (1 for each Chr/Contig going into the Chr8 array). This is important because since there is a break here, the coorindate system changes.


#### **3. Set Up Working *conkord* Directory**

```
# make working directory
mkdir conkord
cd conkord

# ribotin tangles are associated with respective tangles via nodes.txt
scp ribotin-verkko_pig-rDNA_tangle0/consensus.fa conkord/chr8.hap1.consensus.fa
scp ribotin-verkko_pig-rDNA_tangle1/consensus.fa conkord/chr10.hap1.consensus.fa
scp ribotin-verkko_pig-rDNA_tangle2/consensus.fa conkord/chr10.hap2.consensus.fa
scp ribotin-verkko_pig-rDNA_tangle3/consensus.fa conkord/chr8.hap2.consensus.fa

# assembly.haplotype{1,2}.chr{8,10}rDNA.bed files
echo -e 'dam_compressed.k31.hapmer-0000007\t55333000\t55434000' >> assembly.haplotype1.chr8rDNA.bed
echo -e 'dam_compressed.k31.hapmer-0000014\t54086100\t54205400' >> assembly.haplotype1.chr10rDNA.bed
echo -e 'sire_compressed.k31.hapmer-0000432\t0\t224300' >> assembly.haplotype2.chr8rDNA.bed
echo -e 'sire_compressed.k31.hapmer-0000446\t54490200\t54622500' >> assembly.haplotype2.chr8rDNA.bed
echo -e 'sire_compressed.k31.hapmer-0000429\t37666000\t38067500' >> assembly.haplotype2.chr10rDNA.bed

# link verkko assembly and F1 illumina data directories
ln -s project/ruminant_t2t/Pig/illumina_data/F1/ illumina
ln -s /project/ruminant_t2t/Pig/assembly/verkko2.2_hifi-duplex_trio
    
    # Note: the illumina reads need to be in a directory by themselves (no parental data)
    # AND, files need to be named as such: {id}_1.fastq.gz
    # the R1 vs R2 format confuses the string splicing function, use *_1* vs *_2* instead!
    mv Fetus-AM1-PigT2T_S4_L004_R1_001.fastq.gz fetus_1.fastq.gz 
	mv Fetus-AM1-PigT2T_S4_L004_R2_001.fastq.gz fetus_2.fastq.gz

```
___

### **Run the CONKORD pipeline [Pig HxYL T2T Assembly]**

#### **1. Set Up Environment**
```
# navigate to working directory generated previously
cd Pig/assembly/conkord/

# import CONKORD software
git clone https://github.com/borcherm/CONKORD.git
cd CONKORD/

# allocate compute resources
salloc -p priority -q agil -c 96 --mem-per-cpu=3968 --time=4-00:00:00

# activate snakemake env
micromamba activate verkko-v2.2 # verkko has built in snakemake 

# load necessary modules
module load bedtools
module load jellyfish2 # I used jellyfish 2.2.9, old version will throw errors!
module load samtools

# if you dont have the necessary python packages, install them with pip
pip install numpy matplotlib seaborn # needed for Call_Copy_Number_GC_Normalization_Version8.py
```

#### **2. Run CONKORD**

**Run Parameter Specifications:**

* `--no_uniq` Used the default ("on"), as this run was for rDNA features

* `-k 31` Used default k-mer size, in accordance with the literature as a reasonable length

* `-f chr{8,10}.hap{1,2}.consensus.fa` Ribotin consensus.fa output file; fasta of consensus morph for each tangle

* `-bed assembly.haplotype{1,2}.chr{8,10}rDNA.bed` BED file formatted as: Chr/Hapmer# | Start Coordinate rDNA region | End Coordinate rDNA region

* `-r /project/ruminant_t2t/Pig/illumina_data/F1/` Directory containing illumina reads for the F1 individual

* `-g assembly.haplotype{1,2}.fasta` Verkko assembly for each individual haplotype (there will be multiple runs of CONKORD, one for each rDNA tangle on each haplotype)

* `-gzip` Illumina reads gzipped, not needed if not compressed

* `-t 15` Number of threads used on USDA ceres

* `-w_size 30500` Approximate length of pig rDNA morph (get from consensus.fa or reference.fa)
* `--cluster #` Indicate that script is being executed on USDA ceres


**Run CONKORD: Haplotype1-Chr8**
```
python conkord.py \
	--no_uniq \
	-k 31 \
	-bed ../assembly.haplotype1.chr8rDNA.bed \
	-f ../chr8.hap1.consensus.fa \
	-r ../illumina-reads/ \
	-g ../assembly.haplotype1.fa \
	-w_size 30500 \
	-t 15 \
	--cluster \
	--gzip
```
**Run CONKORD: Haplotype1-Chr10**
```
python conkord.py \
	--no_uniq \
	-k 31 \
	-bed ../assembly.haplotype1.chr10rDNA.bed \
	-f ../chr10.hap1.consensus.fa \
	-r ../illumina-reads/ \
	-g ../assembly.haplotype1.fa \
	-w_size 30500 \
	-t 15 \
	--cluster \
	--gzip
```
**Run CONKORD: Haplotype2-Chr8**
```
python conkord.py \
	--no_uniq \
	-k 31 \
	-bed ../assembly.haplotype2.chr8rDNA.bed \
	-f ../chr8.hap2.consensus.fa \
	-r ../illumina-reads/ \
	-g ../assembly.haplotype2.fa \
	-w_size 30500 \
	-t 15 \
	--cluster \
	--gzip
```
**Run CONKORD: Haplotype2-Chr10**
```
python conkord.py \
	--no_uniq \
	-k 31 \
	-bed ../assembly.haplotype2.chr10rDNA.bed \
	-f ../chr10.hap2.consensus.fa \
	-r ../illumina-reads/ \
	-g ../assembly.haplotype2.fa \
	-w_size 30500 \
	-t 15 \
	--cluster \
	--gzip
```
___


### **Output of the CONKORD Pipeline**

**`Conkord` outputs a variety of intermediate files and graphs for the input data at hand, but we are primarily interested in the contents of the `results/` folder. This folder should comprise one file for each succesful run of the pipeline (4 files w/ above code), with the following nomenclature: `Copy_Numbers_nu_(feature-ID}_k{}.tsv`. These files contain the median and mean copy number estimates of the rDNA morphs.
```
cd CONKORD/results/
cat *

# > Sample	                    Median_Haploid_Copy_Number	Median_Diploid_Copy_Numer	Mean_Haploid_Copy_Number	Mean_Diploid_Copy_Number
# > chr10.hap1.consensus_fetus  110                         220                         116	                        233
# > chr10.hap2.consensus_fetus	116	                        232	                        121	                        243
# > chr8.hap1.consensus_fetus	108	                        216	                        112	                        224
# > chr8.hap2.consensus_fetus	109	                        218	                        112	                        225
```

**In the pigT2T assembly, we elected to use the *Median Haploid Copy Number* for our estimate At this point, we feed a patch file back to *verkko* containing the ribotin consensus.fa morphs with the now estimated copy numbers to reolve the tangled region.**
