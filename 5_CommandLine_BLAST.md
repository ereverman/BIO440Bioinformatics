### BLAST and Working with Conda Environments
* Installing software and tools used in bioinformatics projects can be nearly as complex as running the analyses themselves.
* We can use additional tools to help streamline the installation process

## Conda
* powerful package manager and environment manager that you can use at the Anaconda prompt (windows) or teminal (mac)
* Search for anaconda prompt on windows from the start menu
* Verify that conda is installed:
```
conda --version
```

* The primary advantages of conda are:
  * Creating separate environments that cater to specific analyses
  * Avoiding discrepancies between incompatible packages (python 2 vs 3)
  * Creating a snapshot of the versions of software used in an analysis.
* Create a new environment for BLAST:
```
conda create -n blast_env --yes
# provides information on where to find your conda environments
# and how to activate and deactivate the environment
# It is also possible to install packages at the same time as making the environment 
# but this gets a little messy if you need to install via a special channel

# activate and move into the environment:
conda activate blast_env
```
* install BLAST through bioconda
* bioconda: https://bioconda.github.io
* Makes installing bioinformatics packages which tend to have complex dependencies more streamlined 
* Hosted in coordination with GitHub 
* Can browse packages on the website, but the list tends to be exhaustive.

```
# install blast with bioconda:
conda install -c bioconda blast --yes
ls
# you still have access to your directory system and can save files as usual
```

* To check what conda environments exist:
```
conda info --envs
# * indicates what environment you're in
```

* Conda automatically handles dependencies and versions:
```
python --version

conda deactivate # to close and move out of the environment

python --version
```

* Check to see what packages are installed in an environment:
```
conda activate blast-env
conda list
```

* Create a spec list file that can be shared with collaborators or that can be used to recreate the exact environment:
```
conda list --explicit > ~/Desktop/GeneralScripts_Setup/blast-spec-file.txt
cat ~/Desktop/GeneralScripts_Setup/blast-spec-file.txt

conda create -n blast-copy --file ~/Desktop/GeneralScripts_Setup/blast-spec-file.txt
conda install -n blast-copy --file ~/Desktop/GeneralScripts_Setup/blast-spec-file.txt
# Conda doesn't run checks when installing from a spec file
# Method isn't perfect and can be susceptible to differences in operating system
# Not the preferred method for working on HPC because the architecture of the computing systems are so different
```

### Using conda on KU HPC:
* The same major guidelines apply with some upfront modifications
* You need access to where conda environments are stored
```
module load anaconda
# similar to loading a package in R in order to access tools

conda info --envs

# Create a condarc file:
cd 
ls -lha

# Move to scratch
cd scratch
conda create -n blast-env --yes # takes a while for this to complete, complete for prep homework


# Contents of the .condarc file:
# This is a sample .condarc file.
# It adds the r Anaconda.org channel and enables
# the show_channel_urls option.

# channel locations. These override conda defaults, i.e., conda will
# search *only* the channels listed here, in the order given.
# Use "defaults" to automatically include all default channels.
# Non-url channels will be interpreted as Anaconda.org usernames
# (this can be changed by modifying the channel_alias key; see below).
# The default is just 'defaults'.
channels:
  - conda-forge
  - bioconda
  - anaconda
  - defaults

# Show channel URLs when displaying what is going to be downloaded
# and in 'conda list'. The default is False.
show_channel_urls: None

# For more information about this file see:
# https://conda.io/docs/user-guide/configuration/use-condarc.html

envs_dirs:
  - /panfs/pfs.local/work/sjmac/software/conda_envs
  - /panfs/pfs.local/software/install/anaconda/4.3.11/envs

pkgs_dirs:
 - /panfs/pfs.local/work/sjmac/software/conda_envs/pkgs
```


# Blast
* Basic Local Alignment Search Tool
* Identifies regions of similarity between biological sequences.
* Use ProjectTemplate.sh to generate a new directory: Blast_Example
```
sh GeneralScripts_Setup/ProjectTemplate.sh Blast_Example

cd ~/Desktop/Blast_Example/data
ls 

# Download data from NCBI
# These files are mouse and zebrafish RefSeq protein datasets from NCBI

curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/R_norvegicus/mRNA_Prot/rat.1.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/X_tropicalis/mRNA_Prot/frog.1.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.1.protein.faa.gz

# Each file is compressed.
gunzip *.faa.gz

# Preview one of the files:
head mouse.1.protein.faa
# fasta formatted protein sequences
```

### Blast exercise 1:
* Compare the first two protein sequences in mouse.1.protein.faa against the zebrafish protein data set.
```
# First make a database using the zebrafish data:

makeblastdb -in zebrafish.1.protein.faa -dbtype prot
# generates additional files that have been indexed to facilitate faster searches for sequence similarity

# Second make a smaller subset file to hold example data:
head -n 11 mouse.1.protein.faa > mm-test.faa

# Third call blast to do the search:

# View options:
blastp -help

blastp -query mm-test.faa -db zebrafish.1.protein.faa -outfmt 7 -out ../results/mm-test.x.zebrafish.txt
ls ../mm-test.x.zebrafish.txt
less ../mm-test.x.zebrafish.txt








# Fourth, write a loop that will generate full results for each file comparison:

# Run BLASTP
# my solution:
mkdir ../scripts
mkdir reference
mv zebrafish* reference/

for FILE in ~/Desktop/BLAST_PRACTICE/data/*.protein.faa
do 
echo $FILE
SAMPLE=$(basename ${FILE} .1.protein.faa)
blastp -query ${FILE} \
-db ~/Desktop/BLAST_PRACTICE/zebrafish.1.protein.faa \
-outfmt 7 \
-out ~/DESKTOP/BLAST_PRACTICE/results/${SAMPLE}.x.zebrafish.txt
done
```

### Analyzing results:
* Tab-formated data can be easilty read into R:
```
library(ggplot2)

mouse <- read.table("blast_test/results/mm-100.x.zebrafish.txt", header = FALSE)
head(mouse)
colnames(mouse) <- c("query","subject","pct.identity","alignment.length","mismatches",
                     "gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score")


head(mouse)
hist(mouse$evalue)
hist(mouse$pct.identity)

mouse.high.pct.identity <- subset(mouse, pct.identity  > 80)
```
