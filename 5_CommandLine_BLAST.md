change command prompt:
```
PS1='\w $ '

```
### BLAST and Working with Conda Environments
* Installing software and tools used in bioinformatics projects can be nearly as complex as running the analyses themselves.
* We can use additional tools to help streamline the installation process

## Conda
* powerful package manager and environment manager that you can use in gitbash or terminal.


* First, we need to establish a connection between anaconda code and gitbash.

```
cd ~
ls # look for anaconda3/

cd anaconda3/etc/profile.d
ls # We need a copy of the conda.sh script in our ~/.bash_profile file

echo ". ${PWD}/conda.sh" >> ~/.bash_profile

# close and reopen gitbash for the change to take effect
```

* Verify that conda is installed:
```
conda --version
```

* The primary advantages of conda are:
  * Creating separate environments that cater to specific analyses
  * Avoiding discrepancies between incompatible packages (python 2 vs 3)
  * Creating a snapshot of the versions of software used in an analysis.
* Go here to check for packages: https://anaconda.org/anaconda/repo
* There isn't a blast package available for windows (and many other genomics packages) so we will make a practice environment

* Create a new environment for python2.7:


```
conda create --name python2.7 python=2.7

```





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
* Removing conda environments:
```
conda remove --name blast-copy --all

# Prompts for verification
conda info --envs
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

## Conda environment script:
* interacting with hpc can be inefficient
* Using scripts can save time or at least free up your terminal so you can work on other things:

```
#!/bin/bash
#SBATCH --partition=sixhour          # Partition Name (Required)
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=8            # Number of CPU cores per task should match the NTHREADS variable
#SBATCH --mem-per-cpu=60g            # Job memory request
#SBATCH --time=06:00:00              # Time limit hrs:min:sec
#SBATCH --output=/panfs/pfs.local/scratch/sjmac/e284e911/CondaEnv.out  # Error and Out stream
#SBATCH --job-name=CondaEnv


# Usage instructions:
# Edit the last line with unique ENVNAME and PKGNAME
# sbatch CondaEnv.sh

module load anaconda

activate_conda_environment () {
  ENVNAME=${1}
  PKGNAME=${2}
  # expected directory for the conda environment
  ENVDIR=/panfs/pfs.local/work/sjmac/software/conda_envs/${ENVNAME}
  # check if the directory exists; if it doesn't, create the env
  if [ ! -d "${ENVDIR}" ]; then
    # create a conda environment for the application
    # --yes = do not ask for confirmation during installation
    conda create --name ${ENVNAME} --yes
    # install the application in the newly created conda environment
    # --yes = do not ask for confirmation during installation
    conda install --name ${ENVNAME} ${PKGNAME} --yes
  fi
  conda activate ${ENVNAME}
}

activate_conda_environment BlastEnv blast



# Example usage:
sbatch PATH/TO/CondaEnv.sh BlastEnv blast
```



# Blast
* Basic Local Alignment Search Tool
* Identifies regions of similarity between biological sequences.
* Use ProjectTemplate.sh to generate a new directory: BLAST_Practice
```
# Already done:
sh GeneralScripts_Setup/ProjectTemplate.sh BLAST_Practice

cd ~/Desktop/BLAST_Practice/data
ls 

```
* Start an interactive job to run processes line by line
* Activate conda environment once in the interactive job
* Edit your prompt to shorten it
```
srun --time=12:00:00 --ntasks=1 --cpus-per-task=16 --mem=125g --partition=sjmac --pty /bin/bash -l
module load anaconda
conda activate BlastEnv
PS1='\w $ '
```

# Download data from NCBI
# These files are mouse and zebrafish RefSeq protein datasets from NCBI
```
curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/R_norvegicus/mRNA_Prot/rat.1.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/X_tropicalis/mRNA_Prot/frog.1.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.1.protein.faa.gz
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/rat.100.protein.faa.gz
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/mouse.100.protein.faa.gz
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/frog.100.protein.faa.gz

# Each file is compressed.
gunzip *.faa.gz

# Preview one of the files:
head mouse.1.protein.faa
# fasta formatted protein sequences


# Do some clean up:
# Zebrafish is our reference, so move that file to the reference directory.
mv zebrafish.1.protein.faa.gz reference/
```

### Blast exercise 1:
* Compare the first two protein sequences in mouse.1.protein.faa against the zebrafish protein data set.
```
# First make a database using the zebrafish data:
cd reference/

makeblastdb -in zebrafish.1.protein.faa -dbtype prot
# generates additional files that have been indexed to facilitate faster searches for sequence similarity




# Second make a smaller subset file to hold example data:
cd ../
head -n 11 mouse.1.protein.faa > mm-test.faa

# Third call blast to do the search:

# View options:
blastp -help

blastp -query mm-test.faa -db zebrafish.1.protein.faa -out ../results/mm-test.x.zebrafish.txt
# Familiar output, but hard to do much additional analysis with.

blastp -query mm-test.faa -db zebrafish.1.protein.faa -outfmt 6 -out ../results/mm-test.x.zebrafish.txt
ls ../mm-test.x.zebrafish.txt
less ../mm-test.x.zebrafish.txt








# Fourth, write a loop that will generate full results for each file comparison for 100 sequence subset:

# Run BLASTP
# my solution:
mkdir ../scripts
mkdir reference
mv zebrafish* reference/

for FILE in ~/Desktop/BLAST_PRACTICE/data/*100.protein.faa
do 
echo $FILE
SAMPLE=$(basename ${FILE} .100.protein.faa)
blastp -query ${FILE} \
-db ~/Desktop/BLAST_PRACTICE/zebrafish.1.protein.faa \
-outfmt 6 \
-out ~/DESKTOP/BLAST_PRACTICE/results/${SAMPLE}_100.x.zebrafish.txt
done


sh blastloop.sh
```

### Analyzing results:
* Tab-formated data can be easilty read into R:
```
# Analyzing BLASTP output:

# Libraries:----
pacman::p_load(ggplot2)
pacman::p_load(Rmisc)
##
# Datafiles and formatting:----

frog <- read.table("BLAST_Project/results/frog_100.x.zebrafish.txt", header = FALSE)
mouse <- read.table("BLAST_Project/results/mouse_100.x.zebrafish.txt", header = FALSE)
rat <- read.table("BLAST_Project/results/rat_100.x.zebrafish.txt", header = FALSE)

head(frog)
head(mouse)
head(rat)

colnames(frog) <- c("query","subject","pct.identity","alignment.length","mismatches","gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score")
colnames(mouse) <- c("query","subject","pct.identity","alignment.length","mismatches","gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score")
colnames(rat) <- c("query","subject","pct.identity","alignment.length","mismatches","gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score")

head(frog)
head(mouse)
head(rat)

##
# Explore data with simple plots:----
hist(frog$evalue) # Lower numbers indicate more significant simiarity in comparison
hist(frog$pct.identity) # Higher numbers indicate more similarity
hist(frog$bit.score) # Higher numbers indicate more similarity

hist(mouse$evalue) # Lower numbers indicate more significant simiarity in comparison
hist(mouse$pct.identity) # Higher numbers indicate more similarity
hist(mouse$bit.score) # Higher numbers indicate more similarity

hist(rat$evalue) # Lower numbers indicate more significant simiarity in comparison
hist(rat$pct.identity) # Higher numbers indicate more similarity
hist(rat$bit.score) # Higher numbers indicate more similarity



# Some of the comparisons are more similar than others, but the pattern isn't straightforward
# There is a normal looking distribution of values and then a bump at the upper end that 
# is more pronounced for frogs and a little in rat. 


# What questions come from the comparisons:
# Does one species share more similar genes than the others? How do we test this?
# What are the proteins that are most simlar across all species?
# What are the proteins that are uniquely similar between the frog and zebrafish?

##
# Testing the level of similarity between species:----
# pct.identity is possibly useful:

# First need a single dataset:
species.comparison <- data.frame(species = c(rep("frog", nrow(frog)),
                                             rep("mouse",nrow(mouse)),
                                             rep("rat",nrow(rat))),
                                 bit.score = c(frog$bit.score,
                                                  mouse$bit.score,
                                                  rat$bit.score))

head(species.comparison)

sc.aov <- aov(bit.score ~ species, species.comparison)
anova(sc.aov)

sc.means <- summarySE(species.comparison, measurevar = "bit.score", groupvars = "species")
head(sc.means)

ggplot(sc.means, aes(species, bit.score)) +
  geom_errorbar(aes(ymin = bit.score - ci, ymax = bit.score + ci), width = 0.1) +
  geom_point() +
  theme_classic()

ggplot(species.comparison, aes(bit.score)) +
  geom_density(aes(fill = species), alpha = 0.5) +
  theme_classic()


ggplot(species.comparison, aes(bit.score)) +
  geom_density(aes(fill = species), alpha = 0.5) +
  geom_vline(xintercept = 100) +
  theme_classic()

ggplot(subset(species.comparison, bit.score > 100), aes(bit.score)) +
  geom_density(aes(fill = species), alpha = 0.5) +
  theme_classic()

# Is the mouse result biologically significant?
# Some level of similarity is expected. Why?

##
# Identify proteins that have high levels of similarity:----

frog.high.bit.score <- subset(frog, evalue  < 0.05)
mouse.high.bit.score <- subset(mouse, bit.score  > 600)
rat.high.bit.score <- subset(rat, bit.score  > 900)

length(unique(mouse.high.bit.score$query))

```
