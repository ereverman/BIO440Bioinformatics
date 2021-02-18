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
