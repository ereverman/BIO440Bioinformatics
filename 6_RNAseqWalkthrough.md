## RNAseq Pipeline building

Data come from Everman et al. 2021 Genetics
10 samples exposed to control and copper treatments for 9 hours


## Setting up the project
* ProjectTemplate.sh can be used for more than setting up a project directory structure
* We will use it to:
  * Set up the project
  * Act as a source file that will be referenced by other scripts

```
cd Desktop
mkdir RNAseq_Project

cp GeneralScripts_Setup/ProjectTemplate.sh RNAseq_Project/RNAseq_config.sh

```

* We need to make some modifications that will accomplish:
  * Project-specific directories
  * Establishing variables that we will use in additonal scripts

```
#!/usr/bin/env bash


# RNAseq config template

# Instructions:
# Create a local and remote project directory and modify config.sh
# Modify sftp connection file
# Change lines 17 and 18
# ON HPC, in PROJECT_DIR: sh RNAseq_config.sh PROJECT_NAME


# Assign user variable (anchoring reference point)
ERE_DIR=/panfs/pfs.local/scratch/sjmac/e284e911

# Assign Project-specific variables
PROJECT="RNAseq_Project"
PROJECT_DIR=${ERE_DIR}/RNAseq_Project




# Assign PROJECT_DIR subdirectories

# DATA_DIR
DATA_DIR=${PROJECT_DIR}/data
DATA_RAW=${DATA_DIR}/raw
DATA_FILT=${DATA_DIR}/filtered
DATA_REFS=${DATA_DIR}/refs

# RESULTS_DIR
RESULTS_DIR=${PROJECT_DIR}/results
KALLISTO_DIR=${RESULTS_DIR}/kallisto
STRINGTIE_DIR=${RESULTS_DIR}/stringtie
HISAT_DIR=${RESULTS_DIR}/hisat
SAM_DIR=${HISAT_DIR}/sam
BAM_DIR=${HISAT_DIR}/bam
FCOUNTS_DIR=${RESULTS_DIR}/feature-counts

# REPORTS_DIR
REPORTS_DIR=${PROJECT_DIR}/reports
FASTP_REPORTS_DIR=${REPORTS_DIR}/fastp_reports/

# SCRIPTS_DIR
SCRIPTS_DIR=${PROJECT_DIR}/scripts




# Create directories
mkdir -p ${DATA_DIR} \
         ${DATA_RAW} \
         ${DATA_FILT} \
         ${DATA_REFS} \
         ${RESULTS_DIR} \
         ${KALLISTO_DIR} \
         ${STRINGTIE_DIR} \
         ${HISAT_DIR} \
         ${SAM_DIR} \
         ${BAM_DIR} \
         ${FCOUNTS_DIR} \
         ${REPORTS_DIR} \
         ${FASTP_REPORTS_DIR} \
         ${SCRIPTS_DIR} \
```

* Next, run RNAseq_config.sh on HPC
* This order is important because we hard-coded the path to our user directory. The script wont run on our local computers
* The reason for this decision is that the config file serves the dual purpose of setting up and acting as a source file for other scripts.

```
cd RNAseq_Project
sh RNAseq_config.sh

mv RNAseq_config.sh scripts/
```


## Data download:
* refer to the paper methods section on data availablility to figure out where the data are stored.
* SRA PRJNA633166
* Go to NCBI and search SRA for this project number
* Select all 20 files and click "send to" --> file --> Accession List
* Move the SraAccList.txt to your local project directory and upload to remote via atom

```
# sra toolkit is a package that helps streamline the process of getting data from SRA
# unfortunately doesn't play nicely with cluster

# Create a README file
