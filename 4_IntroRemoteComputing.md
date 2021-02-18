### Accessing KU High Performance Computing Cluster:
* Web resources for KU HPC: https://crc.ku.edu/hpc/how-to
* Need the following installed:
  * Cisco AnyConnect
  * Duo Mobile

### Logging in:
* Open Cisco AnyConnect and enter your username and password. In the second password box, enter push.
* Duo Mobile will ask for authorization, select yes.
* Click Accept on Cisco AnyConnect to finally connect to the secure VPN.


### In terminal/git bash:
* Using a secure shell ssh we can establish a secure connection over our unsecured networks to access KU HPC

```
ssh USERNAME@hpc.crc.ku.edu
# enter password, hit enter
```

### Landing screen:
* Maintenance information--any processes or jobs running during maintenance times will be terminated and HPC can't be accessed during those times.
* Partitions--there are time and memory limits of jobs that we can run in specific partitions of HPC
  * sjmac allows us to run large and long jobs that will take longer than 6 hours to complete
  * sixhour allows jobs that will take less than 6 hours. Any that go over the time limit will be terminated mid-job
* Storage--location for where our files are stored. This is group-specific and is limited. 

```
crctool
# Provides additional information for our resources available
# Paths to data storage are useful
# Partitions shows the actual amount of resources available for each of the partitions
```
* Storage Variables--Locations were we can temporarily store data and outputs.
  * Important to treat HPC resources as temporary storage. Files that aren't accessed regularly are eventually purged and are not recoverable
  * Files that take up a ton of space in SCRATCH (where we do most of our work) that we aren't using will take up usable space from other users
  * Work is more permanent storage that won't be purged, but that space is shared among the sjmac user group and we run into similar space issues occasionally.
  * Best practice is to keep active analyses files in HPC and move to permanent storage files that aren't the focus.



### Basic tools:
```
pwd
ls -lha
# Shows files that will be purged over the last several purges
# There are a couple of files to note:
# .bash_profile allows you to set aliases that can be used as shortcuts to frequently used commands
# .condarc (which may not exist yet) that can be used to create customized and modular environments to run analyses
```

```
cd scratch/
l

pwd # Use to verify working directory
readlink -f . # Use to identify a path that can be used to transfer files by translating a symlink path to an absolute path
```

### Moving files from local to remote computer
* First we need to create a file that we want to transfer to our remote computer.
* We will create a ProjectTemplate.sh script that can be customized for new projects
* We want this script to be kept from purges, so we will store it in our work directory 
* Making a directory in work allows us to specify a connection between our local and remote computers
```
cd work/
mkdir GeneralScripts
```


#### Initial set up for interfacing Atom with our remote computer
* This needs to be done once for each new project you start
* Open Atom and navigate to a single untitled script.
* File --> Add Project Folder --> (Navigate to where you want the matching GeneralScripts directory) New Folder --> "GeneralScripts" --> Open
* File --> Open --> GeneralScripts (Side bar should change to show directory and contents)
* Packages --> Remote FTP --> Create SFTP Config File
  * If you don't have the Remote FTP plugin, install from here: https://atom.io/packages/search?q=ftp
  * May need to restart the project
* Once you have a new .ftpconfig file, go to remote work directory and get the path:
```
readlink -f .
```
* Open a second terminal tab or window (this one is linked to your local computer) and make a matching directory and get the path:
```
mkdir GeneralScripts
pwd
```

* Change host, user, pass, promptForPass, remote, local:
```
{
    "protocol": "sftp",
    "host": "hpc.crc.ku.edu",
    "port": 22,
    "user": "e284e911",
    "pass": "",
    "promptForPass": true,
    "remote": "/panfs/pfs.local/work/sjmac/e284e911/GeneralScripts",
    "local": "/Users/e284e911/Desktop/Teaching/WJC_Bioinformatics_SP2021/5_UnixRemoteComputing/GeneralScripts",
    "agent": "",
    "privatekey": "",
    "passphrase": "",
    "hosthash": "",
    "ignorehost": true,
    "connTimeout": 10000,
    "keepalive": 10000,
    "keyboardInteractive": false,
    "keyboardInteractiveForPass": false,
    "remoteCommand": "",
    "remoteShell": "",
    "watch": [],
    "watchTimeout": 500
}
```
* This file allows us to establish a connection.
* Activate the connection:
* Packages --> Remote FTP --> Toggle
* This opens a new tab in Atom. Click the tab and then click connect, enter your KU password. You should now be able to see the contents of your GeneralScripts directory.


### Make a file that allows us to set up project directories:
* Open a new file in Atom in your local directory.
```
# Project directory template

# Instructions:
# Usage: sh ProjectTemplate.sh PROJECT_NAME
# Parent directory does not need to exist prior to running script
# Run script in the directory you want to create the new project directory

# Assign Project Variable:
# ${1} receives the PROJECT_NAME argument
PROJECT_DIR=${1}


# Assign Project Subdirectory Variables:
DATA_DIR=${PROJECT_DIR}/data
RESULTS_DIR=${PROJECT_DIR}/results


# Make the directories:
mkdir -p ${PROJECT_DIR} \
${DATA_DIR} \
${RESULTS_DIR}
```
* Transfer the file by right clicking the file and selecting upload the remote and local computer.
* Check that the file exists in your HPC connection through gitbash

### Downloading a file that you created or accessed on HPC:
* Download my teaching notes from github
* Transfer them to your local computer
```
cd work/GeneralScripts
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/4_IntroRemoteComputing.md
```
* Return to Atom and view the file in the Remote toggle. Right click and download
* Verify the file is in your local directory.



### Running the ProjectTemplate.sh:
* Remember that whenever you set up a new project and you intend to use a remote computer for analyses, it is simplest to set up the project directories on your local and remote computer simultaneously.
* As you get more familiar with your preferred method for organizing data, your ProjectTemplate.sh will grow as you make changes and additional directories.
* Our first case study on using command line BLAST is coming up, so let's make a project directory for that project


* On your local computer:
```
cd Desktop # (or the place you want to store your BLAST directory)

sh PATH/TO/ProjectTemplate.sh BLAST_PRACTICE
tree BLAST_PRACTICE
tree.com //a BLAST_PRACTICE
```
* Make a matching directory in scratch on HPC
```
# Keep in mind that paths are softlinked on HPC, so it's not straight forward to use relative paths.
# Get the absolute path for ProjectTemplate.sh and copy the path

readlink -f ProjectTemplate.sh
cd scratch
sh /panfs/pfs.local/work/sjmac/e284e911/GeneralScripts/ProjectTemplate.sh BLAST_PRACTICE
```
* Finally make the connection between the two project directories through Atom:
* Open a new window
* Open a new project and make a new config file with a different path.




