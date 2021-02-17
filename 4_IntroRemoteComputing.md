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
```
cd work/
mkdir GeneralScripts
```


#### Initial set up for interfacing Atom with our remote computer
* This needs to be done once for each new project you start
* Open Atom and navigate to a single untitled script.
* File --> Add Project Folder --> (Navigate to where you want the matching GeneralScripts directory) New Folder --> "GeneralScripts" --> Open
* File --> Open --> GeneralScripts (Side bar should change to show directory and contents)
* Open a second terminal tab or window (this one is linked to your local computer
