## Objectives:
* How to learn bioinformatics
* Biological databases

## Write code for humans, write data for computers:
* Code should be:
  * Well-commented
  * Modular
  * Human-readable
  
### Example:
* Bioinformatics projects are complex and have several sources of input and output files
* Requires knowing how your project directory is structured
* Many analyses have an expected worflow, making it possible to anticipate what directories you will need


* Many people use project template scripts to generate a project directory structure prior to starting their analyses:
* 2 examples
* Both shell scripts (.sh)
* Two scripts do the exact same thing
  * Left script has fewer lines
  * Right script has more line
  * Right script has comments (#)
  * Right script has sections (variable assignment and implementation)
  * Left script requires each line be changed in order to use for a new project
  * Right script requires one line to be changed in order to use for a new project
  * Right one took longer to build
  * Left one takes longer to use

## Let your computer do the work for you:
* I just showed you a project template script.
* Instead of making the directories as i need them, i can create the major ones i need with one command

```
cd Desktop
less VariantCallingTemplate1.sh
less VariantCallingTemplate2.sh

sh VariantCallingTemplate1.sh
sh VariantCallingTemplate2.sh

tree VariantCalling1
tree VariantCalling2
```


## Use existing libraries or environments when possible:
* Many packages and tools require dependencies (other packages and tools) in order to run
* Installing tools over and over again can take time, sometimes hours
* Environments with packages and tools already loaded can save time 
  * Often re-running analyses through the course of a project
* Anaconda
  * Provides a way to create a specific environment to run analyses
  * Portable between computers and collaborators
  * Preserved versioning used in an analysis
  
## Treat data as read-only
* Designate a raw directory and keep a backup separate from your active analyses directories
* Change the permissions on data files to prevent them from being overwritten or altered

## Document Code and Projects with a README
* Can take many forms
* Sometimes optimal to have multiple readme files
  * One for your "Scripts" directory to define what each script is used for
  * For final data files that are submitted with manuscripts
  






