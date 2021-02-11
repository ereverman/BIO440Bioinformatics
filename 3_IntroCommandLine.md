## Intro to Command Line for Genomics

## What is a shell:
* Using the Shell and terminal program, we move from a GUI (graphical user interface) to CLI (command line interface)
   * Advantages: Allows repetetive tasks to be automated and reduces errors, can be leveraged with scripts to make tasks much faster, can interact with supercomputers and high performance cluster systems like HPC at KU.
   * Disadvantages: Can be a little opaque, especially if you are only used to working with the GUI.
 
 * The Shell:
   * What is it:
      * A program with the primary purpose of reading commands and running other programs
      * There are different shells, but Bash is the default on most systems and is the shell we will be working with.
      * Learning to use the shell is like learning to speak a new language.
      * Interacting with files and running code requires an understanding of the shell grammar and proper words that designate functions.

## Navigating Files and Directories:

Refer to the powerpoint for an introduction to files and directories:
 * Let's take a look at the shell.
 * $ is the prompt symbol. When you see this, the shell is waiting for your input.
 * I asked you to try an figure out where in your terminal opened. There are 2 ways to figure this out:
  * Look at the information immediately preceeding prompt
  * Use a function:
 ```
 pwd
 # print working directory.
 ```
* On your map, you should move your marker to the "ME" spot. pwd is an especially useful tool. More often than not, when a series of tasks errors out, it's because either your position in your directory structure or your file directories are not referenced correctly.

* Now we know where we are, we can look at the contents of the User directory:
```
ls
# lists the contents of directory
```

* Most functions have additional options that make them more user friendly. For example, if we wanted the directory contents organized line by line:
```
ls -l
# lists the contents of directory line by line with a lot more information.

ls -F
# add a marker to file and directory names to indicate what they are
# / indicates a directory
# @ indicates a link
# * indicates executable
# Depending on your shell preferences, the shell might use colors to indicate whether each entry is a file or directory
```


* One of the big mysteries when I started to learn how to use the shell was knowing what tools were available.
* Here is a quick reference that has useful basic functions:
https://github.com/ereverman/July-August_SWCarpentry_Day1/blob/master/command-handout.md
* Here is an example of how to find the options for ls:
http://linuxcommand.org/lc3_man_pages/ls1.html

```
man ls
# to exit, q

ls --help
```

* What happens if you mistype:
```
ks
# Error message: command not found

ls-F
# spaces are important

ls -s
ls -S
# capitalization is important

ls -lh
# can combine flags

```

* Some suggestions for staying motivated and sane while learning to code:
 * Error messages are actually your friend. Get in the habit of looking at them closely. Sometimes, the problem is easily fixed (typo, path issue).
 * Everyone gets error messages, and code doesn't always run right the first time. More often than not, your error isn't new and someone has provided a solution online.
 

## Navigating Files and Directories:

* Let's get some data that we will use to practice on.
https://github.com/kulibraries/ku-swc-files.git

If need to install curl on gitbash:
https://stackoverflow.com/questions/2710748/run-curl-commands-from-windows-console/20598096

```
cd 

curl -O http://swcarpentry.github.io/shell-novice/data/data-shell.zip

# unzip the directory
unzip data-shell.zip

tree data-shell
# Windows: https://superuser.com/questions/531592/how-do-i-add-the-tree-command-to-git-bash-on-windows
tree.com //a
```


* This repository will appear in your downloads directory.
* Open finder/folder and navigate to the downloaded directory.
* Copy or move the directory to your Desktop.

* The ls command can be used to look at the contents of other directories, not just the one you are in

```
# we just did this:
ls

# now look in one of the directories that is listed:
ls Desktop

## EXERCISE:
# list the contents of your downloaded directory (you should see your downloaded directory).
ls Desktop/ku-swc-files-master/

```

* You can continue to look into the contents of other directories as long as you follow the path through your directories.
* Let's update our map: (powerpoint)
* Your marker should still be on the "Me" position because we are still in our designated user directory but we are about to change that.

```
# change directory
cd Desktop
pwd

# to help you keep track of where you are move your marker to the desktop directory.

# change directory to the data-shell directory:
ls # what is in front of us in the tree?
# two ways to accomplish this because our target directory is two steps down the tree.

cd ku-swc-files-master
cd data-shell

cd ku-swc-files-master/data-shell

```

* How do we return to where we came from?

```
# cd can go forward and backward
cd ku-swc-files-master

# if you want to step out one level use a general shortcut:
cd ../
# this will move you to the parent directory
pwd 
# move your marker.
```
* the ../ is present in every directory by default and we can see it with a special flag:
```
ls -F -a

# -a stands for show all. Normally files and directories that begin with . are hidden.

ls -Fa
# flags can be combined.
```
* There is another directory ./ which means the current directory.
* This is a useful shortcut for moving files around.

* cd can also be used to jump across all of your directories to the home directory
```
cd
pwd

## Exercise:
cd to the data-shell directory with one line

cd Desktop/ku-swc-files-master/data-shell
```

### Relative and Absolute paths (ppt)
* Relative paths point to directories starting from your current position.
 * pros: shorter path (less typing), easier to remember what's immediately above and below your current position.
 * cons: writing scripts that use input files and generate output files are fragile with relative paths. You have to be in a specific place for the scripts to run properly.
* Absoulute paths point to directories starting from the home directory.
 * pros: makes scripts more stable if they can find input files and output directories from anywhere in your filing system.
 * cons: they can get very long, especially if you are referencing files that are buried deep in your system. (but there are ways around this)

* Two exercises.


### Nelle Case Study: File organization tips

* Nelle is a marine biologist, studying marine life from the North Pacific Gyre where trash accumulates rapidly. She studies the effects of this toxic dump on proteins in gelatinous marine organisms.
* Her goal is to analyze 1520 samples, and she needs to learn how to automate the task with the shell in order to not spend 12 hours analyzing each file independently.
* Our goal is to help her put together a pipeline that does this.
* Let's look at her data organization:
 
```
cd north-pacific-gyre/2012-07-03
ls

# data are organized according to year-month-day so that they will be in chronological order. If she used month names, the months would be out of order.
# her names are long but informative. 
```
* Descriptive directory and file names are good practice, but they can be annoying to type out if you are moving around a lot.
* Use tab completion to do two things:
 1. Reduce keystrokes/time.
 2. Double check your path.
  * If the directory doesn't exist or you forgot a step in the path, tab completion won't work.
  
```
cd
# go home

cd Desktop/ku-swc-files-master/data-shell/north...

# pressing tab twice will show you what is next in the path.
```

## Working with files and directories:
* Starting with a template for a project gives you some consistent directories for elements that all projects have (like data, figures, output files)
* Common to need to add directories in the middle of a project that are unique.
 * Maybe the presentation of data will be different for different audiences.
* We are going to do three things:
 1. create a directory hierarchy that matches a desired tree structure.
 2. create files in the hierarchy using an editor to copy and rename existing files.
 3. Delete, copy and move specified files and directories

```
pwd
# move your marker to where you are currently

# cd to the data-shell directory and look at what is inside with ls -F (ppt)
```

* Let's add a new empty directory called thesis

```
mkdir thesis
# no output, so need to check if it did anything:
ls -F
# you can also look at the GUI to see if the directory appears. Directories can be made here as well.

ls -F thesis
# nothing is inside yet
```

* Let's create a text file.
```
cd thesis
nano draft.txt
```

* nano is a text editor. 
* Only works with text, not tables or images.
* nano is one of the least complex editors, others prefer different programs for other features.
 * Less complex is a good place to start because most of us don't actually write theses or manuscripts in the command line...
 * Practical uses for nano include making small edits to scripts to fix an error, updating a ReadMe file that you're using as a log, etc.

```
# type a few lines:
nano is a text editor program.
Use the arrow keys to navigate between letters, pointing and clicking won't help you here.
To save, press CTRL+o and Return/Enter to save the file edits to draft.txt
To exit, press CTRL+x

```
* nano automatically creates a file named draft.txt
* We can also use touch to make an empty file 

```
touch my_file.txt
ls -l
# what happened?
```
* File extensions (ppt)

### Moving files and directories:

* Move to your data-shell directory:
```
cd ../data-shell
```
* We made a draft.txt file in the thesis directory that doesn't have a very informative name.
* We can change the name of a file with move.
* The syntax for mv is similar to others that we will use. 
 * command what-we're-moving where-it-goes 
```
mv thesis/draft.txt thesis/quotes.txt
ls thesis
```
* mv is sometimes dangerous. It won't warn you if a file exists with that name, it will overwrite it.
* mv has flags that can help
```
mv -i thesis/my_file.txt thesis/quotes.txt
```

* mv also works on directories
* Now mv the quotes.txt to the current working directory (data-shell)
```
mv thesis/quotes.txt .
# the second argument can be quick-referenced with the . or you can type the absolute path
# you must specify one of the two ways.

ls thesis
ls quotes.txt
```
* We can also move multiple files at once.
```
cd thesis
touch another_file.txt
ls

cd ../
mv thesis/my_file.txt thesis/another_file.txt .
```
* cp works like mv, except it makes a copy of the file instead of moving it.
```
cp quotes.txt thesis/quotations.txt
ls quotes.txt thesis/quotations.txt
```
* cp also works on directories, but you need an additional flag.
```
cp -r thesis thesis_backup
ls
```
* -r stands for recursive and indicates that the directory and everything it contains should be acted on.
```
ls thesis thesis_backup
```

* removing data files and directories is also possible but should be done with extreme caution. 
* rm does not send objects to the trash bin, it deletes them completely and there's no getting them back.

```
rm -i thesis_backup/quotations.txt
rm thesis # error message
rm -r thesis
```

* Move to the data/ and take a quick break
```
ls -lh
# follow the exercise instructions on the powerpoint

cd data/
mkdir backup/
cp amino-acids.txt animals.txt backup/
ls backup/
ls

```

### How to move multiple files at once:
* Often, we need to move similar output files to a new subdirectory to help with organization.
* We can use wildcards to access multiple files at once.
* There are many types of wildcards that are more or less specific to types of characters.
* * is a wild card that matches zero or more characters.

```
pwd
# move to the molecules/ refer to your map if you need to.

cd data-shell/molecules/
ls
ls *.pdb
ls p*.pdb

# ? is a wildcard that matches exactly one character.
ls ?ethane.pdb
ls *ethane.pdb

```
* powerpoint exercise and/or take a break:
```
pwd
ls
touch fructose.dat
touch sucrose.dat
ls
mkdir analyzed
mkdir raw

mv *.dat analyzed/
cp *.pdb raw
ls
ls analyzed/
ls raw/
```

### Pipes and Filters:
* We've been using commands individually to do single tasks.
* Commands can be daisy chained together to save time and lines of code.
```
ls molecules/
```
* wc can be used to get parameters for a file.
```
wc cubane.pdb
# lines, words, characters

wc *.pdb
# Also shows the total of all lines 
```
* wc has flags like other commands we've looked at so far
```
# we want the number of lines in a file
man wc
q

wc -l *.pdb
```
* What happens if you forget to finish a command?
```
wc -l
CTRL+c to escape
```

* We can combine wc with another command to filter files based on a criterion.
* First we will break it into steps:
```
# which files has the fewest lines:
wc -l *.pdb > lengths.txt
# makes a new file:
ls -lh

# What is inside lengths.txt:
cat lengths.txt
# This prints the entire contents to your terminal
# Useful if the file is small, problematic if the file is very large.
less lengths.txt
q
# also works and allows you to scroll page by page.
```
* sort command will sort the contents of a file.
```
man sort
# lots of options, but we want to sort numerically

sort -n lengths.txt
# doesn't change the file, but outputs the contents to terminal.

sort -n lenths.txt > sorted-lengths.txt
head -n 1 sorted-lengths.txt
# we only want the first number of lines. -n flag means different things for different commands.
tail -n 2 sorted-lengths.txt
# works in a similar way
```
* The > operator will create a new file that doesn't exist or it will completely overwrite an existing file.
* Sometimes, we want to append a new line to an existing file.
* Use the >> operator to append:
```
echo hello > testfile01.txt
cat testfile01.txt

echo "hello again" >> testfile01.txt
cat testfile01.txt
```

* Pipes allow us to link together commands to avoid generating intermediate files
```
ls
# we have lengths.txt and sorted-lengths.txt but we really only wanted sorted-lengths.txt

wc -l *.pdb | sort -n
wc -l *.pdb | sort -n | head -n 1

# common use: copying files over, make sure you copied all of them:

ls -l *.pdb | wc -l

```
* Let's work through an example:
```
cd data/
ls
cat animals.txt | head -n 5 | tail -n 3 | sort -r > final.txt
# takes a subset of the data using positional information from head and tail, then sort those remaining lines in reverse order.

```
* We can use cut to obtain sections of a file:
* cut expects files to be separated by a delimiter (e.g. tab separated)

```
cut -d , -f 2 animals.txt
# the file is comma separated, and we want the second column of the animals.txt file.
```
* uniq is another filter like sort that can be used to reduce the complexity of a file.
* We want to know the types of animals that appear in the file:
```
cut -d , -f 2 animals.txt | sort | uniq

cut -d , -f 2 animals.txt | sort | uniq -c
```

* Let's use Nelle's files as a way to demonstrate how these tools can be used:
* Go to the north-pacific-gyre/
```
cd ../north-pacific-gyre/
ls
cd 2012-07-03
ls

wc -l *.txt
# The expected file size is 300 lines

wc -l *.txt | sort -n | head -n 5
# one of the files has fewer than 300 lines

# are any files too big?
wc -l *.txt | sort -n | tail -n 5
# looks good, but there is another problem. Z indicates a file with missing information.

# Find other files with missing information:
ls *Z.txt
```

* Wildcards with the * can be narrowed down.
```
# if you knew you needed to look for files that had a unique identifier followed by A or B:

ls *A.txt *B.txt
ls *[AB].txt
```

## Loops

* Loops are a programming construct that allos us to repeat a command of set of commands for each item in a list.
* Increases productivity through automation.
* Reduces the amount of typing and time interacting with files.

* For this exercise, change directory to creatures
* This directory contains genome files. Our example will show how to loop over these three files, but the same idea can be applied to any number of files.
* take a look at the files
```
head -5 basilisk.dat monotaur.dat unicorn.dat

# first three lines: common name, classification, updated date
# Loops require the structure of the files to be identical. The goal is to the same thing to each of them, and differences in file structure can cause breakdown
```

* Say we want the classification for each species (second line of each file)
* One strategy is to call for the first two lines with head and then the last of those two lines with tail.
* we can do this for each file or we can use a loop
```
for filename in basilisk.dat minotaur.dat unicorn.dat
do
head -n 2 $filename | tail -n 1
done

# a few things to notic: the full prompt is gone for the other lines of the loop.
# you can combine the tools we've been using such as pipes and tab completion with loops.
```
* In our loop, we used the variable name filename, but there isn't anything special about this word except that it can be informative.
```
for temperature in basilisk.dat unicorn.dat
do
head...

# this seems like it will give the temperature of the files, but the output is still what we asked for before.
```
* variable names can be designated with $, but curly brackets are also sometimes used for clarity: ${filename}
* this is useful if the variable is actually part of a longer name:
 * For example I have a list of file names that have a common part and a unique part
 * I define the variable using the sample identifier as the "thing" that changes in the loop while the rest of the name stays the same.
 * I could define the whole file name as the variable, but at other parts of my pipeline, the file extension changes but the sample id stays the same.
 * Using the unique part of the name makes the variable reusable in other parts of my pipeline.
 
* Let's practice some more:
* move to the molecules directory, refer to powerpoint

```
ls

for datafile in *.pdb
do
ls *.pdb
done


# then try
for datafile in *.pdb
do
ls $datafile
done

# the difference is using a variable name operates on one of the files at a time. 
# in the first case because there are six files that end in .pdb, the loop prints all of the file names six times, once for each position in the list of file names.

# wildcards can also be used to limit the files that are being operated on:

for filename in c*
do
ls $filename
done
```

* In reality, these loops are more complicated than needed because they are redundant with the ls command.
* Now let's practice with some other examples

* Save a file in a loop:
```
for alkanes in *.pdb
do
echo $alkanes
cat $alkanes > alkanes.pdb
done

# look at the file.
# what would we need to change to have all of the contents?

for datafile in *.pdb
do
cat $datafle >> all.pdb
done
```

* questions/break?
* let's go back to the creatures directory.
* Goal is to get lines 81-100 of each file.
```
for filename in *.dat
do
echo $filename
head -n 100 $filename | tail -n 20
done

# echo is an important command because we want some kind of separator between the outputs so we know what file contents is which.
```
* Say we want to make a backup copy of our files called original-unicorn.dat
```
for filename in *.dat
do
cp $filename original-$filename
done

# refer to powerpoint

```

* Nelle has a statistical program she wants to use a loop to execute
* Go to the 2012-07-03 directory in north-pacific-gyre
* the program is called goostats
* goostats takes two arguments: the input file and the output file

```
# first, set up getting the input files:
for datafile in NENE*[AB].txt
do
echo $datafile
done

# next, set up for the output files:
for datafile in NENE*[AB].txt
do
echo $datafile stats-$datafile
done

# final command:
for datafile in NENE*[AB].txt
do echo $datafile
bash goostats $datafile stats-$datafile
done

# If you have a long list of data files, you can open another terminal window while the current process is running to check that the proper outputs are being made
```

* echo is a good way to put "training wheels" on your code. It's possible to wreak a lot of havoc if you have an error and it overwrites files incorrectly.
```
# go to molecules

# this is the correct way to use echo
for datafile in *.pdb
do
echo "cat $datafile >> all.pdb"
done

# This doesn't work as well:
for datafile in *.pdb
do
echo cat $datafile >> all.pdb
done
```

* Loops can be nested to increase their power
```
for species in cubane ethane methane
do
for temperature in 25 30 37 40
do
mkdir $species-$temperature
done
done
```
* One final tool that you can use to look at all the commands you've run is to use history:
```
history 

history | tail -n 25

history | grep bash # for a keyword search.
```
