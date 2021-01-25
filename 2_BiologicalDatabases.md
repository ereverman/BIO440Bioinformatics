# Biological Databases:

* First database: in 1965 Margaret Dayhoff at the National Biomedical Research Foundation collected all known protein sequences (65) and published a book 
* Dayhoff's efforts to compile and organize the 65 protein sequences represents one of the first examples of bioinformatics
* Also the one of the first examples of the large impact that computing can have in biology and chemistry
* She also made important advances in biological computing 
  * Interested in comparing protein sequences and examining evolutionary relationships between proteins to understand what proteins did
* Modern datapbases and repositories are built on her first public database



## Timeline:

* Starting off with the Atlas of Protein Sequence and Structure, databases were in print format
* Accumulation of data made this unfeasible
* Computing and the Internet made online repositories possible
* Started off as places to deposit data for others to access, but tools for performing comparisons between data was not far behind
* This timeline shows the 4 major repositories
* There are of course many more.

* These databases aren't fully distinct; they circulate data so that data deposited in one database is available at all of the other ones within 24 hours
* Finding a way to organize data so that it was consistent between the databases was necessary

## Nucleotide Databases:
* Sequence information is held in a flatfile
* Several formats have been used in the past
* The most common and most widely used is the FASTA file

## FASTA format:
* Simplest sequence text file
* > indicates the start of a new sequence entry (there can be many in a file)
* Sample information (can be brief, can contain more information) follows the >
* FASTA format is used for both protein sequences and dna sequences


## Flatfile formats:
* Generally consist of three components
* Header, feature table, sequence

## Go to Euopean Nucleotide Archive
* ebi.ac.uk
* U54469
* Top Selection
* EMBL
* First line:
  * Accession number; sequence version(SV); topology of the DNA molecule (linear); molecule type; ENA class (std = standard annotated and assembled sequence, a way to group sequences together for search queries); taxonomic division (invertebrate); how long the sequence is

* Next section has date information for when the entry was created or last updated.
* Definition lines, describes the type of information contained in the record
* Full taxonomic information for the sequence of interest
* Reference blocks pertain to the publication associated with the sequence data and with the reference for the sequence itself
* ENA-specific cross reference lines

* Feature table: consists of key words that describe the biological property of the sequence
  * Source: range of positions covered by sequence
  * gene: gene sequence starts at postition 80
  * mRNA: there are two transcripts for this gene, join indicates where the exons start and end
  * CDS: sections that ultimately code for the protein. protein id gives the accession number for the protein sequence

* Finally the sequence with some stats on how many of each nucleotide is present
  * 60 bases per row 


# Monday, Jan 25:

## Go to Uniprot
* uniprot.org
* combines resources from Swiss-Prot, TrEMBL, and the Protein Information Reousrce Protein Sequence Database
* three databases:
  * uniprot archive: all publicly available protein sequences compiled from source database
  * uniprotkb: combines entries from swiss prot and trembl
  * UniRef: combines data from uniparc and uniprotkb clusters of proteins based on sequence identity

* UniProtKB:
  * search p09651: human heterogeneous nuclear ribonucleoprotein A1
  * select/deselect areas of interest in index
  * Main window shos basic identifying information about the sequence and whether the sequence has been manually reviewed (STATUS-protein exists)
  * next section contains functional information, gene ontology,links to enzyme and pathway databases
  
* Click subcellular location: takes you down the page to a schematic that shows where the protein is localized
* Find information on disease association and other things...

* Click on Feature Viewer:
  * Summarizes annotations for this sequence. Post translational modification section is expanded
  * Can see protein sequence variants that have been documented so far, and can click on individual ones to get a snapshot of what they do.
  * Variant 351 is assocaited with multiple sclerosis


# Removed from lecture:
## Reference Genome databases:
* ensembl.org
* scroll down to see the "sister sites" with data for other species (non-vertebrate)
* Downloads
* FTP download (File Transfer Protocol)
* Opens a list of organisms with multiple genome and other files associated with each

## FlyBase:
* flybase.org
* Downloads --> current release --> current ftp repository 
* Navigate to the file you need and copy the link for a wget protocol (or curl)




## Online/GUI Bioinformatics Tools:
* Databases often have web-based tools that can perform small-scale bioinformatics analyses
* A slew of other online tools have been developed to fill the needs:
  * to provide users a quick way to access and analyze deposited data
  * to support biologists with diverse levels of computing experience
  * Many are free to use and have a simple graphical user interface
  * Other tools have additional capabilities when used in the command line.
* Typically, online tools are built using the tools we will learn to use in the command line

## BLAST:
https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/#return-footnote-34-1
* Basic Local Alignment Search Tool
* typically uses fasta formatted files
* Looks for matching sequence regions between the query and subject sets:

## Running BLAST:
* BLASTN (Nucleotide BLAST): compares nucleotide query with a database of nucleotide sequences
* Define the query sequence:
  * paste in raw sequence
  * paste in a fasta file
  * enter an NCBI identifier (preferred method to keep track of searches)
* Search AF287139
* Click off align two or more sequences
* Choose the database: Standard database (NR indicates nonredundant, although a little bit of a misnomer because the database is no longer screened against redundant sequences)
* Choose blastn
* BLAST pulls the sequence from NCBI databases to use as the Query sequence
* Next, decide on parameters for your search:
  * Number of results (matches) to output
  * Automatic adjustment for short sequence queries
  * Expected number of change matches. Lower thresholds are more stringent
  * Seed word size (can go as low as 7)
  * Scoring parameters:
    * Reward and Penalty for nucleotide programs like BLASTN
    * reward is for a match, penalty is for a mismatch
    * Can be adjusted based on the expected level of sequence similarity:
    * If sequences are highly conserved may choose to have a higher penalty for mismatches (1/-3)
    * Gap costs: a part of the megablast tool and functions similar to the match/mismatch scores but for allowing gaps in aligned sequences
 * Filters and Masking:
   * Low-complexity regions: masks off segments of the query sequence that are exceptionally common motifs for reducing the number of hits for biologically uninteresting similarities during seeding
   * Species-specific repeats: Masks repeat regions, especially useful if you are working with a species that has a notorious number of repeat regions
   * Mask for lookup table only: includes the "biologically unintersting" regions for the extension phase of the sequence similarity search
   * Mask lowercase letters: FASTA sequences may have lowercase letters which indicate repeat regions such as transposons
  
## BLAST Results:
* Header:
  * RID: can be used to retrieve the results for up to 24 hours
  * Program: Version of BLAST used
  * descriptors of the query sequence
* Descriptions: table of each of the aligned sequences and the score of the alignment (default sorting is E value)
  * Many of the sequences are "predicted"-- means that there isn't any experimental evidence that the sequence produces a protein. this can happen when the genome of a species is sequenced, but gene annotation is only based on computational analyses.
* Graphic Summary: only shows alignments for 100 sequences (there may be more). color coded based on the score of the sequence alignment. can click on the alignment to see the match/mismatch between your query and subject comparison.



# So what does bioinformatics gain us?
* Depending on your source of data (what species, are you generating it yourself or using already deposited data) your interactions with databases and online tools will vary
* Use of tools like BLAST can be great for starting off a new project, getting ideas for question generation and for analyzing results after you have completed your data clean up and larger bioinformatic analyses
  * when you have identified genes that have a snp or that are highly expressed
  * BLAST and GO analyses can be useful in understanding what your data mean in a biological context rather than a computational context.
