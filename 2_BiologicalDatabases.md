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
