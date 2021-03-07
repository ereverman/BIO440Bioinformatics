## Working with fastq files:

* Download some example data:

```
mkdir fastq-example
cd fastq-example

curl -O https://raw.githubusercontent.com/vsbuffalo/bds-files/master/chapter-10-sequence/untreated1_chr4.fq
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/nuccount.py
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/readfq.py

# Counting lines is not as simple as with fasta files:

grep -c "^@" untreated1_chr4.fq
wc -l untreated1_chr4.fq
grep "^@" untreated1_chr4.fq


```
* BIOAWK-- a useful tool that can help explore basic statistics of your fastq files:
* Provides the following fields by file type (a la https://github.com/vsbuffalo/bioawk-tutorial):
  * bed:
     1:chrom 2:start 3:end 4:name 5:score 6:strand 7:thickstart 8:thickend 9:rgb 10:blockcount 11:blocksizes 12:blockstarts
  * sam:
    1:qname 2:flag 3:rname 4:pos 5:mapq 6:cigar 7:rnext 8:pnext 9:tlen 10:seq 11:qual
  * vcf:
    1:chrom 2:pos 3:id 4:ref 5:alt 6:qual 7:filter 8:info
  * gff:
    1:seqname 2:source 3:feature 4:start 5:end 6:score 7:filter 8:strand 9:group 10:attribute
  * fastx:
	1:name 2:seq 3:qual 4:comment

```
conda activate bioawk_env
bioawk -cfastx 'END{print NR}' untreated1_chr4.fq # Counts the numer of records (NR) and reports at the end of the assement

bioawk -cfastx '{print $name, length($seq)}' untreated1_chr4.fq

# Count number of sequences that are shorter than 80 bp
bioawk -cfastx 'BEGIN{ shorter = 0} {if (length($seq) < 80) shorter += 1} END {print "shorter sequences", shorter}' untreated1_chr4.fq
```

* Count Nucleotides in a fastq file:

```
cat nuccount_edit.py

chmod +x nuccount_edit.py

cat untreated1_chr4.fq | ./nuccount_edit.py

```
