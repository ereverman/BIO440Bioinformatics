## Working with fastq files:

* Download some example data:

```
mkdir fastq-example
cd fastq-example

curl -O https://raw.githubusercontent.com/vsbuffalo/bds-files/master/chapter-10-sequence/untreated1_chr4.fq
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/nuccount_edit.py
curl -O https://raw.githubusercontent.com/ereverman/BIO440Bioinformatics/main/readfq_edit.py

# Counting lines is not as simple as with fasta files:

grep -c "^@" untreated1_chr4.fq
wc -l untreated1_chr4.fq
grep "^@" untreated1_chr4.fq

conda activate bioawk_env
bioawk -cfastx 'END{print NR}' untreated1_chr4.fq

# bioawk is a useful tool that can help you explore basic statistics of your fastq files:

bioawk -cfastx '{print $name, length($seq)}' untreated1_chr4.fq

# Count number of sequences that are shorter than 80 bp
bioawk -cfastx 'BEGIN{ shorter = 0} {if (length($seq) < 80) shorter += 1} END {print "shorter sequences", shorter}' untreated1_chr4.fq

```

* Count Nucleotides in a fastq file:

```
cat nuccount.py

chmod +x nuccount.py

cat contam.fastq | ./nuccount.py

```
