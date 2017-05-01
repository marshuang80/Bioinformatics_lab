# Week 3 Lab report: E.coli outbreak investigation (24/04/17)

## Download the squencing datasets: 
```
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R2_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R2_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R2_001.fastq.gz
```

## Run fastQC on files:
```
fastqc -o . ./SRR292678sub_S1_L001_R1_001.fastq ./SRR292678sub_S1_L001_R2_001.fastq
```

## Number of reads:
```
SRR292678: 5499346
SRR292770: 5102041
SRR292862: 5102041

```

## Install and use Jellyfish
```
brew install jellyfish
jellyfish count -m 31 -s 10000000 -o 31 -C SRR292862_S2_L001_R1_001.fastq
jellyfish histo 31 > 31.histo

```

##Make Histogram
```
R

pdf('frag_1_31.pdf')

r31 <- read.table("31.histo")

plot(r31[2:100,], type="l", main="frag_1_31", xlab="number of times kmer appears", ylab="number of k-mers")

dev.off()

quit()
```

## Get SPades
```
brew install SPades
./SPAdes-3.10.1-Darwin/bin/spades.py --test

SPAdes log can be found here: /Users/marshuang/UCSD/BIMM185/week3/spades_test/spades.log

Thank you for using SPAdes!

./SPades-3.10.1-Darwin/bin/spades.py --pe1-1 ./SRR292862_S2_L001_R1_001.fastq --pe1-2 ./SRR292862_S2_L001_R2_001.fastq -o ./

```

## Run SPades on 3 libraries
```
./SPades-3.10.1-Darwin/bin/spades.py --pe1-1 ./SRR292862_S2_L001_R1_001.fastq --pe1-2 ./SRR292862_S2_L001_R2_001.fastq --mp1-1 ./SRR292770_S1_L001_R1_001.fastq --mp1-2 ./SRR292770_S1_L001_R1_001.fastq --mp2-1 ./SRR292678sub_S1_L001_R2_001.fastq --mp2-2 ./SRR292678sub_S1_L001_R2_001.fastq -o ./

```

*SPades REPORT IN FOLDER*

## Install prokka
```

sudo cpan Time::Piece XML::Simple Digest::MD5 Bio::Perl
brew tap homebrew/science

```

## Run Prokka
```

**GET PRAKKA COMMAND**
```


## Prokka result:
```

80    tRNAs
0    rRNAs
1    CRISPRs
5064    CDS
2923    Unique gene codes

```

## RNAmmer
```
6 copies o  16s_rRNA
1529 bp

```
## Blast Result
```
SequenceID: NC_011748.1


```

## Mauve
```
stxB 3483605-3483874
stxA 3483886-3484845

Found Phage proteins around these genes
From 449066 to 487444
```

##Antibiotic resistance detection


## Antibiotic resistance mechanism
```
Two beta-lactamase genes are found:
- bla_1
- bla_2
both of them did not align to any part of the reference
```



