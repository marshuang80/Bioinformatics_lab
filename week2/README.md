# Week2 Lab notebook:Deep sequencing H3N2 to look for mutation (17/04/17)
## 1 Inspect the data from your roommate

Download My roommate's sequence data

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/003/SRR1705853/SRR1705853.fastq.gz
gunzip -d SRR1705853.fastq.gz
mv SRR1705853.fastq roommate.fastq
```

## 2 Align your roommate's data to the reference sequence

Doanload reference file and Index reference with bwa:

```
efetch -db nucleotide -id KF848938.1 -format fasta > KF848938.1.fasta
bwa index reference.fasta
```

Align your roommate’s viral data to the reference sequence and make an mpileup

```
>> bwa mem reference.fasta ~/UCSD/BIMM185/week2/roommate.fastq | samtools view -S -b | samtools sort > roommate.bam

```

Extract all unmapped reads with samtools and count them, determine % mapped

```
samtools view -f4 roommate.bam | wc -l
cat roommate.fastq | wc -l

``` 
**output**:
```
number unmapped: *160*
number of lines in fastq: *1442244*
number of reader: *1442244/4 = 360561*
% unmapped: *0.000437*
% mapped = *1-0.00437 = 0.9995623*
```

## 3. Look for common variants with VarScan

Sort and index the **bam** file, then make a mpileup

```
samtools sort roommate.bam -o roommate_sorted.bam
samtools index roommate_sorted.bam
```

Make a mpileup

```
samtools mpileup -d 1000000 -f reference.fasta roommate_sorted.bam > roommate.mpileup
```

*Given the length of the reference sequence, the number of reads, and the average length, estimate the desired value of d if we want to keep all possible variants.??????*

Run VarScan on the mpileup

```
VarScan mpileup2snp roommate.mpileup --min-var-freq 0.95 --variants --output-vcf > VarScan.vcf
```

**output:**

```
Only SNPs will be reported
Warning: No p-value threshold provided, so p-values will not be calculated
Min coverage:   8
Min reads2:     2
Min var freq:   0.5
Min avg qual:   15
P-value thresh: 0.01
Reading input from roommate.mpileup
1665 bases in pileup file
6 variant positions (6 SNP, 0 indel)
0 were failed by the strand-filter
6 variant positions reported (6 SNP, 0 indel)
```

Pull out the variants in a convenient format:

```
cat VarScan.vcf | awk 'NR>24 {print $1, $2}'
```

**output**

```
KF848938.1 72
KF848938.1 117
KF848938.1 134
KF848938.1 604
KF848938.1 774
KF848938.1 1260
```

Use online sequence editor called WebDSV *http://www.molbiotools.com/WebDSV/index.html* to find the mutations 

```
A72G ACA>ACG Thr24Thr synonymous
C117T -> GCC>GCT Ala39Ala synonymous
G134A -> AGT>AAT Ser45Asn non-synonymous
A604G -> ATC>GTC Ile202Val non-synonymous
T774C TTT>TTC Phe258Phe synonymous
A1260C CTA>CTC Leu420Leu synonymous 
```

Use VarScan again with 0.001 to find rare variants

```
VarScan mpileup2snp roommate.mpileup --min-var-freq 0.001 --variants --output-vcf > VarScan2.vcf
```

**output**

```
Min coverage:   8
Min reads2:     2
Min var freq:   0.001
Min avg qual:   15
P-value thresh: 0.01
Reading input from roommate.mpileup
1665 bases in pileup file
23 variant positions (23 SNP, 0 indel)
1 were failed by the strand-filter
22 variant positions reported (22 SNP, 0 indel)
```

Find frequency of mutation: 

```
cat VarScan2.vcf | awk 'NR>24 {print $4, $2, $5, $10}'|awk -F '[ :]' '{print $1, $2, $3, $10 }'
```

**output**

```
T 38 C 0.42%
A 72 G 99.94%
C 117 T 99.84%
G 134 A 99.89%
A 216 G 0.18%
T 222 C 0.17%
A 254 G 0.23%
A 276 G 0.24%
T 340 C 0.24%
T 409 C 0.19%
A 604 G 99.86%
G 717 A 0.17%
A 744 G 0.18%
T 774 C 99.93%
T 915 C 0.19%
A 1043 G 0.19%
A 1086 G 0.28%
A 1213 G 0.23%
A 1260 C 99.92%
T 1314 C 0.57%
A 1404 G 0.39%
A 1421 G 0.19%
```

download fastq files for 3 controls:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz 
gzip -d SRR1705858.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz 
gzip -d SRR1705859.fastq.gz
wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz 
gzip -d SRR1705860.fastq.az
```
Calcualte how many reads are in each file

```
wc -l SRR1705858
```

SRR1705858.fastq: 1026344
SRR1705859.fastq: 933308
SRR1705860.fastq: 999856
*TODO:* Calculate a rough estimate of the coverage in your samples. 
reference genome length: 1665
read lenght: 151

```
coverage = (read count * read length ) / total genome size.
genome_length = sum([len(a) for a in open("KF848938.1.fasta.1") if a[0] != '>'])
a = [a for a in open("SRR1705858.fastq")]
read_lenght = len(a[1].strip('\n'))
list(map(lambda x: ((x/4)151)/1665,read_count))
```
Coverage: 
SRR1705858.fastq: 23269.96
SRR1705859.fastq: 21160.58
SRR1705860.fastq: 22669.40

look for number of unmapped hits

```
samtools view -f4 SRR1705858.bam | wc -l
```
number of unmapped: 
SRR1705858: 86
SRR1705859: 76
SRR1705860: 76

Use the bwa mem pipeline from last time to align each control fastq file to the reference (KF848938.1.fasta), convert it to a bam file, and sort it.

```
bwa index KF848938.1.fasta
bwa mem KF848938.1.fasta ~/UCSD/BIMM185/week2/SRR1705858.fastq | samtools view -S -b | samtools sort > SRR1705858.bam
bwa mem KF848938.1.fasta ~/UCSD/BIMM185/week2/SRR1705859.fastq | samtools view -S -b | samtools sort > SRR1705859.bam
>> bwa mem KF848938.1.fasta ~/UCSD/BIMM185/week2/SRR1705860.fastq | samtools view -S -b | samtools sort > SRR1705860.bam
```
Run VarScan with a minimum variant frequency of 0.001 (0.1%) on each of the reference
alignments. Be sure to tell VarScan to only output variants and to format the output in the vcf format.
```
samtools sort SRR1705858.bam -o SRR1705858.bam
samtools index SRR1705858.bam
samtools sort SRR1705859.bam -o SRR1705859.bam
samtools index SRR1705859.bam
samtools sort SRR1705860.bam -o SRR1705860.bam
samtools index SRR1705860.bam
samtools mpileup -d 1000000 -f KF848938.1.fasta SRR1705858.bam > SRR1705858.mpileup
samtools mpileup -d 1000000 -f KF848938.1.fasta SRR1705859.bam > SRR1705859.mpileup
samtools mpileup -d 1000000 -f KF848938.1.fasta SRR1705860.bam > SRR1705860.mpileup
VarScan mpileup2snp SRR1705858.mpileup --min-var-freq 0.001 --variants --output-vcf > SRR1705858.vcf
VarScan mpileup2snp SRR1705859.mpileup --min-var-freq 0.001 --variants --output-vcf > SRR1705859.vcf
VarScan mpileup2snp SRR1705860.mpileup --min-var-freq 0.001 --variants --output-vcf > SRR1705860.vcf
```

**SRR1705858:**

Min coverage:   8

Min reads2:     2

Min var freq:   0.001

Min avg qual:   15

P-value thresh: 0.01

Reading input from SRR1705858.mpileup

1665 bases in pileup file

58 variant positions (58 SNP, 0 indel)

1 were failed by the strand-filter

57 variant positions reported (57 SNP, 0 indel)

**SRR17058589:**

Min coverage:   8

Min reads2:     2

Min var freq:   0.001

Min avg qual:   15

P-value thresh: 0.01

Reading input from SRR1705859.mpileup

1665 bases in pileup file

54 variant positions (54 SNP, 0 indel)

2 were failed by the strand-filter

52 variant positions reported (52 SNP, 0 indel)


**SRR1705860:**

Min coverage:   8

Min reads2:     2

Min var freq:   0.001

Min avg qual:   15

P-value thresh: 0.01

Reading input from SRR1705860.mpileup

1665 bases in pileup file

61 variant positions (61 SNP, 0 indel)

0 were failed by the strand-filter

61 variant positions reported (61 SNP, 0 indel)


Use awk to parse each vcf file so that you get three lists containing the reference base, position,alternative base, and frequency. 

```
cat VarScan2.vcf | awk 'NR>24 {print $4, $2, $5, $10}'|awk -F '[ :]' '{print $1, $2, $3, $10 }'
```

find average and std:
```
percentage_58 = np.array([float(a.split()[3].strip("%")) for a in open("SRR1705858.txt")])
percentage_58.mean()
percentage_58.std()
```
SRR1705858

    mean: 0.25649122807017544

    std: 0.071093988397129201

SRR1705859

    mean: 0.2369230769230769

    std: 0.051870343592256042

SRR1705860

    mean: 0.250327868852459

    std: 0.077395454882323786

Roommate

    mean: 27.429090909090913 

    std: 44.377247129693806

find ones in roommate that is outside of 3std 

C 307 T 0.94%

C307T CCG>TCG Pro103Ser nonsynonymous

T 1458 C 0.84% 

T1458C TAT>TAC Tyr486Tyr synonymous
