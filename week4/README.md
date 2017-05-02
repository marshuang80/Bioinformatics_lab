# Week 4. Tardigrades: from gene stealers to space marines (05/01/17)

## Goal: Analyze tardigrade genome and understand how they are able to repair damaged DNA so effectively

## Download the assembled DNA and longest contig
```
wget http://kumamushi.org/data/YOKOZUNA-1.scaffolds.fa
wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/scaffold015.fa
```  

## Download the RepeatMaster and it's library
```
wget http://www.repeatmasker.org/RepeatMasker-open-4-0-7.tar.gz
wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/RepBaseRepeatMaskerEdition-20170127.tar.gz
```

## Download Gene-Mark-ES:
```
wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/GeneMark/gm_et_linux_64.tar.gz
wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/GeneMark/gm_key_64.gz
```

## Install Augustus: 
```
brew install Augustus
```

## Run gm
```
sudo cpan install YAML

```

## Run Augustus
```
augustus --species=human .fa > augustus.abinitio.human.gff 

with human 302
with acyrthosiphon 345

```

## getAnnoFasta
```
perl /usr/local/Cellar/augustus/3.2.2_2/libexec/scripts/getAnnoFasta.pl augustus.abinitio.acyrthosiphon.gff

```

## download mass spectromatry peptide data
```
wget http://public.dobzhanskycenter.ru/mrayko/BIMM185/peptides.fa


```

## WoLF Report
```
g702.t1 details extr: 29, plas: 2, lyso: 1
g1285.t1 details extr: 25, plas: 5, mito: 1, lyso: 1
g3428.t1 details mito: 18, cyto: 11, extr: 2, nucl: 1
g3679.t1 details extr: 26, mito: 2, lyso: 2, plas: 1, E.R.: 1
g4106.t1 details E.R.: 14.5, E.R._golg: 9.5, extr: 7, golg: 3.5, lyso: 3, pero: 2, plas: 1, mito: 1
g5237.t1 details plas: 24, mito: 8
g5502.t1 details extr: 31, lyso: 1
g5503.t1 details extr: 29, plas: 1, mito: 1, lyso: 1
g5510.t1 details plas: 23, mito: 7, E.R.: 1, golg: 1
g5616.t1 details extr: 31, mito: 1
g5641.t1 details extr: 31, lyso: 1
g10513.t1 details nucl: 20, cyto_nucl: 14.5, cyto: 7, extr: 3, E.R.: 1, golg: 1
g10514.t1 details nucl: 19, cyto_nucl: 15, cyto: 9, extr: 3, mito: 1
g12510.t1 details plas: 29, cyto: 3
g12562.t1 details extr: 30, lyso: 2
g13530.t1 details extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5
g14472.t1 details nucl: 28, plas: 2, cyto: 1, cysk: 1
g15153.t1 details extr: 32
g15484.t1 details nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1
```
## CBS Report
```
g702.t1               215          0.054  0.876  0.141   S    2
g1285.t1              193          0.290  0.846  0.028   S    3
g3428.t1              171          0.352  0.047  0.706   _    4
g3679.t1              285          0.136  0.862  0.023   S    2
g4106.t1              642          0.018  0.982  0.039   S    1
g5237.t1              228          0.696  0.067  0.309   M    4
g5502.t1              244          0.266  0.879  0.021   S    2
g5503.t1              215          0.267  0.879  0.021   S    2
g5510.t1              208          0.273  0.137  0.604   _    4
g5616.t1              202          0.032  0.895  0.170   S    2
g5641.t1              192          0.056  0.861  0.122   S    2
g10513.t1             344          0.199  0.026  0.870   _    2
g10514.t1             286          0.088  0.044  0.933   _    1
g12510.t1             513          0.124  0.243  0.753   _    3
g12562.t1             200          0.037  0.873  0.184   S    2
g13530.t1             931          0.038  0.918  0.147   S    2
g14472.t1             445          0.113  0.049  0.910   _    2
g15153.t1             177          0.051  0.845  0.157   S    2
g15484.t1             903          0.095  0.037  0.934   _    1
```


## Extract nuc seq from WoLF
```
#python code
nuc = [a.split() for a in open("WoLF_result.txt") if len(a.split()) > 2]
nuc_seq = [b[0] for b in nuc if b[2] == 'nucl:']
```
## Result
```
g10513.t1
g10514.t1
g14472.t1
g15484.t1
['g10513.t1', 'g10514.t1', 'g14472.t1', 'g15484.t1']
```

## Extract nuc seq from CBS
```
#python script
raw = [a.replace('\t','').split(' ') for a in open("cbs_result.txt")]
seq = [b[4] for b in raw if b[-5] == '_']
```
## Result (not in secretary pathway or membrane)
['g3428.t1',
'g5510.t1',
'g10513.t1',
'g10514.t1',
'g12510.t1',
'g14472.t1',
'g15484.t1']

## Sequence in both CBS and WoLF
```
#python code
overlaps = [a for a in CBS if a in wolf]

#result
['g10513.t1', 'g10514.t1', 'g14472.t1', 'g15484.t1']
```

## Extract seqeuence for overlap genes
```
#read in data form augustus
data = [a.strip('\n') for a in open("augustus.whole.aa") if a != '\n']
#index sequence names
index = [i for i,x in enumerate(data) if x[0] == ">"]
#set padding
index.append(max(index)+10000)
#extract seqeunces
fasta = {data[a]:''.join(data[a+1:b]) for a, b in zip(index[:-1], index[1:])}

#seqs to get
aa = ['>g10513.t1', '>g10514.t1', '>g14472.t1', '>g15484.t1']

#write to file
with open("aa_sequence.txt","w") as o:
    for a in aa:
    o.write(a + '\n')
    o.write(fasta[a] + '\n')
```

## Blast 
```
ID  Description     Query   Coverage    E value     Idnet   Accession
#Non-red
g14472.t1   Damage suppressor protein   100%    0.0 100%    P0DOW4.1
g15484.t1   Vacuolar protein sorting-associated protein 51 homolog  78%  0.0 45% Q155U0.1

#ALL
g3428.t1    Myosin regulatory light chain  91%     8e-65   57%     Q09510.1
g5510.t1    Kita-kyushu lung cancer antigen 1 homolog   11%   2.4     54%     Q4R717.1

#No sig
g10513.t1
g10514.t1
g12510.t1
```


```
