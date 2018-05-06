# compare-sets
The program compares two sets of DNA sequences (provided as FASTA files) and returns the percent of sequences in each set that align to sequences in the other set. It also provides alignment information for each sequence.

## Requirements
The program requires perl, bioperl and MUMmer (http://mummer.sourceforge.net/)

## Usage

compare-sets.pl -1 \<set1-scafs\> -2 \<set2-scafs\> -o \<out-prefix\> [--silent] [-p <% identity threshold>]

Where
- **set1-scafs** is a DNA fasta file (genome, metagenome etc)
- **set2-scafs** is a second DNA fasta file
- **out-prefix** is the prefix for all output files (see below)
- **--silent**: No output at all except for error messages
- **% identity threshold** is the % identity that will be used to determine similarity (0..100, default: 96)

### Output
compare-sets.pl will generate these files:
- **\<out-prefix\>.coors, <out-prefix>.delta** - output for the Mummer programs
- **\<out-prefix\>.report.txt** - alignment report for each sequence. File looks as follows:
```
### GCF_001546235.1_ASM154623v1_genomic.fna ###

NZ_KQ955959.1	152422	90615	59.4
NZ_KQ955951.1	99091	86459	87.2
NZ_KQ955947.1	4015	0	0
NZ_KQ955953.1	2668	0	0
NZ_KQ955960.1	2545	2544	100
:
NZ_KQ955949.1	4452	0	0

1949214/2358436 (82.6%) are covered in GCF_001546235.1_ASM154623v1_genomic.fna

### GCF_002076055.1_Bbif1891B_genomic.fna ###

NZ_NAQG01000030.1	5942	0	0
NZ_NAQG01000014.1	8932	0	0
NZ_NAQG01000012.1	327702	285050	86.9
NZ_NAQG01000029.1	2039	595	29.1
NZ_NAQG01000003.1	8677	8005	92.2
:
NZ_NAQG01000015.1	22174	489	2.2

1977267/2418976 (81.7%) are covered in GCF_002076055.1_Bbif1891B_genomic.fna
```
Fields for each sequence are 
```
<sequence-name> <sequence-length> <Number of bps covered by sequences from the other set> <% of bps covered by sequences from the other file>
```
## Example
Download the example files, change directory to the example directory and use the following command line to generate the example output files:

```
compare-sets.pl -1 GCF_001546235.1_ASM154623v1_genomic.fna -2 GCF_002076055.1_Bbif1891B_genomic.fna -o comparison -p 90
```
