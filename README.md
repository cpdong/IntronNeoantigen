# IntronNeoantigen
A pipeline for calling intron retained derived neoantigens using RNA-Seq data

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

##  Prerequisites
Make sure the following programs are in your `PATH`:
- [Samtools](http://www.htslib.org/)
- [bedtools](http://bedtools.readthedocs.io/)
- [pigz](https://zlib.net/pigz/)
- [Kallisto](https://pachterlab.github.io/kallisto/)
- Python 3.6+
- [arcasHLA](https://github.com/RabadanLab/arcasHLA)
    - remember to fetch IMGT/HLA database before use: arcasHLA reference --update
- [NetMHCpan-4.1](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1)
    - complete the installation following the [INSTALLATION INSTRUCTIONS](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme)


IntronNeoantigen requires the following Python modules:
- [HTSeq](https://pypi.org/project/HTSeq/)
- [pysam](https://pypi.org/project/pysam/)
- [pybedtools](https://daler.github.io/pybedtools/)
- NumPy
- collections
- multiprocessing

##  Getting Started

In order to run , we need to roll back the reference to an earlier version. First, fetch IMGT/HLA database version 3.24.0:


###  Genotyping with arcasHLA: ###
```
/path/to/arcasHLA extract test/test.bam -o test/output --paired -t 8 -v
/path/to/arcasHLA genotype test/output/test.extracted.1.fq.gz test/output/test.extracted.2.fq.gz -g A,B,C -o test/output -t 8 -v
```

### 1. Re-annotation the GTF file: ###
```
./intronneoantigen index -gtf /path/to/gencode.v32.annotation.gtf -out /path/to/index_dir -t 8
```

### 2. Quantification reads fall on intron region: ###
```
./intronneoanitgen calling -b /path/to/sample.bam -g /path/to/index_dir/gencode.v32.annotation.re-annotation.gtf -t 20 -o /path/to/out
```
#### Options: ####
- `-m`          : mode of calling, normal/tumor. If tumor mode choose, a reference of intron retained frequency is required to filter commonly retained introns.
- `-b`          : input file with full path
- `-g`          : re-parsed gtf file from the index step
- `-t`          : number of threads (default: 1)
- `-o`          : output directory (default: `.`)
- `-c`          : minumum number of reads falls on intron region (default: 50)
- `-p`          : pecentage retain index,pri (default: 0.05)
- `-n`          : novel intron retain events minumum number of reads (default: 10)
- `-j`          : anchor length of read fall on the intron region(default:10, optional) 

### 3. Predicting MHC-I binding affinity with NetMHCpan: ###
```
./intronneoanitgen present -genotype /path/to/sample.genotype.json -file /path/to/sample_intron_call.result.txt -fasta /path/to/GRCh38.primary_assembly.genome.fa -gtfindex /path/to/gencode.v32.annotation.re-annotation.gtf -len 8,9,10,11 -thread 20 -outdir /path/to/out
```
#### Options: ####
- `-genotype`      : genotype file generated from arcasHLA(default). If this option choosed, the '- hla' will be invalid.
- `-hla`           : HLA type provide by user, multiple HLA seperate with coma with out space,example: HLA-A01:01,HLA-B02:02
- `-file`          : intron-retain events reported from the step 2
- `-fasta`         : full path of reference genome fasta file(indexed)
- `-gtfindex`      : re-parsed gtf file from the index step
- `-len`           : peptide length for MHC prediction,multiple length(8,9,10,11)separated by commas and without spaces(default: 9)
- `-thread`        : number of threads (default: 8)
- `-outdir`        : output directory
<br>

## Citations ##
Dong C, Liu Y(2020) IntronNeoantigen: identify intron retained derived neoantigens using RNA-Seq data. Prepared for submission.
