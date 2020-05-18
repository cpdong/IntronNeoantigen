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
- [NetMHCpan-4.0](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.0)
    - complete the installation following the [INSTALLATION INSTRUCTIONS](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme)


IntronNeoantigen requires the following Python modules:
- [HTSeq](https://pypi.org/project/HTSeq/)
- [pysam](https://pypi.org/project/pysam/)
- NumPy
- collections
- multiprocessing

##  Getting Started

In order to run , we need to roll back the reference to an earlier version. First, fetch IMGT/HLA database version 3.24.0:
```
./intronneoantigen index /path/to/gencode.v32.annotation.gtf
```
<strong>1. Genotyping with arcasHLA:</strong>
```
/path/to/arcasHLA extract test/test.bam -o test/output --paired -t 8 -v
/path/to/arcasHLA genotype test/output/test.extracted.1.fq.gz test/output/test.extracted.2.fq.gz -g A,B,C -o test/output -t 8 -v
```
<strong>2. Quantification reads fall on intron region:</strong>
```
./intronneoanitgen call test/test.bam -o test/output --paired -t 8 -v
```
<strong>3. Predicting MHC-I binding affinity with NetMHCpan:</strong>
```
./intronneoanitgen netMHCpan test/output/test.extracted.1.fq.gz test/output/test.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o test/output -t 8 -v
```
