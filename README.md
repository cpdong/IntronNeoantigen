# IntronNeoantigen
A pipeline for calling intron retained derived neoantigens using RNA-Seq data

##  Getting Started
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
    - Finish the installation following the [INSTALLATION INSTRUCTIONS](http://www.cbs.dtu.dk/services/doc/netMHCpan-4.0.readme)

IntronNeoantigen requires the following Python modules:
- [HTSeq](https://pypi.org/project/HTSeq/)
- [pysam](https://pypi.org/project/pysam/)
- NumPy
- collections
- multiprocessing
