## GloSED: Global Standardised Soil Eukaryome Dataset [<img src='assets/MycologyMicrobiologyCenter_logo.png' align="right" height="100" alt="Mycology and Microbiology Center">](https://mmc.ut.ee/?lang=en)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17827890.svg)](https://doi.org/10.5281/zenodo.17827890)

The Global Standardised Soil Eukaryome Dataset (GloSED) is a global metabarcoding dataset encompassing the entire spectrum of soil eukaryotes, covering more than 4100 sampling sites in 121 countries and nearly one million taxonomically annotated OTUs. All samples were collected and analysed using standardised field, laboratory, and bioinformatic protocols, and are accompanied by extensive soil and land-cover metadata. This repository hosts the analysis scripts and supporting files associated with the GloSED data descriptor (*Scientific Data*, in preparation) and the Zenodo data release linked above.

## Repository structure

- `bin/`: Analysis scripts  
- `DRI/`: Data Reuse Information tag used to support equitable data reuse and sequence accession number list  
- `assets/`: Auxiliary files for the repository  

## Bioinformatic processing

GloSED was analysed using the fully automated bioinformatics pipeline [NextITS](https://github.com/vmikk/NextITS) (DOI:10.5281/zenodo.15074882) implemented with the workflow manager [Nextflow](https://www.nextflow.io/). Pipeline dependencies are containerised and available on Docker Hub [https://hub.docker.com/r/vmikk/nextits](https://hub.docker.com/r/vmikk/nextits) and Singularity library [https://cloud.sylabs.io/library/vmiks/nextits/nextits](https://cloud.sylabs.io/library/vmiks/nextits/nextits).

To install Nextflow, use the following command:

```bash
curl -s https://get.nextflow.io | bash
```

To install NextITS, run:

```bash
nextflow pull vmikk/NextITS
```

