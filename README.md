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

The pipeline is divided into two main steps:  
1. The first step is executed separately for each sequencing run (multiplexed data) because it involves removal of tag-jumps. This step primarily focuses on sequence quality control, ITS extraction, and chimera removal.  
2. The second step combines data from all runs, performs denoising and sequence clustering.  

The following command was used to run the Step-1 of the pipeline for each sequencing run:

```bash
nextflow run vmikk/NextITS -r main \
  -profile singularity,hpc \
  -resume \
  --input    /path/to/sequencing/data/<RunID>/<RunID>.fastq.gz \
  --barcodes /path/to/sequencing/data/<RunID>/<RunID>_tags.fasta \
  --primer_forward GTACACACCGCCCGTCG \
  --primer_reverse CCTSCSCTTANTDATATGC \
  --lima_barcodetype "dual_symmetric" \
  --lima_minscore 80 \
  --its_region    "full" \
  --chimera_db    "/path/to/Eukaryome_1.9.3_241222_FullITS_100-800.udb" \
  --outdir        "/path/to/step1_results/<RunID>" \
  -work-dir       "/path/to/step1_work_dir/<RunID>"
```

For reference-based chimera detection, the [EUKARYOME](https://eukaryome.org/) database was used as the reference. The file is available on the University of Tartu's ownCloud [https://owncloud.ut.ee/owncloud/s/PkasGDNDimNssm3](https://owncloud.ut.ee/owncloud/s/PkasGDNDimNssm3).

After completion of Step-1 tasks for all sequencing runs, the following command was used to run the Step-2 of the pipeline:

```bash
nextflow run vmikk/NextITS -r main \
  --step "Step2" \
  -profile singularity,hpc \
  --ampliconlen_min 250 \
  --preclustering "unoise" \
  --unoise_alpha 6.0 \
  --unoise_minsize 1 \
  --clustering "vsearch" \
  --otu_id 0.98 \
  --merge_replicates false \
  --max_MEEP 0.6 \
  --max_ChimeraScore 0.6 \
  --data_path "/path/to/step1_results" \
  --outdir    "/path/to/step2_results" \
  -work-dir   "/path/to/step2_work_dir" \
  -resume
```

