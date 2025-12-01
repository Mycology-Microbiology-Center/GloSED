# Data Reuse Information (DRI) for GloSED

## What is DRI?

Data Reuse Information (DRI) is a machine-readable metadata standard proposed by [Hug et al. (2025; DOI:10.1038/s41564-025-02116-2)](https://www.nature.com/articles/s41564-025-02116-2) for equitable sharing of public microbiome data. The DRI tag contains ORCID identifiers of data creators who prefer to be contacted before data reuse, promoting collaboration and proper attribution while facilitating open science.

## DRI implementation in GloSED

GloSED implements DRI to support equitable data reuse while enabling collaboration between data creators and consumers. The presence of DRI indicates that creators welcome contact and collaboration opportunities.

## Files

### `DRI.json`

Machine-readable DRI metadata for the entire GloSED dataset:
- **ORCID identifiers**: Stable, unique identifiers for data creators  
- **Contact information**: Email addresses for direct communication  
- **Dataset details**: Version and descriptive notes  
- **Usage**: Parse programmatically to check reuse requirements  

### `DRI.csv`

Granular DRI mapping for individual samples and projects:
- **Columns**: accession, biosample, bioproject, dri_orcids  
- **Purpose**: Maps specific data entries to their associated DRI ORCIDs  
- **Usage**: Cross-reference with ENA/PlutoF accession numbers  

## DRI ORCIDs

- **0000-0002-1635-1249**: Leho Tedersoo (Project Lead)  
- **0000-0003-2786-2690**: Vladimir Mikryukov (Bioinformatics Lead)  

## How to use DRI

1. **Check for DRI presence**: If DRI files exist, contact the listed creators before reuse  
2. **Contact data creators**: Use provided ORCIDs or email addresses  
3. **Collaborate**: Discuss potential collaboration opportunities  
4. **Cite properly**: Acknowledge data creators in publications  

## Citation

Hug LA, et al. (2025). A roadmap for equitable reuse of public microbiome data. *Nature Microbiology*. DOI:10.1038/s41564-025-02116-2 URL:[https://www.nature.com/articles/s41564-025-02116-2](https://www.nature.com/articles/s41564-025-02116-2)

