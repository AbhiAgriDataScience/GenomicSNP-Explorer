# GenomicSNP-Explorer: Variant Analysis in Maize

## Overview
- This project aims to analyze overlaps between genomic features (e.g., SNPs, coding regions, regulatory elements) to gain biological insights. 
- However, due to the limited availability of whole genome sequencing (WGS) data for drought-resistant and susceptible maize genotypes, the primary objective here is to execute the entire workflow rather than derive biological conclusions.
- This project contains a bioinformatics pipeline for analyzing genomic variants in drought-resistant and susceptible maize genotypes.
- The workflow includes SNP detection, genomic feature annotation, and comparative analysis, using publicly available data from NCBI SRA and MaizeGDB.

## Data Acquisition
To obtain SNP data, the project utilized **SNPVersity 2.0**, a web-based tool for visualizing maize diversity ([MaizeGDB](https://wgs.maizegdb.org/)). The dataset used is **MaizeGDB 2024 (High Quality Dataset)**, specifically the project:

- **Project Name**: Deep DNA resequencing of the association mapping panel ([PRJNA531553](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA531553))
- **Institution**: Huazhong Agricultural University
- **Total Accessions**: 519
- **Other Accessions**: 228

### Selected Genotypes:
| Genotype | Drought Resistance | SRR Accession |
|----------|--------------------|--------------|
| YE478    | Resistant          | SRR8906861   |
| H21      | Susceptible        | SRR8906666   |


## Workflow
The pipeline includes the following steps:

### 1. **Data Download and Preprocessing**
   - Retrieve raw sequencing data (FASTQ) from NCBI using SRR numbers.
   - Perform quality control (QC) using **FastQC** and trim reads using **Trimmomatic**.

### 2. **Genome Alignment**
   - Index the maize reference genome.
   - Align reads to the reference genome using **BWA**.
   - Convert and sort SAM to BAM files using **SAMtools**.

### 3. **Variant Calling**
   - Call SNPs and small indels using **BCFtools**.
   - Filter SNPs based on quality metrics.
   - Identify and extract unique variants

### 4. **Genomic Feature Annotation**
   - Annotate SNPs with gene regions using **BCFTools**.
   - Identify overlaps with regulatory elements and coding regions.

### 5. **Comparative Analysis**
   - Compare SNP distributions between resistant and susceptible genotypes.
   - Generate summary statistics and visualizations.

## Repository Structure
```
GenomicSNP-Explorer/
│── data/                              # Raw and processed data
│── results/                           # Output files, figures, and reports
│── GenomicSNP_analysis.ipynb    # Jupyter notebook for interactive analysis
│── GenomicSNP_pipeline.sh             # Bash script for automated analysis
│── README.md                          # Project documentation
│── requirements.txt                   # Dependencies and installation requirements
```

## Running the Workflow
#### 1. Running the Bash Script (Full Pipeline)
To run the entire pipeline non-interactively, execute:

```bash
bash GenomicSNP_pipeline.sh
```
This will perform data processing, alignment, SNP calling, and annotation.

#### 2. Running the Jupyter Notebook (Step-by-Step Analysis)
To run step-by-step analysis in interactive mode, open the notebook:

```bash
jupyter notebook GenomicSNP_analysis.ipynb
```

## Dependencies
To replicate the analysis, install the required dependencies listed in `requirements.txt`:

```bash
conda create --name genomic_analysis_env python=3.9
conda activate genomic_analysis_env
pip install -r requirements.txt
```

### Key Tools and Packages
- **Python:** Pandas, NumPy, Matplotlib, Biopython
- **Bioinformatics Tools:** sra-tools, FastQC, Trimmomatic, BWA, SAMtools, BCFtools, bedtools, vcftools


## Acknowledgments
- This project uses publicly available datasets from MaizeGDB and NCBI SRA. 
- Special thanks to Huazhong Agricultural University for providing high-quality genomic data.
































README.txt
# Project: Genomic Feature Analysis 
Objective: Analyze overlaps between genomic features (eg. SNPs, coding regions, or regulatory elements) to gain biological insights.

I did not find any SNP databases or VCF files of drought-resistance and non-resistant genotypes in maize.
I checked online sources and found that there were few genotypes assoicated with drought resistance (G10, G123, SC302, CE704, KSC260, SC302, DC370, SC647, KSC704, H082183, 478) and susceptibility (2023, G2, Lv28).
However, I did not find whole genome sequencing of many genotypes and 
only found 478 (resistant), Lv28 and H21 (drought-sensitive)

The main objective here is to perform the whole workflow rather than gain biological insights since 3 genotypes wont be enough to infer biological insights. 

For data acquisition, I used SNPVersity 2.0: a web-based tool for visualizing maize diversity (https://wgs.maizegdb.org/). I used MaizeGDB 2024 (High Quality Dataset) and chose the project 'Deep DNA resequencing of the association mapping panel (PRJNA531553) - "Other" accessions
Huazhong Agricultural University (Reference): 519 total accessions, 228 other accessions' as it had SRR numbers of those 3 genotypes.
Drought resistant (YE478: SRR8906861), Drought susceptible genotypes(H21: SRR8906666 and LV28 (SRR8907091) and 