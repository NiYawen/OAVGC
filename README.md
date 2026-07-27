## OAVGC - Oral and Airway Viral Genome Catalogue

The Oral and Airway Viral Genome Catalogue (OAVGC) comprises 141,459  high-quality viral genomes reconstructed from 19,997 publicly available and 2,673 newly sequenced human oral and airway metagenomes. The OAVGC contains 68,708 viral operational taxonomic units (vOTUs) and 769 families, predominantly bacteriophages and few eukaryotic viruses. The pipeline below handles raw metagenomic data processing, viral discovery, taxonomic classification, host prediction, and functional annotation.
<br>

------

## Figure 1-7

Folders Figure1–Figure7 contain the codes used to generate the main figures and the associated input files.

## [Scripts](scripts)

The Scripts folder includes shell scripts for OAVGC construction.

### 1. Software Versions & Dependencies

The following software versions were used in the development and execution of the OAVGC pipeline:

| Category                       | Software          | Version            |
| :----------------------------- | :---------------- | :----------------- |
| **Quality Control & Assembly** | Fastp             | v0.23.4            |
|                                | BBTools (bbduk)   | v39.00             |
|                                | Bowtie2           | v2.4.4             |
|                                | Megahit           | v1.2.9             |
| **Viral Identification**       | CheckV            | v0.7.0             |
|                                | DeepVirFinder     | v1.0               |
|                                | VIBRANT           | v1.2.1             |
|                                | geNomad           | v1.7.4             |
|                                | hmmsearch (HMMER) | v3.3.2             |
| **Taxonomy**                   | BLASTN            | v2.12+             |
|                                | vConTACT3         | v3.1.6             |
|                                | ViPTreeGen        | v1.1.2             |
|                                | iTOL              | v7 (Web Interface) |
| **Host & Lifestyle**           | MinCED            | v0.4.2             |
|                                | IPEV              | v4                 |
|                                | BACPHLIP          | v0.9.3             |
| **Annotation**                 | Prodigal-gv       | v2.11.0            |
|                                | eggNOG-mapper     | v2.1.12            |
|                                | DIAMOND           | v2.0.13            |
| **Profiling**                  | Kraken 2          | v2.1.3             |
|                                | Bracken           | v2.8               |

### 2. Script Overview

Descriptions of scripts are shown below:

| Script Name                  | Function Description                                         |
| :--------------------------- | :----------------------------------------------------------- |
| `flow_fastp_bbduk_rmhost.sh` | **Pre-processing:** Quality control (`fastp`), low-complexity filtering (`bbduk`), and human/phiX174 removal (`Bowtie2`). |
| `flow_megahit.sh`            | **Assembly:** De novo assembly with optimized k-list (21, 41, 61, 81, 101, 121, 141). |
| `flow.virus.find.sh`         | **Discovery:** Integrates`flow.virus.{dvf,ckv,vib,busco,geNomad}.sh`  for viral identification |
| `flow.virus.dvf.sh`          | **Identification:** DeepVirFinder scoring (Threshold: Score >0.90, P <0.01). |
| `flow.virus.vib.sh`          | **Identification:** Automated viral identification using `VIBRANT`. |
| `flow.virus.genomad.sh`      | **Decontamination:** Removing non-viral elements (e.g., plasmids) via `geNomad`. |
| `flow.virus.busco.sh`        | **Purity Check:** Bacterial marker filtering using `hmmsearch` (Threshold: BUSCO ratio <5%). |
| `flow.virus.ckv.sh`          | **Quality Control:** Completeness estimation and final viral genome selection (`CheckV`). |
| `flow.virus.clu.sh`          | **Clustering:** vOTU generation (95% nucleotide identity across 85% of the genome length) using in-house `BLASTN` logic. |
| `flow.virus.host.sh`         | **Host Prediction:** Combined CRISPR-spacer matching (`MinCED`) and sequence similarity. |

## 3. Data Availability

The viral genomes data and associated databases for OAVGC are stored at **Zenodo**:[https://doi.org/10.5281/zenodo.18896747](https://doi.org/10.5281/zenodo.18896747).

| Description                                          | Size     | Filename          |
| :--------------------------------------------------- | :------- | :---------------- |
| **Kraken & Bracken Database**                        | 42.5 GB  | `OAVGC.v1.tar.gz` |
| **All Viral Genomes (high-quality, n=141,459)**      | 2.2 GB   | `virus.fna.gz`    |
| **Viral Genome Quality Assessment (CheckV results)** | 2.2 MB   | `virus.ckv.gz`    |
| **vOTU Representative Genomes (n=68,708)**           | 1.0 GB   | `votu.fa.gz`      |
| **vOTU Protein Sequences**                           | 764 MB   | `votu.faa.gz`     |
| **vOTU Coding Sequences (CDS)**                      | 1.1 GB   | `votu.ffn.gz`     |
| **vOTU Gene Annotations (GFF format)**               | 175.6 MB | `votu.gff.gz`     |
| **vOTU Clustering Information**                      | 4.9 MB   | `votu.clu.tsv`    |

## 4. Kraken2 & Bracken Database Structure

The `OAVGC.v1.tar.gz` archive contains a pre-built library for taxonomical profiling. After extraction, the directory structure is as follows:

```text  
OAVGC/  
├── hash.k2d                 # Kraken2 hash table  
├── opts.k2d                 # Kraken2 options  
├── taxo.k2d                 # Kraken2 taxonomy information  
├── seqid2taxid.map          # Mapping file for sequences to TaxIDs  
├── database.kraken          # Intermediate kraken DB file  
├── database150mers.kraken   # Bracken 150bp k-mer counts  
├── database150mers.kmer_distrib # Bracken k-mer distribution file  
└── taxonomy/                # NCBITaxonomy-style directory  
    ├── db.accession2taxid  
    ├── names.dmp  
    ├── nodes.dmp  
    └── prelim_map.txt  
```

## 5. Implementation Summary

### Data Pre-processing

Raw reads were processed with `fastp` (parameters adjusted based on read length ≤100bp or >100bp). Human contaminants were filtered against the **CHM13v2.0** genome.

### Viral Selection Strategy

Only contigs ≥5 kbp were retained for viral identification. A consensus approach was adopted: contigs were retained if identified by `CheckV`, `DeepVirFinder`, or `VIBRANT`, followed by strict removal of plasmid-like (`geNomad`) and bacterial-like (`BUSCO`) sequences.

### Functional Annotation

Protein-coding genes were predicted using `Prodigal-gv`. Functional assignments were performed using `DIAMOND` (against KEGG, CAZy, CARD, and VFDB) and `eggNOG-mapper` (for evolutionary orthologs).

### Profiling

A custom `Kraken2` library was constructed, incorporating OAVGC vOTUs, OAPGC prokaryotic species, and reference fungi. `Bracken` was used for relative abundance estimation and normalization.

